# =============================================================================
# 11_count_distribution.R
# -----------------------------------------------------------------------------
# Which statistical distribution describes a night's catch?
#
# WHY IT IS A CONDITIONAL QUESTION
#   The marginal distribution of nightly counts is not a useful target: the
#   nights differ enormously in what they should hold, from a locality-year
#   totalling a hundred moths to one totalling seventy thousand, and over the
#   season from the first nights of flight to the peak. The question worth
#   asking is conditional -- given what a night is expected to hold, how are
#   the counts scattered around it?
#
#   So every candidate here shares one mean structure, the `reg+` predictor of
#   R/10 (the regional BioSIM curve with its empirical correction, plus a
#   locality-year level), and they differ only in the distribution of the
#   count around that mean.
#
# THE CANDIDATES
#   poisson        variance = mu; the reference, certain to fail
#   nbinom1        variance = mu * (1 + a); overdispersion proportional to mu
#   nbinom2        variance = mu + mu^2 / theta; the usual negative binomial
#   poisson-ln     Poisson with a lognormal nightly multiplier (an
#                  observation-level random effect): the classical model for
#                  counts of individuals arriving at a trap whose catchability
#                  varies from night to night
#   zinb           nbinom2 with a zero-inflation component, for nights that
#                  produce a structural zero rather than a small count
#
# WHAT IS REPORTED
#   AIC and BIC, the fitted dispersion, the observed and expected number of
#   zero-catch nights, and the empirical mean-variance relationship with the
#   slope each family implies (1 for Poisson and nbinom1, 2 for nbinom2).
#
# INPUT   data/processed/nightly_count_model.rds (the nights table of R/10)
# OUTPUT  data/processed/count_distribution.rds  list(fits, table, mv, summary)
#
# REQUIREMENTS  dplyr, splines, glmmTMB
#
# RUNTIME  a few minutes; the zero-inflated and lognormal fits are the slow
#          ones.
#
# USAGE  source("R/11_count_distribution.R")   # from the project root
# =============================================================================

library(dplyr)
library(splines)
library(glmmTMB)


## ---- 1. Configuration -------------------------------------------------------

out_file <- file.path("data", "processed", "count_distribution.rds")

## Degrees of freedom of the correction spline, as in R/10.
df_corr <- 4

reuse <- TRUE


## ---- 2. The nights, with the mean structure of `reg+` -----------------------

nights <- readRDS(file.path("data", "processed", "nightly_count_model.rds"))$nights %>%
  mutate(site_year = factor(site_year),
         obs_id    = factor(row_number()))

form <- moths ~ log_reg + ns(delta_reg, df = df_corr) + (1 | site_year)


## ---- 3. Fit the candidates ---------------------------------------------------
## Each is the same mean structure with a different conditional distribution.
## The locality-year level is a random intercept rather than a fixed effect,
## which keeps the models comparable and the fits quick.

fit_or_null <- function(label, ...) {
  f <- try(glmmTMB(...), silent = TRUE)
  if (inherits(f, "try-error") || !is.finite(stats::AIC(f))) {
    warning("fit failed: ", label, call. = FALSE)
    return(NULL)
  }
  f
}

if (reuse && file.exists(out_file)) {
  previous <- readRDS(out_file)
  fits <- previous$fits
} else {
  fits <- list(
    poisson = fit_or_null("poisson", form, data = nights, family = poisson()),
    nbinom1 = fit_or_null("nbinom1", form, data = nights, family = nbinom1()),
    nbinom2 = fit_or_null("nbinom2", form, data = nights, family = nbinom2()),
    `poisson-ln` = fit_or_null(
      "poisson-ln",
      update(form, . ~ . + (1 | obs_id)), data = nights, family = poisson()),
    zinb = fit_or_null("zinb", form, data = nights, family = nbinom2(),
                       ziformula = ~ 1)
  )
  fits <- fits[!vapply(fits, is.null, logical(1))]
}


## ---- 4. Compare ---------------------------------------------------------------
## Zero-catch nights and the upper tail are checked by simulating from each
## fitted model rather than analytically: the families differ in where their
## random effects sit, so their conditional means are not on equal footing,
## and only a draw from the whole fitted model is comparable across them.

n_sim <- 100

sim_check <- function(fit) {
  s <- try(stats::simulate(fit, nsim = n_sim, seed = 1978), silent = TRUE)
  if (inherits(s, "try-error")) return(c(zeros = NA_real_, max = NA_real_))
  m <- as.matrix(as.data.frame(s))
  c(zeros = mean(colSums(m == 0)), max = stats::median(apply(m, 2, max)))
}

sims <- vapply(fits, sim_check, numeric(2))

comparison <- tibble::tibble(
  family    = names(fits),
  parameters = vapply(fits, function(f) attr(stats::logLik(f), "df"), numeric(1)),
  logLik    = vapply(fits, function(f) as.numeric(stats::logLik(f)), numeric(1)),
  AIC       = vapply(fits, stats::AIC, numeric(1)),
  BIC       = vapply(fits, stats::BIC, numeric(1)),
  dispersion = vapply(fits, function(f) tryCatch(sigma(f), error = function(e) NA_real_),
                      numeric(1)),
  zeros_simulated = sims["zeros", ],
  max_simulated   = sims["max", ]
) %>%
  mutate(dAIC = AIC - min(AIC)) %>%
  arrange(AIC)

zeros_observed <- sum(nights$moths == 0)
max_observed   <- max(nights$moths)

## Variance components of the lognormal fit, if it converged: the nightly
## multiplier is the quantity of interest.
nightly_sd <- tryCatch(
  sqrt(glmmTMB::VarCorr(fits[["poisson-ln"]])$cond$obs_id[1]),
  error = function(e) NA_real_)


## Where each family puts its mass: the share of nights in each count band,
## observed and simulated, which is where the differences between the
## quadratic-variance families show up.

bands <- c(-1, 0, 9, 99, 999, 9999, Inf)
band_labels <- c("0", "1-9", "10-99", "100-999", "1,000-9,999", "10,000+")

band_share <- function(x) {
  as.numeric(table(cut(x, bands, labels = band_labels)) / length(x))
}

band_table <- bind_rows(
  tibble::tibble(source = "observed", band = band_labels,
                 share = band_share(nights$moths)),
  bind_rows(lapply(names(fits), function(nm) {
    s <- try(stats::simulate(fits[[nm]], nsim = n_sim, seed = 1978), silent = TRUE)
    if (inherits(s, "try-error")) return(NULL)
    m <- as.matrix(as.data.frame(s))
    tibble::tibble(source = nm, band = band_labels,
                   share = rowMeans(apply(m, 2, band_share)))
  }))
) %>%
  mutate(band = factor(band, levels = band_labels))


## ---- 5. The empirical mean-variance relationship ------------------------------
## Nights are binned by the fitted mean of the best-fitting model and the
## variance within each bin is compared with it. A slope of 1 on the log scale
## is Poisson-like (variance proportional to the mean), 2 is negative-binomial
## NB2-like (standard deviation proportional to the mean).

## The mean is taken from the nbinom2 fit, whose fitted values condition on
## the locality-year level only. Using the Poisson-lognormal's fitted values
## would condition on its nightly random effect as well -- the very
## dispersion this diagnostic exists to display -- and the points would
## collapse onto the Poisson line.
mu_fit <- if (!is.null(fits[["nbinom2"]])) fits[["nbinom2"]] else fits[[1]]

mv <- nights %>%
  mutate(mu = stats::fitted(mu_fit),
         bin = dplyr::ntile(mu, 20)) %>%
  group_by(bin) %>%
  summarise(n = n(), mean_mu = mean(mu), mean_y = mean(moths),
            var_y = stats::var(moths), .groups = "drop") %>%
  filter(n > 5, var_y > 0)

mv_slope <- stats::coef(stats::lm(log(var_y) ~ log(mean_y), data = mv))[[2]]


## ---- 6. Write and report -------------------------------------------------------

summary_list <- list(
  n_nights       = nrow(nights),
  zeros_observed = zeros_observed,
  max_observed   = max_observed,
  nightly_sd     = nightly_sd,
  comparison     = comparison,
  mv_slope       = mv_slope,
  best           = comparison$family[1],
  var_ratio      = stats::var(nights$moths) / mean(nights$moths)
)

dir.create(file.path("data", "processed"), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(fits = fits, comparison = comparison, mv = mv,
             band_table = band_table, summary = summary_list), out_file)

message(sprintf(
  paste0("Conditional distribution of a night's catch, %d nights ",
         "(%d zero-catch, raw variance/mean %.0f).\n%s\n",
         "  best: %s | empirical variance-mean slope %.2f ",
         "(1 = Poisson-like, 2 = NB2-like)\n",
         "  nightly multiplier: lognormal SD %.2f, so a one-sigma night is a ",
         "factor of %.0f\n",
         "  written to %s"),
  nrow(nights), zeros_observed, summary_list$var_ratio,
  paste(sprintf("    %-11s %2.0f par | AIC %8.0f | dAIC %7.0f | zeros %5.0f (obs %d) | max %9.0f (obs %d)",
                comparison$family, comparison$parameters, comparison$AIC,
                comparison$dAIC, comparison$zeros_simulated, zeros_observed,
                comparison$max_simulated, max_observed),
        collapse = "\n"),
  summary_list$best, mv_slope, nightly_sd, exp(nightly_sd), out_file
))
