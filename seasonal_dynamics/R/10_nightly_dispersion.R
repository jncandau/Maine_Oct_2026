# =============================================================================
# seasonal_dynamics/R/10_nightly_dispersion.R
# -----------------------------------------------------------------------------
# A nightly dispersion term for the seasonal predictor (log, "Next steps"):
# the curve in R/07 and R/functions/predict_season.R gives the *expected*
# share of a season's catch on each night, with no account of how far a real
# night can scatter around that expectation. The resampling used elsewhere in
# this analysis (R/04, R/05, R/06) treats moths as independent draws, i.e. a
# night's count as Binomial(season_total, share) -- the daily_dynamics count
# model instead found nights to be far more variable than that (NB2
# theta ~ 0.19 around a within-season mean curve; data/processed/count_distribution.rds).
# This script asks the same question at the scale of this analysis: around
# the *predicted* nightly share, not the observed one.
#
# MODEL  a night's count is Beta-Binomial(n = season_total, mean = share,
#        precision = theta): the natural two-parameter family whose
#        between-night variance is
#          Var(Y_d) = n p_d (1 - p_d) (n + theta) / (1 + theta)
#        which is the multinomial variance as theta -> Inf, and inflates it
#        as theta falls. It is the Dirichlet-multinomial seen one night at a
#        time, so the shares still sum to season_total exactly.
#
# WHICH SHARE  the leave-one-year-out predictive share already stored by R/07
#   (seasonal_predictor.rds$nights, the "best" row of its score table): this
#   is the out-of-sample prediction predict_season() would itself produce,
#   so theta describes the dispersion a user of predict_season() will
#   actually meet, not an in-sample fit.
#
# ESTIMATION  theta is the only free parameter (the mean curve is fixed), so
#   it is found by a 1-D likelihood search rather than a regression package;
#   this avoids adding a new dependency for a single scalar.
#
# INPUT   data/processed/seasonal/seasonal_predictor.rds  (R/07: nights, the
#         leave-one-year-out predictive share per night)
#         data/processed/seasonal/season_percentiles.rds  (season_total)
#
# OUTPUT  data/processed/seasonal/nightly_dispersion.rds
#         list(theta, loglik, fit, table, summary)
#
# REQUIREMENTS  dplyr
#
# USAGE  source("seasonal_dynamics/R/10_nightly_dispersion.R")   # from the project root
# =============================================================================

library(dplyr)

out_dir  <- file.path("data", "processed", "seasonal")
out_file <- file.path(out_dir, "nightly_dispersion.rds")


## ---- 1. Inputs ---------------------------------------------------------------

predictor <- readRDS(file.path(out_dir, "seasonal_predictor.rds"))
season    <- readRDS(file.path(out_dir, "season_percentiles.rds")) %>%
  select(flightID, season_total)

nights <- predictor$nights %>%
  distinct(flightID, doy, moths, share) %>%
  inner_join(season, by = "flightID")

stopifnot(!anyNA(nights$share), !anyNA(nights$season_total),
          all(nights$share > 0), all(nights$share <= 1))


## ---- 2. Beta-binomial log-likelihood, theta the only free parameter -----------

## lchoose(n, y) is dropped: it does not depend on theta and cancels in any
## comparison between theta values.
bb_loglik <- function(theta, y, n, p) {
  a <- theta * p
  b <- theta * (1 - p)
  sum(lbeta(y + a, n - y + b) - lbeta(a, b))
}

## theta -> Inf is the multinomial (no extra dispersion); theta -> 0 is
## maximal clumping. optimize() on a log-theta grid keeps the search well
## behaved over several orders of magnitude.
nll <- function(log_theta) -bb_loglik(exp(log_theta), nights$moths, nights$season_total, nights$share)
fit <- optimize(nll, interval = log(c(1e-2, 1e6)))
theta_hat <- exp(fit$minimum)

## Likelihood-ratio check against the plain multinomial (theta = Inf): twice
## the gap to the multinomial log-likelihood, compared to chi-sq(1).
ll_mult  <- bb_loglik(1e8, nights$moths, nights$season_total, nights$share)
ll_fit   <- -fit$objective
lrt_stat <- 2 * (ll_fit - ll_mult)
lrt_p    <- stats::pchisq(lrt_stat, df = 1, lower.tail = FALSE)


## ---- 3. What this theta means in practice -------------------------------------

## The variance-inflation factor (n + theta) / (1 + theta) at a few typical
## season totals, and in the pooled data.
vif <- function(n, theta) (n + theta) / (1 + theta)
season_grid <- c(100, 300, 1000, 3000, 10000)
vif_table <- tibble::tibble(
  season_total = season_grid,
  vif = vif(season_grid, theta_hat),
  sd_inflation = sqrt(vif)
)

## Standardised residuals under the fitted model and under the plain
## multinomial, for a quick observed-vs-expected check (mean 0, SD 1 if the
## model is right).
resid_check <- nights %>%
  mutate(
    mu          = season_total * share,
    var_mult    = season_total * share * (1 - share),
    var_bb      = var_mult * vif(season_total, theta_hat),
    z_mult      = (moths - mu) / sqrt(pmax(var_mult, 1e-6)),
    z_bb        = (moths - mu) / sqrt(pmax(var_bb, 1e-6))
  )

by_season <- resid_check %>%
  group_by(flightID) %>%
  summarise(season_total = first(season_total),
            rmse_z_mult = sqrt(mean(z_mult^2)),
            rmse_z_bb   = sqrt(mean(z_bb^2)), .groups = "drop")


## ---- 4. Write and report -------------------------------------------------------

summary_list <- list(
  n_nights    = nrow(nights),
  n_seasons   = n_distinct(nights$flightID),
  theta       = theta_hat,
  loglik      = ll_fit,
  loglik_mult = ll_mult,
  lrt_stat    = lrt_stat,
  lrt_p       = lrt_p,
  rmse_z_mult_pooled = sqrt(mean(resid_check$z_mult^2)),
  rmse_z_bb_pooled   = sqrt(mean(resid_check$z_bb^2)),
  vif_table   = vif_table
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(theta = theta_hat, loglik = ll_fit, by_season = by_season,
             vif_table = vif_table, summary = summary_list),
        out_file)

cat(sprintf(paste0(
  "Nightly dispersion: %d nights, %d locality-years.\n",
  "  beta-binomial theta = %.2f (plain multinomial is theta = Inf)\n",
  "  likelihood-ratio vs. multinomial: stat %.0f, df 1, p %s\n",
  "  RMS standardised residual: multinomial %.2f, beta-binomial %.2f ",
  "(1.0 is correctly scaled)\n",
  "  variance-inflation (n + theta)/(1 + theta) at a few season sizes:\n%s\n",
  "  written to %s\n"),
  summary_list$n_nights, summary_list$n_seasons, theta_hat,
  lrt_stat, if (lrt_p < 1e-10) "< 1e-10" else sprintf("%.2g", lrt_p),
  summary_list$rmse_z_mult_pooled, summary_list$rmse_z_bb_pooled,
  paste(sprintf("    n = %5d : SD x %.1f", vif_table$season_total, vif_table$sd_inflation),
        collapse = "\n"),
  out_file))
