# =============================================================================
# seasonal_dynamics/R/10_population_classes.R
# -----------------------------------------------------------------------------
# Can the season totals of Maine light traps be split into two classes, low
# and high population level? Read on ALL records of the selection window
# (1968-1989), before the 99-moth threshold of R/01, because that threshold
# selects on the very quantity in question.
#
#   SCALE   log10 of the season total. Seasons with no moth at all (a point
#           mass that no continuous distribution describes) are set aside and
#           reported separately; the tests are repeated on log10(total + 1)
#           with them included.
#   TESTS   (a) Hartigan's dip test of unimodality;
#           (b) Gaussian mixtures with 1 to 4 components, equal ("E") or
#               unequal ("V") variances, compared by BIC;
#           (c) parametric bootstrap likelihood-ratio test of 1 against 2
#               and 2 against 3 components (`n_boot` resamples).
#   CLASSES the two-component fit: means, SDs, proportions, the boundary
#           where a season is equally likely to be low or high, and the share
#           of seasons classified with posterior probability >= 0.9.
#   WHERE THE LEVEL LIES  variance of log10(total + 1) between years, between
#           localities and within locality-year (mixed model), and the share
#           of high seasons in each year: a population level should be a
#           property of the year and place, not of the trap-night.
#   EFFORT  a small total can mean a short trapping season rather than few
#           moths: the tests are repeated on records with a daily series of
#           at least `min_nights` nights.
#
# INPUT   data/processed/trap_season_totals.rds   (R/01)
#
# OUTPUT  data/processed/seasonal/population_classes.rds
#         list(data, tests, bic, fit2, by_year, variance, summary)
#
# REQUIREMENTS  dplyr, mclust, diptest, lme4
#
# USAGE  source("seasonal_dynamics/R/10_population_classes.R")   # from the project root
# =============================================================================

library(dplyr)
library(mclust)
library(lme4)

out_dir    <- file.path("data", "processed", "seasonal")
n_boot     <- 999
min_nights <- 20
set.seed(20261002)


## ---- 1. Season totals ---------------------------------------------------------------

tot <- readRDS(file.path("data", "processed", "trap_season_totals.rds")) %>%
  mutate(log_total = ifelse(season_total > 0, log10(season_total), NA_real_),
         log_total1 = log10(season_total + 1))
stopifnot(!anyNA(tot$season_total))


## ---- 2. Tests of one against two classes ------------------------------------------------

tests_on <- function(y, label) {
  y <- y[!is.na(y)]
  dt  <- diptest::dip.test(y)   # not `dip`: tibble() would mask it with its column
  bic <- mclust::mclustBIC(y, G = 1:4, modelNames = c("E", "V"), verbose = FALSE)
  best <- summary(bic, y)
  lrt  <- mclust::mclustBootstrapLRT(y, modelName = "V", maxG = 2, nboot = n_boot,
                                     verbose = FALSE)
  list(
    tests = tibble::tibble(set = label, n = length(y),
                           dip = unname(dt$statistic), dip_p = dt$p.value,
                           best_G = best$G, best_model = best$modelName,
                           lrt_1v2 = lrt$obs[1], lrt_1v2_p = lrt$p.value[1],
                           lrt_2v3 = lrt$obs[2], lrt_2v3_p = lrt$p.value[2]),
    bic = tibble::as_tibble(as.data.frame(unclass(bic)), rownames = "G") %>%
      mutate(set = label, G = as.integer(G)) %>%
      tidyr::pivot_longer(c(E, V), names_to = "model", values_to = "BIC")
  )
}

sets <- list(
  `positive totals, log10`            = tot$log_total,
  `all totals, log10(total + 1)`      = tot$log_total1,
  `daily series >= 20 nights, log10`  = tot$log_total[tot$has_daily & tot$n_nights >= min_nights]
)
res   <- lapply(names(sets), function(nm) tests_on(sets[[nm]], nm))
tests <- bind_rows(lapply(res, `[[`, "tests"))
bic   <- bind_rows(lapply(res, `[[`, "bic"))


## ---- 3. The two-class fit on the positive totals -------------------------------------------

y   <- tot$log_total[!is.na(tot$log_total)]
m2  <- mclust::Mclust(y, G = 2, modelNames = "V", verbose = FALSE)
pr  <- m2$parameters
ord <- order(pr$mean)                                  # component 1 = low
fit2 <- tibble::tibble(class = c("low", "high"),
                       mean_log = unname(pr$mean[ord]),
                       sd_log = sqrt(unname(pr$variance$sigmasq[ord])),
                       prop = unname(pr$pro[ord])) %>%
  mutate(median_moths = 10^mean_log)

## Boundary: where the posterior of the high class is 0.5, between the means.
post_high <- function(x) {
  d <- sapply(seq_along(ord), function(k)
    pr$pro[ord[k]] * stats::dnorm(x, pr$mean[ord[k]], sqrt(pr$variance$sigmasq[ord[k]])))
  d[2] / sum(d)
}
boundary_log <- stats::uniroot(function(x) post_high(x) - 0.5,
                               range(fit2$mean_log))$root

tot <- tot %>%
  mutate(p_high = ifelse(is.na(log_total), 0, sapply(log_total, function(x)
                    if (is.na(x)) 0 else post_high(x))),
         class = ifelse(season_total == 0, "none",
                        ifelse(p_high >= 0.5, "high", "low")))
certain_share <- mean(pmax(tot$p_high, 1 - tot$p_high)[tot$season_total > 0] >= 0.9)


## ---- 4. Where the population level lies ---------------------------------------------------

m_var <- lmer(log_total1 ~ 1 + (1 | year) + (1 | locality), data = tot, REML = TRUE)
vc <- as.data.frame(lme4::VarCorr(m_var))
variance <- tibble::tibble(component = recode(vc$grp, year = "year", locality = "locality",
                                              Residual = "locality-year"),
                           var = vc$vcov) %>%
  mutate(share = var / sum(var), sd = sqrt(var))

by_year <- tot %>%
  group_by(year) %>%
  summarise(n = n(), n_zero = sum(season_total == 0), n_high = sum(class == "high"),
            share_high = n_high / n, median_total = stats::median(season_total),
            .groups = "drop")


## ---- 5. Write and report ----------------------------------------------------------------

summary_list <- list(
  n = nrow(tot), n_zero = sum(tot$season_total == 0), n_positive = length(y),
  n_daily = sum(tot$has_daily), n_total_only = sum(!tot$has_daily),
  n_retained = sum(tot$retained), threshold = 99, n_boot = n_boot,
  min_nights = min_nights, boundary_log = boundary_log,
  boundary_moths = 10^boundary_log, certain_share = certain_share,
  n_high = sum(tot$class == "high"), n_low = sum(tot$class == "low"),
  high_retained = sum(tot$class == "high" & tot$season_total > 99),
  low_above_99 = sum(tot$class == "low" & tot$season_total > 99),
  cor_nights = with(filter(tot, has_daily, season_total > 0),
                    stats::cor(log_total, log10(n_nights)))
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(data = tot, tests = tests, bic = bic, fit2 = fit2, by_year = by_year,
             variance = variance, summary = summary_list),
        file.path(out_dir, "population_classes.rds"))

s <- summary_list
cat(sprintf("Population classes: %d records 1968-1989 (%d with a daily series, %d totals only); %d with no moth, %d retained by R/01.\n",
            s$n, s$n_daily, s$n_total_only, s$n_zero, s$n_retained))
cat("  tests:\n")
for (i in seq_len(nrow(tests))) cat(sprintf(
  "    %-34s n %3d | dip %.3f p %.3f | BIC best G %d (%s) | LRT 1v2 %6.1f p %.3f | 2v3 %5.1f p %.3f\n",
  tests$set[i], tests$n[i], tests$dip[i], tests$dip_p[i], tests$best_G[i], tests$best_model[i],
  tests$lrt_1v2[i], tests$lrt_1v2_p[i], tests$lrt_2v3[i], tests$lrt_2v3_p[i]))
cat("  BIC (positive totals):\n")
print(bic %>% filter(set == names(sets)[1]) %>% select(G, model, BIC) %>%
        tidyr::pivot_wider(names_from = model, values_from = BIC))
cat("  two-class fit (positive totals, log10):\n")
for (i in seq_len(nrow(fit2))) cat(sprintf(
  "    %-4s mean %.2f (%.0f moths) | SD %.2f | proportion %.2f\n",
  fit2$class[i], fit2$mean_log[i], fit2$median_moths[i], fit2$sd_log[i], fit2$prop[i]))
cat(sprintf("    boundary %.2f (%.0f moths); %.0f%% of seasons classified with p >= 0.9\n",
            s$boundary_log, s$boundary_moths, 100 * s$certain_share))
cat(sprintf("    low %d, high %d; high seasons above 99 moths %d; low seasons above 99 moths %d\n",
            s$n_low, s$n_high, s$high_retained, s$low_above_99))
cat(sprintf("  cor(log total, log nights), daily records: %.2f\n", s$cor_nights))
cat("  variance of log10(total + 1):\n")
for (i in seq_len(nrow(variance))) cat(sprintf("    %-14s SD %.2f | share %.2f\n",
  variance$component[i], variance$sd[i], variance$share[i]))
cat("  by year:\n")
for (i in seq_len(nrow(by_year))) cat(sprintf(
  "    %d  n %2d | zero %2d | high %2d (%.2f) | median %7.0f\n", by_year$year[i], by_year$n[i],
  by_year$n_zero[i], by_year$n_high[i], by_year$share_high[i], by_year$median_total[i]))
cat("  written to data/processed/seasonal/population_classes.rds\n")
