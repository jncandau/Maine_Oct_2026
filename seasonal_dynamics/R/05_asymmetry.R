# =============================================================================
# seasonal_dynamics/R/05_asymmetry.R
# -----------------------------------------------------------------------------
# Asymmetry of the flight season: choose a measure, then describe it -- is a
# season consistently skewed, and is the variation between seasons a property
# of years, of localities, or sampling noise?
#
# INPUT   data/processed/seasonal/season_percentiles.rds  (01)
#         data/processed/trap_daily.rds                    (R/01, for the
#                                                           resampling)
#
# STEP 1 -- THE MEASURE
#   The asymmetry of 01, (p95 - p50) - (p50 - p05), is in days, so it grows
#   with the width of the season and would inherit the catch-size effect of
#   04. Two scale-free measures are computed instead, both bounded in [-1, 1],
#   positive when the decline is longer than the rise:
#     kelly   = (p90 + p10 - 2 p50) / (p90 - p10)   uses the central 80%
#     bowley  = (p75 + p25 - 2 p50) / (p75 - p25)   uses the central 50%
#   Neither uses p05 or p95, the percentiles most exposed to the four
#   right-censored records and to single stray nights in the tails. Kelly is
#   the primary measure; Bowley, which rests on a central 50% that is often
#   only three or four nights wide, is reported for comparison.
#
# STEP 2 -- DESCRIPTION
#   (a) Location: is the typical season skewed? Intercept of a mixed model
#       with crossed random effects of year and locality, so that seasons
#       sharing a year or a place are not counted as independent evidence.
#   (b) Noise: every record's asymmetry is an estimate from a finite catch.
#       Its sampling error is measured by drawing `n_boot` catches of the
#       record's OWN season total from its OWN nightly shares (multinomial)
#       and recomputing the measure. The mean sampling variance, set against
#       the observed variance, gives the share of the between-season
#       variation that is real (the reliability of the measure).
#   (c) Structure: the variance of the measure split into year, locality and
#       residual components (same mixed model), and its correlation with
#       season width and with catch size.
#
# ANALYSIS SET   records usable for phenology (R/01) with at least 3 positive
#                nights; a season concentrated on one night has no shape.
#
# OUTPUT  data/processed/seasonal/asymmetry.rds
#         list(data, location, components, reliability, correlations, summary)
#
# REQUIREMENTS  dplyr, lme4
#
# RUNTIME  about a minute.
#
# USAGE  source("seasonal_dynamics/R/05_asymmetry.R")   # from the project root
# =============================================================================

library(dplyr)
library(lme4)

source(file.path("R", "functions", "phenology.R"))

out_dir <- file.path("data", "processed", "seasonal")
n_boot  <- 200
set.seed(1978)


## ---- 1. The measures --------------------------------------------------------------

skew_from <- function(p10, p25, p50, p75, p90) {
  c(kelly  = (p90 + p10 - 2 * p50) / (p90 - p10),
    bowley = (p75 + p25 - 2 * p50) / (p75 - p25))
}

season <- readRDS(file.path(out_dir, "season_percentiles.rds")) %>%
  filter(usable_phenology, n_positive >= 3) %>%
  mutate(kelly  = (p90 + p10 - 2 * p50) / (p90 - p10),
         bowley = (p75 + p25 - 2 * p50) / (p75 - p25),
         log_total = log10(season_total))

stopifnot(all(is.finite(season$kelly)), all(is.finite(season$bowley)))


## ---- 2a/2c. Location and variance components ---------------------------------------

mixed <- function(v) {
  f <- stats::as.formula(paste(v, "~ 1 + (1 | year) + (1 | locality)"))
  lmer(f, data = season, REML = TRUE)
}

## A variance component estimated at zero gives lme4's "singular fit"
## message; that is a result (no detectable year or locality variance), not a
## failure, and is reported in the variance table.
fits <- suppressMessages(list(kelly = mixed("kelly"), bowley = mixed("bowley")))

location <- bind_rows(lapply(names(fits), function(v) {
  m  <- fits[[v]]
  est <- unname(lme4::fixef(m)[1])
  se  <- unname(sqrt(diag(as.matrix(stats::vcov(m))))[1])
  tibble::tibble(measure = v, mean_raw = mean(season[[v]]),
                 median_raw = stats::median(season[[v]]),
                 share_positive = mean(season[[v]] > 0),
                 intercept = est, lo = est - 1.96 * se, hi = est + 1.96 * se)
}))

components <- bind_rows(lapply(names(fits), function(v) {
  vc <- as.data.frame(lme4::VarCorr(fits[[v]]))
  tibble::tibble(measure = v, component = vc$grp, variance = vc$vcov) %>%
    mutate(share = variance / sum(variance))
}))


## ---- 2b. Sampling noise: resample each record from itself -------------------------

nights <- readRDS(file.path("data", "processed", "trap_daily.rds")) %>%
  filter(status == "counted", flightID %in% season$flightID) %>%
  arrange(flightID, doy)

boot_one <- function(id, total) {
  n_id <- nights[nights$flightID == id, ]
  prob <- n_id$moths / sum(n_id$moths)
  draws <- stats::rmultinom(n_boot, round(total), prob)
  sk <- apply(draws, 2, function(x) {
    q <- unname(curve_percentiles(n_id$doy, x, c(0.10, 0.25, 0.50, 0.75, 0.90)))
    skew_from(q[1], q[2], q[3], q[4], q[5])
  })
  tibble::tibble(flightID = id,
                 kelly_boot_mean  = mean(sk["kelly", ],  na.rm = TRUE),
                 kelly_se         = stats::sd(sk["kelly", ],  na.rm = TRUE),
                 bowley_boot_mean = mean(sk["bowley", ], na.rm = TRUE),
                 bowley_se        = stats::sd(sk["bowley", ], na.rm = TRUE))
}

boot <- bind_rows(Map(boot_one, season$flightID, season$season_total))
season <- left_join(season, boot, by = "flightID")

## Reliability: the share of the observed variance that is not sampling
## noise. The bootstrap variance is an estimate of each record's measurement
## variance; averaged over records it is the noise part of the total.
reliability <- bind_rows(lapply(c("kelly", "bowley"), function(v) {
  obs_var   <- stats::var(season[[v]])
  noise_var <- mean(season[[paste0(v, "_se")]]^2)
  ## Same split for the residual of the mixed model: how much of what is
  ## left after year and locality is noise.
  resid_var <- components$variance[components$measure == v &
                                     components$component == "Residual"]
  tibble::tibble(measure = v, observed_var = obs_var, noise_var = noise_var,
                 reliability = 1 - noise_var / obs_var,
                 residual_var = resid_var,
                 residual_noise_share = noise_var / resid_var,
                 median_se = stats::median(season[[paste0(v, "_se")]]),
                 ## Bias of the measure from finite catches: resampled mean
                 ## minus the observed value, averaged.
                 mean_bias = mean(season[[paste0(v, "_boot_mean")]] - season[[v]]))
}))

correlations <- bind_rows(lapply(c("kelly", "bowley"), function(v) {
  tibble::tibble(measure = v,
                 with = c("width90", "log10 season total", "latitude",
                          "median catch date", "asymmetry in days (01)"),
                 r = c(stats::cor(season[[v]], season$width90),
                       stats::cor(season[[v]], season$log_total),
                       stats::cor(season[[v]], season$latitude),
                       stats::cor(season[[v]], season$p50),
                       stats::cor(season[[v]], season$asymmetry)))
}))
## For reference: the day-based measure of 01 against width, the reason it
## was set aside.
r_days_width <- stats::cor(season$asymmetry, season$width90)


## ---- 3. Write and report --------------------------------------------------------------

summary_list <- list(n = nrow(season), n_boot = n_boot,
                     n_years = n_distinct(season$year),
                     n_localities = n_distinct(season$locality),
                     r_days_width = r_days_width)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(data = season, location = location, components = components,
             reliability = reliability, correlations = correlations,
             summary = summary_list),
        file.path(out_dir, "asymmetry.rds"))

cat(sprintf("Asymmetry of the flight season: %d locality-years, %d years, %d localities.\n",
            nrow(season), summary_list$n_years, summary_list$n_localities))
cat(sprintf("  day-based asymmetry (01) vs width90: r = %+.2f\n", r_days_width))
for (i in seq_len(nrow(location))) cat(sprintf(
  "  %-6s raw mean %+.3f median %+.3f, %.0f%% positive | mixed-model mean %+.3f [%+.3f, %+.3f]\n",
  location$measure[i], location$mean_raw[i], location$median_raw[i],
  100 * location$share_positive[i], location$intercept[i], location$lo[i], location$hi[i]))
for (v in c("kelly", "bowley")) {
  cc <- filter(components, measure == v)
  cat(sprintf("  %-6s variance: %s\n", v, paste(sprintf("%s %.0f%%", cc$component,
                                                         100 * cc$share), collapse = " | ")))
}
for (i in seq_len(nrow(reliability))) cat(sprintf(
  "  %-6s observed var %.4f | sampling var %.4f | reliability %.2f | median se %.3f | noise / residual %.2f | bias %+.3f\n",
  reliability$measure[i], reliability$observed_var[i], reliability$noise_var[i],
  reliability$reliability[i], reliability$median_se[i],
  reliability$residual_noise_share[i], reliability$mean_bias[i]))
for (v in c("kelly", "bowley")) {
  cc <- filter(correlations, measure == v)
  cat(sprintf("  %-6s r with: %s\n", v, paste(sprintf("%s %+.2f", cc$with, cc$r),
                                              collapse = " | ")))
}
cat("  written to data/processed/seasonal/asymmetry.rds\n")
