# =============================================================================
# 06_pupal_lag_test.R
# -----------------------------------------------------------------------------
# Would a temperature-driven pupal stage improve the pupation-to-catch step?
#
# THE QUESTION
#   bayessbw predicts pupation; the traps see flight. R/05 bridges the two with
#   a constant: the median observed interval between predicted pupation and the
#   observed median catch (12.9 d). A constant is a poor model of a
#   poikilotherm's development, so the obvious improvement is to accumulate
#   degree-days from pupation and predict the catch median when a threshold is
#   reached. This script tests whether that helps.
#
#     A  predicted catch median = pupation + lag                 (1 parameter)
#     B  predicted catch median = the day K degree-days above    (2 parameters)
#        base T have accumulated since pupation
#
#   Both are fitted on the comparable locality-years, so they are compared by
#   leave-one-year-out cross-validation as well as in sample.
#
# THE CEILING, AND WHERE THE VARIANCE IS
#   Two diagnostics decide whether any version of B could work:
#     - the best achievable fit of the observed interval to post-pupation
#       temperature, using a smooth spline instead of a degree-day form;
#     - a decomposition of the interval's variance into locality, year and
#       residual.
#   A development model can only exploit temperature. If the interval is a
#   property of the site rather than of the temperature it experienced, no
#   development function will recover it.
#
# INPUT   data/processed/model_vs_observed.rds
#         data/processed/weather_hourly.rds
# OUTPUT  data/processed/pupal_lag_test.rds  list(rows, grid_fit, summary)
#
# RUNTIME  about a minute: the degree-day predictions are computed once for
#          the whole (base, K) grid and the cross-validation then only selects
#          from that array.
#
# USAGE  source("R/06_pupal_lag_test.R")   # from the project root
# =============================================================================

library(dplyr)
library(lubridate)


## ---- 1. Configuration -------------------------------------------------------

## Degree-day grid. The base runs from 0 degC (no threshold) to 16 degC, past
## any plausible developmental zero for spruce budworm; K is in degree-days.
bases <- seq(0, 16, by = 0.5)
Ks    <- seq(20, 500, by = 5)

## Window over which "post-pupation temperature" is averaged for the ceiling
## diagnostic, in days after the predicted pupation date.
ceiling_window <- 14

out_file <- file.path("data", "processed", "pupal_lag_test.rds")


## ---- 2. Data: one daily temperature series per locality-year ----------------

mvo <- readRDS(file.path("data", "processed", "model_vs_observed.rds")) %>%
  filter(comparable)

daily <- readRDS(file.path("data", "processed", "weather_hourly.rds")) %>%
  select(Location, Year, Weather) %>%
  rowwise() %>%
  mutate(tmean = list(
    Weather %>%
      mutate(doy = yday(datetime_est)) %>%
      group_by(doy) %>%
      summarise(t = mean(Temp), .groups = "drop") %>%
      arrange(doy)
  )) %>%
  ungroup() %>%
  select(Location, Year, tmean)

rows <- mvo %>%
  select(Location, Year, Latitude, obs_p50, obs_p05, obs_p95,
         sim_p05, sim_p95, season_total,
         pup = bayes_pupation_p50, pup_biosim = biosim_pupation_p50) %>%
  left_join(daily, by = c("Location", "Year")) %>%
  mutate(lag_obs = obs_p50 - pup)

stopifnot(nrow(rows) == nrow(mvo), !any(vapply(rows$tmean, is.null, logical(1))))


## ---- 3. Degree-day predictions over the whole grid --------------------------
## For a given base the accumulation is monotone, so every K can be read off
## one cumulative curve. The crossing day is interpolated, which keeps the
## prediction continuous in K.

dd_curve <- function(tm, start, base) {
  keep <- tm$doy > start
  list(doy = tm$doy[keep], cum = cumsum(pmax(tm$t[keep] - base, 0)))
}

dd_dates <- function(curve, Ks) {
  i <- findInterval(Ks, curve$cum) + 1      # first day at or above each K
  ok <- i <= length(curve$cum)
  out <- rep(NA_real_, length(Ks))
  lo <- ifelse(i[ok] == 1, 0, curve$cum[pmax(i[ok] - 1, 1)])
  out[ok] <- curve$doy[i[ok]] - 1 + (Ks[ok] - lo) / (curve$cum[i[ok]] - lo)
  out
}

## pred[row, K, base]
pred <- array(NA_real_, dim = c(nrow(rows), length(Ks), length(bases)))
for (b in seq_along(bases)) {
  for (i in seq_len(nrow(rows))) {
    pred[i, , b] <- dd_dates(dd_curve(rows$tmean[[i]], rows$pup[i], bases[b]), Ks)
  }
}

## Best (base, K) on a subset of rows, by root mean squared error.
fit_grid <- function(idx) {
  err <- apply(pred[idx, , , drop = FALSE], c(2, 3),
               function(p) sqrt(mean((rows$obs_p50[idx] - p)^2)))
  j <- which(err == min(err, na.rm = TRUE), arr.ind = TRUE)[1, ]
  list(base = bases[j[2]], K = Ks[j[1]], rmse = err[j[1], j[2]])
}

rmse <- function(x) sqrt(mean(x^2))
all_rows <- seq_len(nrow(rows))

fit_B <- fit_grid(all_rows)
lag_A <- median(rows$lag_obs)

rows <- rows %>%
  mutate(pred_A = pup + lag_A,
         pred_B = pred[cbind(all_rows,
                             match(fit_B$K, Ks),
                             match(fit_B$base, bases))])


## ---- 4. Leave-one-year-out cross-validation ---------------------------------
## Both models are refitted without the held-out year. The degree-day fit only
## selects a cell of the precomputed array, so this is cheap.

rows$cv_A <- NA_real_
rows$cv_B <- NA_real_
for (y in unique(rows$Year)) {
  te <- rows$Year == y
  tr <- which(!te)
  rows$cv_A[te] <- rows$pup[te] + median(rows$lag_obs[tr])
  f <- fit_grid(tr)
  rows$cv_B[te] <- pred[cbind(which(te), match(f$K, Ks), match(f$base, bases))]
}


## ---- 5. Could any temperature function do better? ---------------------------

rows <- rows %>%
  rowwise() %>%
  mutate(t_post = mean(tmean$t[tmean$doy > pup & tmean$doy <= pup + ceiling_window]),
         t_july = mean(tmean$t[tmean$doy >= 182 & tmean$doy <= 212])) %>%
  ungroup()

spline_fit <- lm(lag_obs ~ splines::ns(t_post, 5), data = rows)

## Variance of the interval by locality and year. Type-I sums of squares with
## locality first: localities and years are unbalanced here, so this is a
## description of where the variation sits, not a test.
var_parts <- anova(lm(lag_obs ~ factor(Location) + factor(Year), data = rows))
var_share <- 100 * var_parts[["Sum Sq"]] / sum(var_parts[["Sum Sq"]])
names(var_share) <- rownames(var_parts)


## ---- 6. What the catch curve looks like where the interval is long ----------
## If the long intervals are immigration, the catch should carry weight after
## the local cohort's modelled flight window rather than simply being shifted.

rows <- rows %>%
  mutate(asymmetry     = (obs_p95 - obs_p50) - (obs_p50 - obs_p05),
         share_late    = NA_real_,
         share_early   = NA_real_)

trap_daily <- readRDS(file.path("data", "processed", "trap_daily.rds")) %>%
  filter(status == "counted")

for (i in seq_len(nrow(rows))) {
  td <- trap_daily[trap_daily$locality == rows$Location[i] &
                     trap_daily$year == rows$Year[i], ]
  tot <- sum(td$moths)
  rows$share_late[i]  <- if (tot > 0) sum(td$moths[td$doy > rows$sim_p95[i]]) / tot else NA_real_
  rows$share_early[i] <- if (tot > 0) sum(td$moths[td$doy < rows$sim_p05[i]]) / tot else NA_real_
}


## ---- 7. Write and report ----------------------------------------------------

summary_list <- list(
  n            = nrow(rows),
  lag_A        = lag_A,
  fit_B        = fit_B,
  rmse_A       = rmse(rows$obs_p50 - rows$pred_A),
  rmse_B       = rmse(rows$obs_p50 - rows$pred_B),
  cv_rmse_A    = rmse(rows$obs_p50 - rows$cv_A),
  cv_rmse_B    = rmse(rows$obs_p50 - rows$cv_B),
  r_lag_tpost  = cor(rows$lag_obs, rows$t_post),
  r_lag_lat    = cor(rows$lag_obs, rows$Latitude),
  r_lag_tjuly  = cor(rows$lag_obs, rows$t_july),
  spline_r2    = summary(spline_fit)$r.squared,
  spline_sd    = sd(residuals(spline_fit)),
  sd_lag       = sd(rows$lag_obs),
  sd_site      = sd(residuals(lm(lag_obs ~ factor(Location), data = rows))),
  sd_site_year = sd(residuals(lm(lag_obs ~ factor(Location) + factor(Year), data = rows))),
  var_share    = var_share,
  r_late_lag   = cor(rows$share_late, rows$lag_obs),
  r_late_lat   = cor(rows$share_late, rows$Latitude)
)

dir.create(file.path("data", "processed"), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(rows = select(rows, -tmean), grid_fit = fit_B, summary = summary_list),
        out_file)

message(sprintf(
  paste0("Pupation-to-catch interval, %d comparable locality-years.\n",
         "  A fixed lag %.1f d      : RMSE %.2f d in sample, %.2f d leave-one-year-out\n",
         "  B degree-days base %.1f C, K = %g: RMSE %.2f d in sample, %.2f d LOYO\n",
         "  interval vs post-pupation temperature: r = %+.2f ",
         "(spline ceiling R2 = %.3f)\n",
         "  interval sd %.2f d; with a per-site constant %.2f d; site + year %.2f d\n",
         "  variance share: locality %.1f%%, year %.1f%%, residual %.1f%%\n",
         "  written to %s"),
  summary_list$n, lag_A, summary_list$rmse_A, summary_list$cv_rmse_A,
  fit_B$base, fit_B$K, summary_list$rmse_B, summary_list$cv_rmse_B,
  summary_list$r_lag_tpost, summary_list$spline_r2,
  summary_list$sd_lag, summary_list$sd_site, summary_list$sd_site_year,
  var_share[1], var_share[2], var_share[3],
  out_file
))
