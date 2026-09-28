# =============================================================================
# 10_nightly_count_model.R
# -----------------------------------------------------------------------------
# A predictive model of the nightly catch curve.
#
# WHAT IS PREDICTED, AND WHAT IS NOT
#   The target is the SHAPE of a locality-year's flight season: the share of
#   its season total caught on each trapped night. Abundance is not predicted
#   -- nothing in this project forecasts population density -- so every model
#   below carries one free level per locality-year, which absorbs the total
#   and leaves only timing to be explained. A Poisson likelihood with such a
#   level per group is equivalent to a multinomial over that group's nights,
#   which is why the fitted values can be read directly as predicted shares.
#
# THE MODELS
#   Each is a Poisson regression of the nightly count with a locality-year
#   level, differing in what drives the shape:
#     clim    a seasonal curve only, ns(day of year), fitted from other years:
#             the climatological baseline a predictor must beat
#     raw     BioSIM's own site curve used as is, no fitted parameters
#     site    the same curve with a free exponent
#     site+   plus a spline correction in days from the modelled median,
#             which lets the fit shift and reshape the modelled curve
#     reg+    the same, driven by ONE regional BioSIM curve per year
#             (network centroid), which R/09 found predicts better
#     reg+lat plus a latitude interaction, to test once more whether
#             site position adds anything after calibration
#     pup+    the same calibration applied to bayessbw's predicted PUPATION
#             distribution instead of BioSIM's flight curve, site by site
#     pupreg+ and to one regional bayessbw run per year
#     reg-med the regional curve stripped out entirely: only the spline in
#             days from its median survives, so the development model
#             contributes a date and nothing else. The gap between this and
#             `reg+` is what the modelled curve's SHAPE is worth.
#
#   The two bayessbw specifications need no lag constant: the correction is a
#   spline in days from the driving curve's own median, so the interval
#   between pupation and the catch is fitted rather than assumed.
#
# HOW THEY ARE SCORED
#   Leave-one-year-out: every model is fitted without the held-out year and
#   predicts the shares of its locality-years. Because the locality-year level
#   cancels when shares are normalised, a held-out locality-year needs no
#   fitted level -- only the covariates. Two metrics:
#     log-lik per night  multinomial log-likelihood of the observed counts
#                        given the predicted shares, per trapped night
#     p50 error          days between the predicted and observed median of
#                        the season's catch
#
#   Nightly weather is deliberately absent. It would be the natural next
#   covariate, and the log records what it could add, but it is not needed to
#   answer how well the curve can be predicted from a development model.
#
# INPUT   data/processed/model_vs_observed.rds, trap_daily.rds,
#         biosim_daily.rds
# OUTPUT  data/interim/regional_daily.rds        regional curve per year
#         data/processed/nightly_count_model.rds list(nights, scores, fits,
#                                                     summary)
#
# RUNTIME  a few minutes, most of it the regional BioSIM requests (cached).
#
# USAGE  source("R/10_nightly_count_model.R")   # from the project root
# =============================================================================

library(dplyr)
library(lubridate)
library(splines)
library(sbwFieldPheno)   # the fitted bayessbw posterior; see R/04


## ---- 1. Configuration -------------------------------------------------------

reg_file   <- file.path("data", "interim",   "regional_daily.rds")
pupreg_file <- file.path("data", "interim",   "regional_pupal.rds")
out_file   <- file.path("data", "processed", "nightly_count_model.rds")

## Degrees of freedom for the two splines: the climatological seasonal curve
## and the correction in days from the modelled median.
df_clim <- 5
df_corr <- 4

## Regional bayessbw run: same settings as R/04, and the same 1 March to
## 31 August window, so that the regional run differs from the site-specific
## one only in its weather.
colony        <- "NB"
window_months <- 3:8
n_post        <- 100
individuals   <- 100
seed          <- 1978

## Days on which every driving curve is evaluated.
doy_grid <- 120:280

## Regional BioSIM point: replicates per year, averaged.
n_reps       <- 5
biosim_url   <- "http://repicea.dynu.net:80/BioSIM/BioSimWeather"
biosim_parms <- "ApplyMortality:0*ApplyAdultMortality:0"

reuse <- TRUE


## ---- 2. Data ----------------------------------------------------------------

cmp <- readRDS(file.path("data", "processed", "model_vs_observed.rds")) %>%
  filter(comparable)

nights <- readRDS(file.path("data", "processed", "trap_daily.rds")) %>%
  filter(status == "counted") %>%
  select(Location = locality, Year = year, doy, moths) %>%
  semi_join(cmp, by = c("Location", "Year"))

biosim_daily <- readRDS(file.path("data", "processed", "biosim_daily.rds")) %>%
  select(Location, Year, doy, flight_site = Flight_mf)

## The same regional point as R/09: the centroid of the trap network.
centroid <- cmp %>%
  distinct(Location, Latitude, Longitude) %>%
  summarise(Latitude = mean(Latitude), Longitude = mean(Longitude))

years <- sort(unique(cmp$Year))


## ---- 3. The regional curve: one BioSIM run per year -------------------------

biosim_daily_curve <- function(lat, lon, year) {
  u <- sprintf(paste0("%s?lat=%f&long=%f&elev=NaN&from=%d&to=%d",
                      "&model=Spruce_Budworm_Biology&Parameters=%s"),
               biosim_url, lat, lon, year, year, biosim_parms)
  x <- try(utils::read.csv(url(u), skip = 1), silent = TRUE)
  if (inherits(x, "try-error")) return(NULL)
  tibble::tibble(
    Year = year,
    doy  = as.integer(format(as.Date(sprintf("%d-%02d-%02d", x$Year, x$Month,
                                             x$Day)), "%j")),
    flight = x$MaleFlight + x$FemaleFlight
  )
}

if (reuse && file.exists(reg_file)) {
  regional_daily <- readRDS(reg_file)
} else {
  regional_daily <- bind_rows(lapply(years, function(yr) {
    reps <- bind_rows(lapply(seq_len(n_reps), function(i)
      biosim_daily_curve(centroid$Latitude, centroid$Longitude, yr)))
    reps %>%
      group_by(Year, doy) %>%
      summarise(flight_reg = mean(flight), .groups = "drop")
  }))
  dir.create(file.path("data", "interim"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(regional_daily, reg_file)
}


## ---- 3b. Pupation curves from bayessbw --------------------------------------
## A predicted pupation distribution is a simulated population of dates. The
## underlying distribution is continuous, so the 100 individual dates of the
## median population are smoothed with a Gaussian kernel and evaluated on
## whole days, then scaled to percent so that the floor below is comparable
## with BioSIM's percentage curves.

pupal_curve <- function(days, grid = doy_grid) {
  d <- stats::density(days, from = min(grid), to = max(grid), n = length(grid))
  tibble::tibble(doy = grid, pct = 100 * d$y / sum(d$y))
}

pupal_population <- readRDS(file.path("data", "processed", "pupal_population.rds"))

pupal_site <- pupal_population %>%
  semi_join(cmp, by = c("Location", "Year")) %>%
  group_by(Location, Year) %>%
  reframe(pupal_curve(pupation_jday)) %>%
  rename(flight_pup = pct)

p50_pup <- pupal_population %>%
  semi_join(cmp, by = c("Location", "Year")) %>%
  group_by(Location, Year) %>%
  summarise(p50_pup = stats::median(pupation_jday), .groups = "drop")

## The regional run: bayessbw driven by the hourly temperature averaged over
## every locality, exactly as in R/09 (that script kept only the median, so
## the full population is recomputed here).
median_population <- function(m) apply(apply(m, 2, sort), 1, median)

regional_weather <- function(yr, weather_hourly) {
  weather_hourly %>%
    filter(Year == yr) %>%
    select(Weather) %>%
    tidyr::unnest(Weather) %>%
    mutate(Date = as.Date(datetime_est), Hour = hour(datetime_est)) %>%
    filter(month(Date) %in% window_months) %>%
    group_by(Date, Hour) %>%
    summarise(Temp = mean(Temp), .groups = "drop") %>%
    arrange(Date, Hour)
}

if (reuse && file.exists(pupreg_file)) {
  regional_pupal <- readRDS(pupreg_file)
} else {
  weather_hourly <- readRDS(file.path("data", "processed", "weather_hourly.rds"))
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  regional_pupal <- bind_rows(lapply(years, function(yr) {
    m <- dev_days(weather = regional_weather(yr, weather_hourly),
                  sbwcolony = colony, period = "hour", stage = "Pupa",
                  ecdf = FALSE, n.post = n_post, individuals = individuals)
    pop <- median_population(m)
    pupal_curve(pop) %>% mutate(Year = yr, p50_pupreg = stats::median(pop))
  }))
  dir.create(file.path("data", "interim"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(regional_pupal, pupreg_file)
}

regional_pupal_curve <- select(regional_pupal, Year, doy, flight_pupreg = pct)
p50_pupreg <- distinct(regional_pupal, Year, p50_pupreg)


## ---- 4. One row per trapped night -------------------------------------------
## Each night carries the modelled flight percentage from both curves and its
## distance from each curve's median. A small floor keeps log() finite on
## nights the model puts outside the flight period.

floor_pct <- 1e-3

curve_p50 <- function(doy, value) {
  keep <- value > 0
  if (!any(keep)) return(NA_real_)
  cum <- cumsum(value[keep]) / sum(value[keep])
  stats::approx(c(0, cum), c(min(doy[keep]) - 1, doy[keep]),
                xout = 0.5, ties = "ordered")$y
}

curve_pct1 <- function(doy, value, prob) {
  keep <- value > 0
  if (!any(keep)) return(NA_real_)
  cum <- cumsum(value[keep]) / sum(value[keep])
  stats::approx(c(0, cum), c(min(doy[keep]) - 1, doy[keep]),
                xout = prob, ties = "ordered")$y
}

p50_site <- biosim_daily %>%
  group_by(Location, Year) %>%
  summarise(p50_site = curve_p50(doy, flight_site), .groups = "drop")

p50_reg <- regional_daily %>%
  group_by(Year) %>%
  summarise(p50_reg = curve_p50(doy, flight_reg), .groups = "drop")

nights <- nights %>%
  left_join(biosim_daily, by = c("Location", "Year", "doy")) %>%
  left_join(regional_daily, by = c("Year", "doy")) %>%
  left_join(pupal_site, by = c("Location", "Year", "doy")) %>%
  left_join(regional_pupal_curve, by = c("Year", "doy")) %>%
  left_join(p50_site, by = c("Location", "Year")) %>%
  left_join(p50_reg, by = "Year") %>%
  left_join(p50_pup, by = c("Location", "Year")) %>%
  left_join(p50_pupreg, by = "Year") %>%
  left_join(select(cmp, Location, Year, Latitude), by = c("Location", "Year")) %>%
  mutate(site_year   = paste(Location, Year),
         ## Nights outside the evaluation grid get the floor, as do nights the
         ## curves put outside their season.
         log_site    = log(pmax(flight_site, floor_pct)),
         log_reg     = log(pmax(flight_reg,  floor_pct)),
         log_pup     = log(pmax(coalesce(flight_pup,    0), floor_pct)),
         log_pupreg  = log(pmax(coalesce(flight_pupreg, 0), floor_pct)),
         delta_site  = doy - p50_site,
         delta_reg   = doy - p50_reg,
         delta_pup   = doy - p50_pup,
         delta_pupreg = doy - p50_pupreg,
         lat_c       = Latitude - mean(Latitude))

stopifnot(!any(is.na(nights$flight_site)), !any(is.na(nights$flight_reg)),
          !any(is.na(nights$delta_pup)),   !any(is.na(nights$delta_pupreg)))


## ---- 5. Model specifications -------------------------------------------------
## Spline bases are built once over all nights so that a fitted model and a
## held-out prediction share the same knots. The bases use only the timing
## covariates, never the counts.

basis_clim <- ns(nights$doy,        df = df_clim)
basis_corr <- ns(nights$delta_site, df = df_corr)
basis_creg <- ns(nights$delta_reg,  df = df_corr)
basis_cpup <- ns(nights$delta_pup,    df = df_corr)
basis_cpre <- ns(nights$delta_pupreg, df = df_corr)

covariates <- list(
  clim      = basis_clim,
  raw       = NULL,                                   # offset only
  site      = cbind(log_site = nights$log_site),
  `site+`   = cbind(log_site = nights$log_site, basis_corr),
  `reg+`    = cbind(log_reg  = nights$log_reg,  basis_creg),
  `reg+lat` = cbind(log_reg  = nights$log_reg,  basis_creg,
                    lat = nights$lat_c,
                    basis_creg * nights$lat_c),
  `pup+`    = cbind(log_pup    = nights$log_pup,    basis_cpup),
  `pupreg+` = cbind(log_pupreg = nights$log_pupreg, basis_cpre),
  `reg-med` = basis_creg
)

## Predicted shares within each locality-year, from the covariate part only.
shares_from <- function(eta, group) {
  s <- exp(eta - stats::ave(eta, group, FUN = max))   # scale-free within group
  s / stats::ave(s, group, FUN = sum)
}

fit_predict <- function(name, train, test) {
  X <- covariates[[name]]
  if (is.null(X)) {                      # `raw`: BioSIM's curve, unfitted
    return(shares_from(nights$log_site[test], nights$site_year[test]))
  }
  X <- as.matrix(X)
  fit <- stats::glm(nights$moths[train] ~ X[train, , drop = FALSE] +
                      factor(nights$site_year[train]),
                    family = stats::poisson())
  b <- stats::coef(fit)[seq_len(ncol(X) + 1)]          # intercept + covariates
  b <- b[-1]
  b[is.na(b)] <- 0
  eta <- as.vector(X[test, , drop = FALSE] %*% b)
  shares_from(eta, nights$site_year[test])
}


## ---- 6. Leave-one-year-out scoring -------------------------------------------

score_rows <- function(p, idx) {
  d <- nights[idx, ] %>% mutate(p_hat = p)
  d %>%
    group_by(site_year, Location, Year) %>%
    summarise(
      nights_n = n(),
      total    = sum(moths),
      loglik   = sum(moths * log(pmax(p_hat, 1e-12))) / n(),
      obs_p50  = curve_p50(doy, moths),
      pred_p50 = curve_p50(doy, p_hat),
      obs_p05  = curve_pct1(doy, moths, 0.05),
      obs_p95  = curve_pct1(doy, moths, 0.95),
      pred_p05 = curve_pct1(doy, p_hat, 0.05),
      pred_p95 = curve_pct1(doy, p_hat, 0.95),
      .groups  = "drop"
    ) %>%
    mutate(err_p50 = pred_p50 - obs_p50)
}

cv <- bind_rows(lapply(names(covariates), function(nm) {
  bind_rows(lapply(years, function(y) {
    test  <- which(nights$Year == y)
    train <- which(nights$Year != y)
    score_rows(fit_predict(nm, train, test), test) %>% mutate(model = nm)
  }))
}))

scores <- cv %>%
  group_by(model) %>%
  summarise(
    `locality-years` = n(),
    `log-lik per night` = sum(loglik * nights_n) / sum(nights_n),
    `mean |p50 error|`  = mean(abs(err_p50)),
    `RMSE p50`          = sqrt(mean(err_p50^2)),
    `mean p50 bias`     = mean(err_p50),
    .groups = "drop"
  ) %>%
  arrange(`RMSE p50`) %>%
  mutate(model = factor(model, levels = names(covariates))) %>%
  arrange(model)


## ---- 7. The fitted correction, for the log ----------------------------------
## The `reg+` model refitted on every year, so the shape correction it applies
## to the regional curve can be drawn.

X_full <- as.matrix(covariates[["reg+"]])
fit_full <- stats::glm(nights$moths ~ X_full + factor(nights$site_year),
                       family = stats::poisson())
b_full <- stats::coef(fit_full)[2:(ncol(X_full) + 1)]
b_full[is.na(b_full)] <- 0

correction <- tibble::tibble(delta = seq(min(nights$delta_reg),
                                         max(nights$delta_reg), by = 0.5)) %>%
  mutate(log_adjust = as.vector(
    predict(basis_creg, newx = delta) %*% b_full[-1]))

## How much does the correction move and reshape the modelled curve? The
## comparison is made on the nights actually trapped, where the fit is
## supported: `raw` is BioSIM's curve as delivered and `reg+` the calibrated
## regional prediction, each reduced to the percentiles of the shares it puts
## on those nights. Evaluating the fitted curve on a wider day grid would only
## show the correction spline extrapolating past the data.

reshape <- cv %>%
  filter(model %in% c("raw", "reg+")) %>%
  select(site_year, Year, model, pred_p05, pred_p50, pred_p95) %>%
  tidyr::pivot_wider(names_from = model,
                     values_from = c(pred_p05, pred_p50, pred_p95)) %>%
  mutate(shift     = .data[["pred_p50_reg+"]] - pred_p50_raw,
         raw_width = pred_p95_raw - pred_p05_raw,
         cal_width = .data[["pred_p95_reg+"]] - .data[["pred_p05_reg+"]],
         widening  = cal_width - raw_width)


## ---- 8. Write and report ----------------------------------------------------

summary_list <- list(
  n_nights   = nrow(nights),
  n_ly       = n_distinct(nights$site_year),
  n_years    = length(years),
  zero_share = mean(nights$moths == 0),
  centroid   = centroid,
  scores     = scores,
  exponent   = unname(b_full[1]),
  shift      = mean(reshape$shift, na.rm = TRUE),
  shift_sd   = stats::sd(reshape$shift, na.rm = TRUE),
  raw_width  = mean(reshape$raw_width, na.rm = TRUE),
  cal_width  = mean(reshape$cal_width, na.rm = TRUE),
  obs_width  = mean(cv$obs_p95[cv$model == 'reg+'] -
                    cv$obs_p05[cv$model == 'reg+'], na.rm = TRUE),
  best       = as.character(scores$model[which.min(scores$`RMSE p50`)])
)

dir.create(file.path("data", "processed"), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(nights = nights, cv = cv, scores = scores,
             correction = correction, reshape = reshape,
             summary = summary_list),
        out_file)

message(sprintf(
  paste0("Nightly count model: %d trapped nights, %d locality-years, %d years",
         " (%.0f%% of nights caught nothing).\n%s\n",
         "  best on median timing: %s | fitted exponent on the regional curve %.2f\n",
         "  on the trapped nights the correction moves the predicted median %+.1f d\n",
         "    (sd %.1f) and changes the predicted central 90%% from %.1f to %.1f d,\n",
         "    against %.1f d observed\n",
         "  written to %s"),
  nrow(nights), summary_list$n_ly, length(years), 100 * summary_list$zero_share,
  paste(sprintf("    %-8s log-lik/night %8.2f | mean |p50 err| %5.2f d | RMSE %5.2f d | bias %+5.2f d",
                scores$model, scores$`log-lik per night`,
                scores$`mean |p50 error|`, scores$`RMSE p50`,
                scores$`mean p50 bias`),
        collapse = "\n"),
  summary_list$best, summary_list$exponent,
  summary_list$shift, summary_list$shift_sd,
  summary_list$raw_width, summary_list$cal_width, summary_list$obs_width,
  out_file
))
