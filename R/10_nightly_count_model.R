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
library(splines)


## ---- 1. Configuration -------------------------------------------------------

reg_file   <- file.path("data", "interim",   "regional_daily.rds")
out_file   <- file.path("data", "processed", "nightly_count_model.rds")

## Degrees of freedom for the two splines: the climatological seasonal curve
## and the correction in days from the modelled median.
df_clim <- 5
df_corr <- 4

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

p50_site <- biosim_daily %>%
  group_by(Location, Year) %>%
  summarise(p50_site = curve_p50(doy, flight_site), .groups = "drop")

p50_reg <- regional_daily %>%
  group_by(Year) %>%
  summarise(p50_reg = curve_p50(doy, flight_reg), .groups = "drop")

nights <- nights %>%
  left_join(biosim_daily, by = c("Location", "Year", "doy")) %>%
  left_join(regional_daily, by = c("Year", "doy")) %>%
  left_join(p50_site, by = c("Location", "Year")) %>%
  left_join(p50_reg, by = "Year") %>%
  left_join(select(cmp, Location, Year, Latitude), by = c("Location", "Year")) %>%
  mutate(site_year   = paste(Location, Year),
         log_site    = log(pmax(flight_site, floor_pct)),
         log_reg     = log(pmax(flight_reg,  floor_pct)),
         delta_site  = doy - p50_site,
         delta_reg   = doy - p50_reg,
         lat_c       = Latitude - mean(Latitude))

stopifnot(!any(is.na(nights$flight_site)), !any(is.na(nights$flight_reg)))


## ---- 5. Model specifications -------------------------------------------------
## Spline bases are built once over all nights so that a fitted model and a
## held-out prediction share the same knots. The bases use only the timing
## covariates, never the counts.

basis_clim <- ns(nights$doy,        df = df_clim)
basis_corr <- ns(nights$delta_site, df = df_corr)
basis_creg <- ns(nights$delta_reg,  df = df_corr)

covariates <- list(
  clim      = basis_clim,
  raw       = NULL,                                   # offset only
  site      = cbind(log_site = nights$log_site),
  `site+`   = cbind(log_site = nights$log_site, basis_corr),
  `reg+`    = cbind(log_reg  = nights$log_reg,  basis_creg),
  `reg+lat` = cbind(log_reg  = nights$log_reg,  basis_creg,
                    lat = nights$lat_c,
                    basis_creg * nights$lat_c)
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


## ---- 8. Write and report ----------------------------------------------------

summary_list <- list(
  n_nights   = nrow(nights),
  n_ly       = n_distinct(nights$site_year),
  n_years    = length(years),
  zero_share = mean(nights$moths == 0),
  centroid   = centroid,
  scores     = scores,
  exponent   = unname(b_full[1]),
  best       = as.character(scores$model[which.min(scores$`RMSE p50`)])
)

dir.create(file.path("data", "processed"), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(nights = nights, cv = cv, scores = scores,
             correction = correction, summary = summary_list),
        out_file)

message(sprintf(
  paste0("Nightly count model: %d trapped nights, %d locality-years, %d years",
         " (%.0f%% of nights caught nothing).\n%s\n",
         "  best on median timing: %s | fitted exponent on the regional curve %.2f\n",
         "  written to %s"),
  nrow(nights), summary_list$n_ly, length(years), 100 * summary_list$zero_share,
  paste(sprintf("    %-8s log-lik/night %8.2f | mean |p50 err| %5.2f d | RMSE %5.2f d | bias %+5.2f d",
                scores$model, scores$`log-lik per night`,
                scores$`mean |p50 error|`, scores$`RMSE p50`,
                scores$`mean p50 bias`),
        collapse = "\n"),
  summary_list$best, summary_list$exponent, out_file
))
