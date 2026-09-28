# =============================================================================
# 12_flight_timing.R
# -----------------------------------------------------------------------------
# The timing and spread of nightly flight, from the inputs a forecast has.
#
# THE FORECASTING PROBLEM
#   At prediction time three things are available for a locality-year:
#     1. a distribution of developmental dates -- pupation from bayessbw, or
#        adult emergence from BioSIM;
#     2. hourly temperature at the locality;
#     3. the season's total catch.
#   Nothing else: no trap data from the year being predicted.
#
#   The season total is an input, so abundance is not modelled. What is
#   modelled is how that total is distributed over the nights: a share per
#   night, which is the timing and the shape of the flight season, and the
#   scatter of the count around it.
#
# WHAT IS NEW HERE, RELATIVE TO R/10
#   * BioSIM's EMERGENCE series (MaleEmergence + FemaleEmergence), which is
#     the distribution of adult dates. It sums to 100% of the cohort, unlike
#     the `Adults` stock kept in R/03, which is re-counted daily and sums to
#     about 1,500%. It is the BioSIM input that corresponds to what bayessbw
#     supplies for pupation.
#   * Hourly temperature as a covariate: the mean temperature of the flight
#     hours of each night (see `flight_hours`), which no model here has used.
#   * A negative binomial (NB2) likelihood throughout, which R/11 identified
#     as the right conditional distribution.
#   * The spread of the predicted season is scored as well as its median,
#     because a flight window that is centred but too wide is a poor forecast.
#
# HOW IT IS SCORED
#   Leave-one-year-out. The locality-year level cancels when the predicted
#   shares are normalised, so a held-out locality-year needs only its inputs.
#     median error   days between predicted and observed median catch date
#     width error    days between predicted and observed central 90%
#     log score      multinomial log-likelihood per moth
#     coverage       share of nights inside the 90% NB2 predictive interval,
#                    with the season total supplied
#
# INPUT   data/processed/nightly_count_model.rds  (nights, with the BioSIM and
#                                                  bayessbw curves of R/10)
#         data/processed/weather_hourly.rds
# OUTPUT  data/interim/regional_emergence.rds
#         data/processed/flight_timing.rds  list(nights, cv, scores, summary)
#
# REQUIREMENTS  dplyr, lubridate, splines, glmmTMB, network access to the
#               BioSIM Web API
#
# RUNTIME  about ten minutes: the emergence requests are cached, the
#          leave-one-year-out NB2 fits are the slow part.
#
# USAGE  source("R/12_flight_timing.R")   # from the project root
# =============================================================================

library(dplyr)
library(lubridate)
library(splines)
library(glmmTMB)


## ---- 1. Configuration -------------------------------------------------------

emerg_file <- file.path("data", "interim",   "regional_emergence.rds")
out_file   <- file.path("data", "processed", "flight_timing.rds")

## Hours of the night over which temperature is averaged. Spruce budworm fly
## at and after dusk, so the window runs from early evening to after midnight;
## hour 0-3 of a night belongs to the trap-night that began the evening before.
flight_hours <- c(20, 21, 22, 23, 0, 1, 2)

## Degrees of freedom of the shape correction.
df_corr <- 4

## Regional BioSIM point and replicate count, as in R/09.
n_reps       <- 5
biosim_url   <- "http://repicea.dynu.net:80/BioSIM/BioSimWeather"
biosim_parms <- "ApplyMortality:0*ApplyAdultMortality:0"

reuse <- TRUE


## ---- 2. The nights, with the drivers already built in R/10 ------------------

nights <- readRDS(file.path("data", "processed",
                            "nightly_count_model.rds"))$nights

years <- sort(unique(nights$Year))

## The regional point of R/09: the centroid of the trap network. The nights
## table carries latitude only, so the coordinates come from the comparison
## table.
centroid <- readRDS(file.path("data", "processed", "model_vs_observed.rds")) %>%
  filter(comparable) %>%
  distinct(Location, Latitude, Longitude) %>%
  summarise(Latitude = mean(Latitude), Longitude = mean(Longitude))


## ---- 3. BioSIM emergence at the regional point ------------------------------
## The distribution of adult dates: a genuine daily flux that sums to the whole
## cohort, so its cumulative curve is a distribution over the season.

emergence_curve <- function(lat, lon, year) {
  u <- sprintf(paste0("%s?lat=%f&long=%f&elev=NaN&from=%d&to=%d",
                      "&model=Spruce_Budworm_Biology&Parameters=%s"),
               biosim_url, lat, lon, year, year, biosim_parms)
  x <- try(utils::read.csv(url(u), skip = 1), silent = TRUE)
  if (inherits(x, "try-error")) return(NULL)
  tibble::tibble(
    Year = year,
    doy  = as.integer(format(as.Date(sprintf("%d-%02d-%02d", x$Year, x$Month,
                                             x$Day)), "%j")),
    emergence = x$MaleEmergence + x$FemaleEmergence
  )
}

if (reuse && file.exists(emerg_file)) {
  regional_emergence <- readRDS(emerg_file)
} else {
  regional_emergence <- bind_rows(lapply(years, function(yr) {
    bind_rows(lapply(seq_len(n_reps),
                     function(i) emergence_curve(centroid$Latitude,
                                                 centroid$Longitude, yr))) %>%
      group_by(Year, doy) %>%
      summarise(emergence = mean(emergence), .groups = "drop")
  }))
  dir.create(file.path("data", "interim"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(regional_emergence, emerg_file)
}


## ---- 4. Nightly temperature --------------------------------------------------
## The mean temperature of the flight hours of each trap-night. Hours after
## midnight are assigned to the night that began the previous evening.

weather_hourly <- readRDS(file.path("data", "processed", "weather_hourly.rds"))

nightly_temp <- weather_hourly %>%
  select(Location, Year, Weather) %>%
  tidyr::unnest(Weather) %>%
  mutate(hour = hour(datetime_est),
         ## a night is labelled by the day its evening falls on
         night = as.Date(datetime_est) - (hour < 12)) %>%
  filter(hour %in% flight_hours) %>%
  group_by(Location, Year, night) %>%
  summarise(temp_night = mean(Temp), n_hours = n(), .groups = "drop") %>%
  filter(n_hours >= 5) %>%
  mutate(doy = yday(night)) %>%
  select(Location, Year, doy, temp_night)


## ---- 5. Assemble ------------------------------------------------------------

curve_pct <- function(doy, value, prob) {
  keep <- value > 0
  if (!any(keep)) return(NA_real_)
  cum <- cumsum(value[keep]) / sum(value[keep])
  stats::approx(c(0, cum), c(min(doy[keep]) - 1, doy[keep]),
                xout = prob, ties = "ordered")$y
}

p50_emerg <- regional_emergence %>%
  group_by(Year) %>%
  summarise(p50_emerg = curve_pct(doy, emergence, 0.5), .groups = "drop")

floor_pct <- 1e-3

dat <- nights %>%
  left_join(select(regional_emergence, Year, doy, emergence),
            by = c("Year", "doy")) %>%
  left_join(p50_emerg, by = "Year") %>%
  left_join(nightly_temp, by = c("Location", "Year", "doy")) %>%
  mutate(log_emerg   = log(pmax(coalesce(emergence, 0), floor_pct)),
         delta_emerg = doy - p50_emerg,
         temp_c      = temp_night - mean(temp_night, na.rm = TRUE))

stopifnot(!any(is.na(dat$delta_emerg)))
temp_missing <- sum(is.na(dat$temp_night))
dat <- filter(dat, !is.na(temp_night))


## ---- 6. Specifications -------------------------------------------------------
## Every specification predicts the share of the season's catch on a night.
## `driver` names the developmental distribution used; `+ temp` adds the
## nightly flight-hour temperature.

bases <- list(
  emerg = ns(dat$delta_emerg,  df = df_corr),
  flight = ns(dat$delta_reg,   df = df_corr),
  pup    = ns(dat$delta_pup,   df = df_corr)
)

specs <- list(
  `emergence`          = cbind(d = dat$log_emerg, bases$emerg),
  `emergence + temp`   = cbind(d = dat$log_emerg, bases$emerg,
                               t1 = dat$temp_c, t2 = dat$temp_c^2),
  `flight (R/10 best)` = cbind(d = dat$log_reg,   bases$flight),
  `flight + temp`      = cbind(d = dat$log_reg,   bases$flight,
                               t1 = dat$temp_c, t2 = dat$temp_c^2),
  `pupation`           = cbind(d = dat$log_pup,   bases$pup),
  `pupation + temp`    = cbind(d = dat$log_pup,   bases$pup,
                               t1 = dat$temp_c, t2 = dat$temp_c^2),
  `temperature only`   = cbind(t1 = dat$temp_c, t2 = dat$temp_c^2,
                               ns(dat$doy, df = 5))
)


## ---- 7. Fit and score leave-one-year-out -------------------------------------

shares_from <- function(eta, group) {
  s <- exp(eta - stats::ave(eta, group, FUN = max))
  s / stats::ave(s, group, FUN = sum)
}

fit_fold <- function(X, train, test) {
  d <- data.frame(y = dat$moths, sy = factor(dat$site_year), X)
  cols <- setdiff(names(d), c("y", "sy"))
  f <- stats::as.formula(paste("y ~", paste(cols, collapse = " + "),
                               "+ (1 | sy)"))
  fit <- try(glmmTMB(f, data = d[train, ], family = nbinom2()), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  b <- glmmTMB::fixef(fit)$cond[-1]
  b[is.na(b)] <- 0
  eta <- as.vector(as.matrix(d[test, cols, drop = FALSE]) %*% b)
  list(shares = shares_from(eta, dat$site_year[test]), theta = sigma(fit))
}

score_fold <- function(pred, idx) {
  d <- dat[idx, ] %>% mutate(p_hat = pred$shares)
  d %>%
    group_by(site_year, Location, Year) %>%
    summarise(
      nights_n = n(), total = sum(moths),
      loglik   = sum(moths * log(pmax(p_hat, 1e-12))),
      obs_p50  = curve_pct(doy, moths, 0.5),
      pred_p50 = curve_pct(doy, p_hat, 0.5),
      obs_w    = curve_pct(doy, moths, 0.95) - curve_pct(doy, moths, 0.05),
      pred_w   = curve_pct(doy, p_hat, 0.95) - curve_pct(doy, p_hat, 0.05),
      ## coverage of the 90% predictive interval, the season total supplied
      covered  = mean(moths >= stats::qnbinom(0.05, mu = p_hat * sum(moths),
                                              size = pred$theta) &
                      moths <= stats::qnbinom(0.95, mu = p_hat * sum(moths),
                                              size = pred$theta)),
      .groups = "drop"
    ) %>%
    mutate(err_p50 = pred_p50 - obs_p50, err_w = pred_w - obs_w)
}

## Per-night predictions are kept as well as the summaries: the width of a
## predicted season has to be compared with a realisation of it, not with a
## mean curve, which needs the nightly shares and the fitted dispersion.
folds <- lapply(names(specs), function(nm) {
  X <- as.data.frame(specs[[nm]])
  names(X) <- make.names(names(X), unique = TRUE)
  out <- lapply(years, function(y) {
    test <- which(dat$Year == y); train <- which(dat$Year != y)
    p <- fit_fold(X, train, test)
    if (is.null(p)) return(NULL)
    list(score = score_fold(p, test) %>% mutate(spec = nm),
         nights = dat[test, c("site_year", "Year", "doy", "moths")] %>%
           mutate(spec = nm, share = p$shares, theta = p$theta))
  })
  out <- out[!vapply(out, is.null, logical(1))]
  list(score  = bind_rows(lapply(out, `[[`, "score")),
       nights = bind_rows(lapply(out, `[[`, "nights")))
})

cv         <- bind_rows(lapply(folds, `[[`, "score"))
pred_night <- bind_rows(lapply(folds, `[[`, "nights"))

scores <- cv %>%
  group_by(spec) %>%
  summarise(`locality-years` = n(),
            `median error, d`   = mean(abs(err_p50)),
            `RMSE median, d`    = sqrt(mean(err_p50^2)),
            `width error, d`    = mean(err_w),
            `|width error|, d`  = mean(abs(err_w)),
            `log score per moth` = sum(loglik) / sum(total),
            `90% coverage`       = mean(covered),
            .groups = "drop") %>%
  arrange(`median error, d`)


## ---- 7b. Is the predicted season really too wide? ---------------------------
## The observed width is measured on ONE realisation of a very overdispersed
## count, and a single spiky realisation is narrower than the smooth mean
## curve behind it. The honest comparison simulates realisations from the
## fitted predictive distribution and measures their widths the same way the
## observed width is measured.

n_draw <- 100

width_of <- function(doy, value) {
  keep <- value > 0
  if (sum(keep) < 2) return(NA_real_)
  cum <- cumsum(value[keep]) / sum(value[keep])
  q <- stats::approx(c(0, cum), c(min(doy[keep]) - 1, doy[keep]),
                     xout = c(0.05, 0.95), ties = "ordered")$y
  q[2] - q[1]
}

set.seed(1978)
width_check <- pred_night %>%
  group_by(spec, site_year) %>%
  group_modify(function(g, key) {
    total <- sum(g$moths)
    mu    <- g$share * total
    sims  <- replicate(n_draw,
                       width_of(g$doy, stats::rnbinom(nrow(g), mu = mu,
                                                      size = g$theta[1])))
    tibble::tibble(
      obs_width       = width_of(g$doy, g$moths),
      mean_width      = width_of(g$doy, mu),          # the smooth mean curve
      sim_width_mean  = mean(sims, na.rm = TRUE),
      sim_width_lo    = stats::quantile(sims, 0.05, na.rm = TRUE, names = FALSE),
      sim_width_hi    = stats::quantile(sims, 0.95, na.rm = TRUE, names = FALSE)
    )
  }) %>%
  ungroup() %>%
  mutate(inside = obs_width >= sim_width_lo & obs_width <= sim_width_hi)

width_summary <- width_check %>%
  group_by(spec) %>%
  summarise(`observed width, d`      = mean(obs_width),
            `mean-curve width, d`    = mean(mean_width),
            `simulated width, d`     = mean(sim_width_mean),
            `inside the 90% band`    = mean(inside),
            .groups = "drop") %>%
  arrange(desc(`inside the 90% band`))


## ---- 8. Write and report -----------------------------------------------------

summary_list <- list(
  n_nights      = nrow(dat),
  n_ly          = n_distinct(dat$site_year),
  temp_missing  = temp_missing,
  scores        = scores,
  width_summary = width_summary,
  best_timing   = scores$spec[1],
  best_width    = scores$spec[which.min(abs(scores$`|width error|, d`))],
  obs_width     = mean(cv$obs_w[cv$spec == scores$spec[1]])
)

dir.create(file.path("data", "processed"), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(dat = dat, cv = cv, scores = scores, pred_night = pred_night,
             width_check = width_check, width_summary = width_summary,
             regional_emergence = regional_emergence,
             summary = summary_list), out_file)

message(sprintf(
  paste0("Flight timing and spread, %d nights from %d locality-years",
         " (%d nights dropped for missing temperature).\n%s\n",
         "  best median timing: %s | observed width %.1f d\n",
         "  width of a REALISATION against the observed width:\n%s\n",
         "  written to %s"),
  nrow(dat), summary_list$n_ly, temp_missing,
  paste(sprintf("    %-19s median %5.2f d | width %+5.2f d | log score %6.3f | coverage %4.0f%%",
                scores$spec, scores$`median error, d`, scores$`width error, d`,
                scores$`log score per moth`, 100 * scores$`90% coverage`),
        collapse = "\n"),
  summary_list$best_timing, summary_list$obs_width,
  paste(sprintf("    %-19s mean curve %5.1f d | simulated %5.1f d | observed %5.1f d | inside %3.0f%%",
                width_summary$spec, width_summary$`mean-curve width, d`,
                width_summary$`simulated width, d`,
                width_summary$`observed width, d`,
                100 * width_summary$`inside the 90% band`),
        collapse = "\n"),
  out_file
))
