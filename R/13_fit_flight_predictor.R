# =============================================================================
# 13_fit_flight_predictor.R
# -----------------------------------------------------------------------------
# Fit and store the nightly flight predictor, and calibrate its width.
#
# WHAT THIS PRODUCES
#   A single saved object holding everything needed to predict the nightly
#   flight of a locality-year that was never trapped, from the three inputs a
#   forecast has: a distribution of developmental dates, hourly temperature,
#   and the season's total catch. `R/functions/predict_flight.R` consumes it.
#
#   Three driver variants are fitted, so the predictor can be used with
#   whichever input is to hand:
#     flight     BioSIM's daily flight percentage
#     emergence  BioSIM's daily adult emergence percentage
#     pupation   bayessbw's pupation-date distribution
#
# THE SHARPENING PARAMETER
#   @sec-forecast shows every driver produces a season about two and a half
#   days too broad. The fitted shape is a log-share curve; raising it to a
#   power kappa > 1 sharpens it around its own peak without moving that peak.
#   kappa is chosen out of sample by stacking: for each year, kappa is the
#   value that scores best on the OTHER years' held-out predictions, so no
#   year contributes to the kappa used on itself.
#
#   Two criteria are reported, because they disagree: the log score is the
#   proper one for nightly counts, and width matching is what makes the
#   predicted season the right length. The saved object carries both.
#
# INPUT   data/processed/flight_timing.rds  (per-night held-out shares from R/12)
# OUTPUT  data/processed/flight_predictor.rds
#           list(models, kappa, calibration, summary)
#
# REQUIREMENTS  dplyr, splines, glmmTMB
#
# RUNTIME  a couple of minutes: kappa reuses the predictions of R/12 and only
#          three final fits are made.
#
# USAGE  source("R/13_fit_flight_predictor.R")   # from the project root
# =============================================================================

library(dplyr)
library(splines)
library(glmmTMB)


## ---- 1. Configuration -------------------------------------------------------

out_file <- file.path("data", "processed", "flight_predictor.rds")

## Sharpening factors searched. 1 leaves the fitted shape untouched.
kappa_grid <- seq(1, 2.5, by = 0.05)

## The specification each driver is fitted with. Temperature is included: it
## improves the nightly log score (@sec-forecast) and is an available input.
df_corr <- 4

drivers <- list(
  flight    = list(curve = "flight_reg", median = "p50_reg",   spec = "flight + temp"),
  emergence = list(curve = "emergence",  median = "p50_emerg", spec = "emergence + temp"),
  pupation  = list(curve = "flight_pup", median = "p50_pup",   spec = "pupation + temp")
)


## ---- 2. Data ----------------------------------------------------------------

timing <- readRDS(file.path("data", "processed", "flight_timing.rds"))
dat    <- timing$dat
pred   <- timing$pred_night          # held-out shares, one row per night per spec


## ---- 3. Choose the sharpening factor by stacking ----------------------------
## Sharpening acts on the held-out shares themselves, so no model is refitted:
## share_kappa is proportional to share^kappa, renormalised within the
## locality-year.

sharpen <- function(share, group, kappa) {
  s <- exp(kappa * log(pmax(share, 1e-12)))
  s / stats::ave(s, group, FUN = sum)
}

season_width <- function(doy, value) {
  keep <- value > 0
  if (sum(keep) < 2) return(NA_real_)
  cum <- cumsum(value[keep]) / sum(value[keep])
  q <- stats::approx(c(0, cum), c(min(doy[keep]) - 1, doy[keep]),
                     xout = c(0.05, 0.95), ties = "ordered")$y
  q[2] - q[1]
}

## Score every (spec, kappa, year) once: the log score per moth and the mean
## signed width error over that year's locality-years.
grid_scores <- bind_rows(lapply(names(drivers), function(nm) {
  sp <- drivers[[nm]]$spec
  p  <- filter(pred, spec == sp)
  bind_rows(lapply(kappa_grid, function(k) {
    p %>%
      mutate(sk = sharpen(share, site_year, k)) %>%
      group_by(Year, site_year) %>%
      summarise(total = sum(moths),
                loglik = sum(moths * log(pmax(sk, 1e-12))),
                w_err  = season_width(doy, sk) - season_width(doy, moths),
                .groups = "drop") %>%
      group_by(Year) %>%
      summarise(driver = nm, kappa = k,
                log_score = sum(loglik) / sum(total),
                width_err = mean(w_err, na.rm = TRUE),
                .groups = "drop")
  }))
}))

## Stacking: for each year, the kappa that is best on the other years.
pick_kappa <- function(g, y, criterion) {
  others <- filter(g, Year != y) %>%
    group_by(kappa) %>%
    summarise(log_score = mean(log_score),
              width_err = mean(width_err), .groups = "drop")
  if (criterion == "log_score") {
    others$kappa[which.max(others$log_score)]
  } else {
    others$kappa[which.min(abs(others$width_err))]
  }
}

kappa_by_year <- bind_rows(lapply(names(drivers), function(nm) {
  g <- filter(grid_scores, driver == nm)
  tibble::tibble(driver = nm, Year = sort(unique(g$Year))) %>%
    rowwise() %>%
    mutate(kappa_logscore = pick_kappa(g, Year, "log_score"),
           kappa_width    = pick_kappa(g, Year, "width")) %>%
    ungroup()
}))

## Held-out skill under each choice, against kappa = 1.
calibrated <- bind_rows(lapply(names(drivers), function(nm) {
  sp <- drivers[[nm]]$spec
  p  <- filter(pred, spec == sp) %>%
    left_join(filter(kappa_by_year, driver == nm), by = "Year")
  criteria <- c(none = "one", logscore = "kappa_logscore", width = "kappa_width")
  bind_rows(lapply(names(criteria), function(cn) {
    which_k <- criteria[[cn]]
    k <- if (which_k == "one") rep(1, nrow(p)) else p[[which_k]]
    p %>%
      mutate(sk = sharpen(share, site_year, k[1]),
             kappa_used = k) %>%
      group_by(site_year, Year) %>%
      summarise(total = sum(moths),
                kappa = first(kappa_used),
                loglik = sum(moths * log(pmax(sk, 1e-12))),
                med_err = { p50 <- function(v) {
                              keep <- v > 0
                              cum <- cumsum(v[keep]) / sum(v[keep])
                              stats::approx(c(0, cum),
                                            c(min(doy[keep]) - 1, doy[keep]),
                                            xout = 0.5, ties = "ordered")$y }
                            p50(sk) - p50(moths) },
                w_err = season_width(doy, sk) - season_width(doy, moths),
                .groups = "drop") %>%
      summarise(driver = nm, criterion = cn,
                kappa = mean(kappa),
                `median error, d` = mean(abs(med_err)),
                `width error, d`  = mean(w_err, na.rm = TRUE),
                `log score per moth` = sum(loglik) / sum(total))
  }))
}))


## ---- 4. Final fits on every year --------------------------------------------
## The coefficients the predictor will use. The spline bases are kept in the
## saved object so that a new locality-year is put on exactly the same basis.

fit_driver <- function(nm) {
  d <- drivers[[nm]]
  delta <- dat$doy - dat[[d$median]]
  basis <- ns(delta, df = df_corr)
  X <- data.frame(
    d  = log(pmax(dat[[d$curve]], 1e-3)),
    b  = unclass(basis),
    t1 = dat$temp_c,
    t2 = dat$temp_c^2
  )
  names(X) <- make.names(names(X), unique = TRUE)
  frame <- data.frame(y = dat$moths, sy = factor(dat$site_year), X)
  form <- stats::as.formula(paste("y ~", paste(names(X), collapse = " + "),
                                  "+ (1 | sy)"))
  fit <- glmmTMB(form, data = frame, family = nbinom2())
  list(driver = nm,
       terms  = names(X),
       coef   = glmmTMB::fixef(fit)$cond[-1],
       theta  = sigma(fit),
       basis  = basis,
       temp_mean = mean(dat$temp_night),
       curve_floor = 1e-3,
       kappa_logscore = mean(kappa_by_year$kappa_logscore[kappa_by_year$driver == nm]),
       kappa_width    = mean(kappa_by_year$kappa_width[kappa_by_year$driver == nm]))
}

models <- lapply(names(drivers), fit_driver)
names(models) <- names(drivers)


## ---- 5. Write and report ----------------------------------------------------

summary_list <- list(
  n_nights   = nrow(dat),
  n_ly       = n_distinct(dat$site_year),
  kappa_grid = range(kappa_grid),
  calibrated = calibrated,
  kappa_by_year = kappa_by_year
)

dir.create(file.path("data", "processed"), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(models = models, calibration = calibrated,
             kappa_by_year = kappa_by_year, grid_scores = grid_scores,
             summary = summary_list), out_file)

message(sprintf(
  paste0("Flight predictor fitted on %d nights from %d locality-years.\n",
         "  sharpening chosen out of sample, kappa searched over %.2f-%.2f\n%s\n",
         "  written to %s"),
  nrow(dat), n_distinct(dat$site_year), min(kappa_grid), max(kappa_grid),
  paste(sprintf("    %-10s %-9s kappa %.2f | median %5.2f d | width %+5.2f d | log score %6.3f",
                calibrated$driver, calibrated$criterion, calibrated$kappa,
                calibrated$`median error, d`, calibrated$`width error, d`,
                calibrated$`log score per moth`),
        collapse = "\n"),
  out_file
))
