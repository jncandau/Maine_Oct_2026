# =============================================================================
# R/functions/predict_flight.R
# -----------------------------------------------------------------------------
# Predict the nightly flight of a locality-year that has not been trapped.
#
# INPUTS, which are the three a forecast has
#   dev_curve     a distribution of developmental dates: data frame with
#                 `doy` and `pct`. Either BioSIM's daily flight or adult
#                 emergence percentage, or bayessbw's pupation distribution.
#                 `pct` need not sum to 100; only its shape is used.
#   hourly_temp   hourly temperature at the locality: data frame with
#                 `datetime` (local standard time) and `Temp` in degrees C.
#   season_total  the number of moths caught over the season.
#
# OUTPUT
#   One row per night requested: the share of the season expected that night,
#   the expected count, and a negative binomial predictive interval.
#
# WHAT THE PREDICTION IS CONDITIONAL ON
#   The shares are normalised over the nights supplied in `nights_doy`, so
#   they answer "of the season's catch, how much falls on each of these
#   nights". Supplying a window shorter than the flight season inflates every
#   share, exactly as a trap that operates for part of the season would.
#
#   `kappa` sharpens the fitted shape around its own peak. "width" uses the
#   value calibrated in R/13 to make the predicted season the right length;
#   "logscore" uses the value that maximises the nightly log score, which is
#   1 (no sharpening); a number applies it directly.
#
# USAGE
#   source("daily_dynamics/R/functions/predict_flight.R")
#   out <- predict_flight(dev_curve    = emergence_1978,
#                         hourly_temp  = weather_1978,
#                         season_total = 12000,
#                         driver       = "emergence")
# =============================================================================


#' Mean temperature of the flight hours of each night
#'
#' @param hourly_temp data frame with `datetime` and `Temp`.
#' @param flight_hours hours of the night that flight is averaged over;
#'   hours after midnight belong to the night that began the evening before.
#' @return data frame with `doy` and `temp_night`.
nightly_flight_temp <- function(hourly_temp,
                                flight_hours = c(20, 21, 22, 23, 0, 1, 2)) {
  stopifnot(all(c("datetime", "Temp") %in% names(hourly_temp)))
  h <- as.integer(format(hourly_temp$datetime, "%H"))
  night <- as.Date(hourly_temp$datetime) - (h < 12)
  keep <- h %in% flight_hours
  agg <- stats::aggregate(list(temp_night = hourly_temp$Temp[keep]),
                          by = list(night = night[keep]), FUN = mean)
  data.frame(doy = as.integer(format(agg$night, "%j")),
             temp_night = agg$temp_night)
}


#' Predict nightly flight from a developmental distribution, temperature and
#' the season total
#'
#' @param dev_curve data frame with `doy` and `pct`.
#' @param hourly_temp data frame with `datetime` and `Temp`.
#' @param season_total total moths caught over the season.
#' @param driver which fitted variant to use: "flight", "emergence" or
#'   "pupation". It must match what `dev_curve` holds.
#' @param nights_doy days of year to predict; defaults to every night the
#'   temperature series covers within the developmental curve's season.
#' @param kappa "width", "logscore", or a number.
#' @param interval width of the predictive interval.
#' @param model_file the object written by R/13_fit_flight_predictor.R.
#' @return data frame: doy, temp_night, share, expected, lower, upper.
predict_flight <- function(dev_curve, hourly_temp, season_total,
                           driver = c("flight", "emergence", "pupation"),
                           nights_doy = NULL,
                           kappa = "width",
                           interval = 0.90,
                           model_file = file.path("data", "processed",
                                                  "flight_predictor.rds")) {
  driver <- match.arg(driver)
  ## The stored shape is a natural-spline basis, and predict() dispatches on
  ## it only once the splines namespace is loaded -- which is not guaranteed
  ## in a fresh session that merely sources this file.
  if (!requireNamespace("splines", quietly = TRUE)) {
    stop("the splines package is required", call. = FALSE)
  }
  stopifnot(all(c("doy", "pct") %in% names(dev_curve)),
            is.numeric(season_total), length(season_total) == 1,
            season_total >= 0)

  fitted_all <- readRDS(model_file)
  m <- fitted_all$models[[driver]]
  if (is.null(m)) stop("no fitted model for driver '", driver, "'", call. = FALSE)

  k <- if (identical(kappa, "width")) m$kappa_width
       else if (identical(kappa, "logscore")) m$kappa_logscore
       else as.numeric(kappa)

  ## The curve's own median is the reference the shape correction is measured
  ## from, exactly as in fitting.
  cv <- dev_curve[order(dev_curve$doy), ]
  pos <- cv$pct > 0
  cum <- cumsum(cv$pct[pos]) / sum(cv$pct[pos])
  p50 <- stats::approx(c(0, cum), c(min(cv$doy[pos]) - 1, cv$doy[pos]),
                       xout = 0.5, ties = "ordered")$y

  temp <- nightly_flight_temp(hourly_temp)
  if (is.null(nights_doy)) {
    nights_doy <- sort(intersect(temp$doy, cv$doy[pos]))
  }
  nights_doy <- sort(unique(nights_doy))
  if (!length(nights_doy)) stop("no nights to predict", call. = FALSE)

  d <- data.frame(doy = nights_doy)
  d$pct        <- stats::approx(cv$doy, cv$pct, xout = d$doy, rule = 2)$y
  d$temp_night <- stats::approx(temp$doy, temp$temp_night, xout = d$doy,
                                rule = 2)$y
  d$delta      <- d$doy - p50

  X <- cbind(d  = log(pmax(d$pct, m$curve_floor)),
             b  = unclass(stats::predict(m$basis, newx = d$delta)),
             t1 = d$temp_night - m$temp_mean,
             t2 = (d$temp_night - m$temp_mean)^2)
  colnames(X) <- m$terms
  stopifnot(identical(colnames(X), names(m$coef)) ||
              length(m$coef) == ncol(X))

  eta <- as.vector(X %*% as.numeric(m$coef))
  sh  <- exp(k * (eta - max(eta)))
  d$share    <- sh / sum(sh)
  d$expected <- d$share * season_total

  a <- (1 - interval) / 2
  d$lower <- stats::qnbinom(a,     mu = d$expected, size = m$theta)
  d$upper <- stats::qnbinom(1 - a, mu = d$expected, size = m$theta)

  ## Attributes are set after the subset: `[.data.frame` drops them.
  out <- d[, c("doy", "temp_night", "share", "expected", "lower", "upper")]
  attr(out, "driver")    <- driver
  attr(out, "kappa")     <- k
  attr(out, "theta")     <- m$theta
  attr(out, "curve_p50") <- p50
  out
}
