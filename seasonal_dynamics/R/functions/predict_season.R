# =============================================================================
# seasonal_dynamics/R/functions/predict_season.R
# -----------------------------------------------------------------------------
# Predict the nightly distribution of a season's light-trap catch from the
# seasonal model (seasonal_dynamics/R/07_seasonal_predictor.R).
#
#   predict_season(driver_p50, driver_p50_network, season_total, nights_doy,
#                  driver = c("BioSIM flight", "BioSIM pupation",
#                             "bayessbw pupation"))
#
# ARGUMENTS
#   driver_p50          the driver's median date (day of year) at the
#                       locality. For the adopted driver, "BioSIM flight", this
#                       is the median of BioSIM's flight series (R/03,
#                       series Flight_mf).
#   driver_p50_network  the mean of the same date over the trap network that
#                       year. It carries most of the skill (the year-level
#                       slope is about 0.9, the within-year slope about 0.1;
#                       see the analysis log).
#   season_total        number of moths caught over the season
#   nights_doy          the nights (day of year) for which shares are wanted;
#                       shares are normalised over these nights, so they
#                       answer "of this season's catch, how much falls on
#                       each of these nights"
#   driver              which development-model date is supplied. BioSIM
#                       flight is the best centre (log section "A better
#                       centre"); the two pupation drivers are kept for
#                       comparison.
#   model_file          the fitted model, data/processed/seasonal/seasonal_predictor.rds
#
# VALUE  data frame: doy, share, expected (share x season_total), with
#        attributes centre (predicted median catch date), width90 (predicted
#        central 90% of the season), sigma (uncertainty of the centre, d) and
#        driver.
#
# THE CURVE  a logistic season with the predicted centre and width, averaged
#   over a normal error of SD sigma on the centre, plus a small uniform floor
#   (eps) so that no night has a share of zero. No prediction interval for a
#   night's count is provided: the seasonal model has no nightly dispersion
#   term.
# =============================================================================

predict_season <- function(driver_p50, driver_p50_network, season_total,
                           nights_doy,
                           driver = c("BioSIM flight", "BioSIM pupation",
                                      "bayessbw pupation"),
                           model_file = file.path("data", "processed", "seasonal",
                                                  "seasonal_predictor.rds")) {
  driver <- match.arg(driver)
  stopifnot(length(driver_p50) == 1, length(driver_p50_network) == 1,
            length(season_total) == 1, season_total > 0,
            length(nights_doy) >= 2)
  m <- readRDS(model_file)$fits[[driver]]

  centre <- unname(m$centre["(Intercept)"] + m$centre["xY"] * driver_p50_network +
                     m$centre["xW"] * (driver_p50 - driver_p50_network))
  width  <- max(unname(m$width["(Intercept)"] +
                         m$width["log_total"] * log10(season_total)), 1)
  scale  <- width / (2 * stats::qlogis(0.95))

  z  <- seq(-3, 3, by = 0.25)
  wz <- stats::dnorm(z) / sum(stats::dnorm(z))
  d  <- sort(unique(nights_doy))
  share <- Reduce(`+`, Map(function(zz, ww) {
    mu <- centre + m$sigma * zz
    ww * (stats::plogis(d, mu, scale) - stats::plogis(d - 1, mu, scale))
  }, z, wz))
  share <- share / sum(share)
  share <- (1 - m$eps) * share + m$eps / length(share)
  share <- share / sum(share)

  out <- data.frame(doy = d, share = share, expected = share * season_total)
  attr(out, "centre")  <- centre
  attr(out, "width90") <- width
  attr(out, "sigma")   <- unname(m$sigma)
  attr(out, "driver")  <- driver
  out
}
