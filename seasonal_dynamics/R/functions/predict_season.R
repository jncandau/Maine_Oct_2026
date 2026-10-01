# =============================================================================
# seasonal_dynamics/R/functions/predict_season.R
# -----------------------------------------------------------------------------
# Predict the nightly distribution of a season's light-trap catch from the
# seasonal model (seasonal_dynamics/R/07_seasonal_predictor.R), with an
# optional nightly dispersion interval (seasonal_dynamics/R/10_nightly_dispersion.R).
#
#   predict_season(driver_p50, driver_p50_network, season_total, nights_doy,
#                  driver = c("BioSIM flight", "BioSIM pupation",
#                             "bayessbw pupation"),
#                  interval = c("none", "nightly"))
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
#   interval            "none" (default): share and expected count only.
#                       "nightly": also a prediction interval for a single
#                       night's count, from the beta-binomial dispersion
#                       fitted in R/10 (theta; multinomial is theta = Inf).
#                       Nights are NOT independent of each other under this
#                       interval (they are draws from one Dirichlet-multinomial
#                       season), so summing per-night intervals overstates the
#                       uncertainty of a multi-night total.
#   level               interval coverage, used only when interval = "nightly"
#                       (default 0.90, i.e. a 5th-95th percentile interval)
#   nsim                Monte Carlo draws used for the "nightly" interval
#                       (default 4000). Exact beta-binomial quantiles are not
#                       used because season_total ranges up to ~70,000, making
#                       exact tail summation slow; simulation (Beta(theta p,
#                       theta (1-p)) then Binomial(n, P)) is cheap and avoids
#                       the normal approximation's poor fit at low counts.
#   model_file          the fitted model, data/processed/seasonal/seasonal_predictor.rds
#   dispersion_file     the fitted theta, data/processed/seasonal/nightly_dispersion.rds
#                       (only read when interval = "nightly")
#
# VALUE  data frame: doy, share, expected (share x season_total), and, when
#        interval = "nightly", sd, lwr, upr for a single night's count. The
#        attributes are centre (predicted median catch date), width90
#        (predicted central 90% of the season), sigma (uncertainty of the
#        centre, d), driver, and, when interval = "nightly", theta.
#
# THE CURVE  a logistic season with the predicted centre and width, averaged
#   over a normal error of SD sigma on the centre, plus a small uniform floor
#   (eps) so that no night has a share of zero.
# =============================================================================

predict_season <- function(driver_p50, driver_p50_network, season_total,
                           nights_doy,
                           driver = c("BioSIM flight", "BioSIM pupation",
                                      "bayessbw pupation"),
                           interval = c("none", "nightly"),
                           level = 0.90,
                           nsim = 4000,
                           model_file = file.path("data", "processed", "seasonal",
                                                  "seasonal_predictor.rds"),
                           dispersion_file = file.path("data", "processed", "seasonal",
                                                       "nightly_dispersion.rds")) {
  driver   <- match.arg(driver)
  interval <- match.arg(interval)
  stopifnot(length(driver_p50) == 1, length(driver_p50_network) == 1,
            length(season_total) == 1, season_total > 0,
            length(nights_doy) >= 2, level > 0, level < 1)
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

  if (interval == "nightly") {
    theta <- readRDS(dispersion_file)$theta
    n     <- season_total
    alo   <- (1 - level) / 2
    ## Var(Y_d) = n p (1-p) (n + theta) / (1 + theta): the Dirichlet-multinomial
    ## seen one night at a time (R/10). Quantiles are simulated rather than
    ## taken from a normal approximation, which fits poorly when n * share is
    ## small (most nights, away from the peak).
    out$sd <- sqrt(n * share * (1 - share) * (n + theta) / (1 + theta))
    bounds <- vapply(share, function(p) {
      draws <- stats::rbinom(nsim, size = n,
                             prob = stats::rbeta(nsim, theta * p, theta * (1 - p)))
      stats::quantile(draws, c(alo, 1 - alo), names = FALSE)
    }, numeric(2))
    out$lwr <- bounds[1, ]
    out$upr <- bounds[2, ]
    attr(out, "theta") <- theta
  }

  attr(out, "centre")  <- centre
  attr(out, "width90") <- width
  attr(out, "sigma")   <- unname(m$sigma)
  attr(out, "driver")  <- driver
  out
}
