# =============================================================================
# seasonal_dynamics/R/functions/seasonal_curve.R
# -----------------------------------------------------------------------------
# The symmetric season curve of the Seasonal dynamics predictor, shared by
# seasonal_dynamics/R/07_seasonal_predictor.R and 08_centre_options.R.
#
#   night_shares(doy, mu, w, shape, sigma = 0, eps = 0.001)
#     share of the season on each night `doy`, for a curve with median `mu`
#     and central 90% `w` (days). The share of night d is F(d) - F(d - 1), the
#     convention of R/functions/phenology.R. With sigma > 0 the curve is
#     averaged over a normal error of SD sigma on its centre (fixed grid of
#     standard-normal points). A uniform floor eps keeps every night above
#     zero; shares are normalised over the nights supplied.
# =============================================================================

## Cumulative distributions with median mu, scaled so that the central 90%
## spans w.
season_shapes <- list(
  normal   = function(x, mu, w) stats::pnorm(x, mu, w / (2 * stats::qnorm(0.95))),
  logistic = function(x, mu, w) stats::plogis(x, mu, w / (2 * stats::qlogis(0.95))),
  laplace  = function(x, mu, w) {
    b <- w / (2 * log(10))
    ifelse(x < mu, 0.5 * exp((x - mu) / b), 1 - 0.5 * exp(-(x - mu) / b))
  }
)

## Integration grid for the centre error.
season_z_grid <- seq(-3, 3, by = 0.25)
season_z_wt   <- stats::dnorm(season_z_grid) / sum(stats::dnorm(season_z_grid))

night_shares <- function(doy, mu, w, shape, sigma = 0, eps = 0.001) {
  F <- season_shapes[[shape]]
  one <- function(m) F(doy, m, w) - F(doy - 1, m, w)
  s <- if (sigma > 0) {
    Reduce(`+`, Map(function(z, wt) wt * one(mu + sigma * z),
                    season_z_grid, season_z_wt))
  } else one(mu)
  s <- s / sum(s)
  s <- (1 - eps) * s + eps / length(s)
  s / sum(s)
}
