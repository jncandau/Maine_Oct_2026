# =============================================================================
# R/functions/phenology.R
# -----------------------------------------------------------------------------
# Helpers shared by the observed (trap) and modelled (BioSIM, bayessbw) sides
# of the analysis, so that a percentile means the same thing everywhere.
#
# USAGE  source("R/functions/phenology.R")   # from the project root
# =============================================================================


## Percentiles of a phenological curve ----------------------------------------
##
## `doy`   day of year
## `value` a flux on that day: moths caught, or the percentage of the cohort
##         in flight. Not a cumulative quantity.
## `probs` the percentiles wanted, e.g. c(0.05, 0.5, 0.95)
##
## The curve is treated as a sample of a continuous distribution over the
## season: values are normalised to sum to one, cumulated, and inverted by
## linear interpolation. The cumulative curve is anchored at zero on the day
## before the first observation, so that p05 can fall before the first day
## with a catch instead of being clamped to it.
##
## Days with NA or negative values are dropped -- for trap data that means
## unrecorded nights are ignored rather than read as zeros, which is why the
## minimum-effort rule in 01_read_trap_data.R matters: percentiles are only
## identifiable when the sampled window brackets the flight period.
##
## Returns a named numeric vector (p05, p50, ...), all NA if the curve is
## empty or sums to zero.
curve_percentiles <- function(doy, value, probs) {
  out <- setNames(rep(NA_real_, length(probs)), sprintf("p%02d", round(probs * 100)))
  keep <- !is.na(value) & value >= 0
  if (!any(keep) || sum(value[keep]) == 0) return(out)
  d <- doy[keep]; v <- value[keep]
  o <- order(d); d <- d[o]; v <- v[o]
  cum <- cumsum(v) / sum(v)
  x <- c(min(d) - 1, d)          # anchor the curve at zero before the first day
  y <- c(0, cum)
  ok <- c(TRUE, diff(y) > 0)     # approx() needs a strictly increasing curve
  setNames(stats::approx(y[ok], x[ok], xout = probs, ties = "ordered")$y, names(out))
}
