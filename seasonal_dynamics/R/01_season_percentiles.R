# =============================================================================
# seasonal_dynamics/R/01_season_percentiles.R
# -----------------------------------------------------------------------------
# Reduce every Maine light-trap record to one row: the timing, spread and shape
# of its flight season. This is the response table of the Seasonal dynamics
# analysis, which models the season as a whole rather than night by night.
#
# INPUT   data/processed/trap_daily.rds    (R/01_read_trap_data.R)
#         data/processed/trap_records.rds  (R/01_read_trap_data.R)
#         Both come from the shared data layer, so the rules applied there --
#         each year's season window from its first to its last capture, "-"
#         nights inside it read as zeros, "m" nights left missing -- hold here
#         exactly as in the Daily dynamics analysis. Nothing in this script
#         reads an output of daily_dynamics/.
#
# WHAT IS COMPUTED, per locality-year
#   p05 ... p95   day of year by which that share of the season's catch was
#                 taken, from the cumulative catch over the counted nights
#                 (shared helper R/functions/phenology.R, so the convention is
#                 the one used everywhere in the project)
#   width90       p95 - p05, the central 90% of the season, in days
#   width50       p75 - p25, the central 50%
#   asymmetry     (p95 - p50) - (p50 - p05): positive when the decline is
#                 longer than the rise
#   peak_doy, peak_share   the night with the largest catch and its share
#   n_positive    nights with at least one moth: a season concentrated on very
#                 few nights has a well-defined median but an unreliable width
#
# CENSORING
#   "m" nights are not counted, so a percentile can be biased when a censored
#   night carried part of the season. After the season window is applied every
#   remaining censored night is right-censored (late tail). `n_cens_right` is
#   carried so that tail statistics can be checked against it; nothing is
#   imputed here.
#
# OUTPUT  data/processed/seasonal/season_percentiles.rds
#         one row per locality-year (226 rows)
#
# REQUIREMENTS  dplyr
#
# USAGE  source("seasonal_dynamics/R/01_season_percentiles.R")   # from the project root
# =============================================================================

library(dplyr)

source(file.path("R", "functions", "phenology.R"))


## ---- 1. Configuration ----------------------------------------------------------

probs   <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)
out_dir <- file.path("data", "processed", "seasonal")


## ---- 2. Inputs ------------------------------------------------------------------

trap_daily   <- readRDS(file.path("data", "processed", "trap_daily.rds"))
trap_records <- readRDS(file.path("data", "processed", "trap_records.rds"))


## ---- 3. Season statistics -------------------------------------------------------

season_stats <- function(doy, moths) {
  p <- curve_percentiles(doy, moths, probs)
  i <- which.max(moths)
  tibble::tibble(
    !!!as.list(p),
    width90     = unname(p["p95"] - p["p05"]),
    width50     = unname(p["p75"] - p["p25"]),
    asymmetry   = unname((p["p95"] - p["p50"]) - (p["p50"] - p["p05"])),
    peak_doy    = doy[i],
    peak_share  = moths[i] / sum(moths),
    n_positive  = sum(moths > 0)
  )
}

season <- trap_daily %>%
  filter(status == "counted") %>%
  arrange(flightID, doy) %>%
  group_by(flightID) %>%
  reframe(season_stats(doy, moths)) %>%
  left_join(
    trap_records %>%
      select(flightID, locality, latitude, longitude, year,
             season_total, season_first, season_last,
             n_nights, n_zero_filled, n_missing, n_cens_right,
             usable_phenology),
    by = "flightID"
  ) %>%
  relocate(flightID, locality, latitude, longitude, year) %>%
  arrange(year, locality)

stopifnot(nrow(season) == nrow(trap_records),
          !any(is.na(season$p50)))


## ---- 4. Write and report ------------------------------------------------------

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(season, file.path(out_dir, "season_percentiles.rds"))

cat(sprintf(
  paste0("Season percentiles: %d locality-years, %d localities, %d years ",
         "(%d usable for phenology).\n",
         "  median catch date: median DOY %.1f (range %.0f-%.0f)\n",
         "  central 90%%: median %.1f d | central 50%%: median %.1f d | ",
         "asymmetry: median %+.1f d\n",
         "  peak night carries a median %.0f%% of the season; ",
         "%d record(s) with fewer than 3 positive nights\n",
         "  written to %s\n"),
  nrow(season), n_distinct(season$locality), n_distinct(season$year),
  sum(season$usable_phenology),
  median(season$p50), min(season$p50), max(season$p50),
  median(season$width90), median(season$width50), median(season$asymmetry),
  100 * median(season$peak_share), sum(season$n_positive < 3),
  file.path(out_dir, "season_percentiles.rds")))
