# =============================================================================
# 05_compare_model_observed.R
# -----------------------------------------------------------------------------
# Compare the modelled phenology with the observed light-trap catches, one row
# per locality-year.
#
# WHAT IS COMPARED
#   observed   percentiles of the nightly catch curve (R/01_read_trap_data.R)
#   BioSIM     percentiles of the daily flight curve, males + females, and of
#              the pupal curve (R/03_biosim_phenology.R)
#   bayessbw   percentiles of the simulated pupation dates
#              (R/04_pupal_phenology.R)
#
#   All three are percentiles of a curve over the season, computed by the same
#   helper, so that p50 always means "half of the season's total has happened".
#   For the traps the curve is moths per night; for BioSIM it is the percentage
#   of the cohort in flight (or pupating) per day; for bayessbw the "curve" is
#   the simulated population itself, already summarised in 04.
#
# WHAT IS NOT COMPARABLE, AND WHY THE FLAGS TRAVEL WITH THE NUMBERS
#   Two flags decide which rows carry a meaningful comparison, and both are
#   kept in the output rather than filtered out:
#     usable_phenology  the trapping season brackets the flight period
#                       (>= 15 nights, >= 70% coverage, <= 5% of the catch on
#                       the first or last counted night)
#     flight_complete   >= 90% of the modelled cohort reached flight; the
#                       others left part of the cohort undeveloped on 31 Dec
#   `comparable` is the conjunction. Headline offsets are computed on those
#   rows only, but every row is written out.
#
#   A caught moth is not a flying moth: the trap samples adults that are
#   active, close enough, and attracted to light, whereas BioSIM's flight
#   curve is the percentage of the cohort in the flight state. The comparison
#   is therefore about timing, not about numbers.
#
# INPUT   data/processed/trap_daily.rds, trap_records.rds
#         data/processed/biosim_flight.rds, biosim_phenology.rds
#         data/processed/pupal_phenology.rds
# OUTPUT  data/processed/model_vs_observed.rds  one row per locality-year
#
# USAGE  source("daily_dynamics/R/05_compare_model_observed.R")   # from the project root
# =============================================================================

library(dplyr)

source(file.path("R", "functions", "phenology.R"))


## ---- 1. Configuration -------------------------------------------------------

probs <- c(0.05, 0.25, 0.50, 0.75, 0.95)

out_file <- file.path("data", "processed", "model_vs_observed.rds")


## ---- 2. Observed percentiles from the nightly catches -----------------------
## Counted nights only: `status == "missing"` marks nights with no record, and
## curve_percentiles() must not read them as zero catches.

trap_daily   <- readRDS(file.path("data", "processed", "trap_daily.rds"))
trap_records <- readRDS(file.path("data", "processed", "trap_records.rds"))

counted <- filter(trap_daily, status == "counted")

obs_peak <- counted %>%
  group_by(Location = locality, Year = year) %>%
  summarise(obs_peak_doy = doy[which.max(moths)],
            obs_peak_n   = max(moths),
            .groups = "drop")

obs_pct <- counted %>%
  group_by(Location = locality, Year = year) %>%
  group_modify(~ tibble::as_tibble(as.list(curve_percentiles(.x$doy, .x$moths, probs)))) %>%
  ungroup() %>%
  rename_with(~ paste0("obs_", .x), matches("^p[0-9]{2}$"))

observed <- obs_peak %>%
  left_join(obs_pct, by = c("Location", "Year")) %>%
  mutate(obs_duration = obs_p95 - obs_p05)


## ---- 3. The modelled sides --------------------------------------------------

biosim_flight <- readRDS(file.path("data", "processed", "biosim_flight.rds")) %>%
  select(Location, Year, Latitude, Longitude, flight_complete,
         sim_peak_doy = peak_doy,
         sim_p05 = p05, sim_p25 = p25, sim_p50 = p50, sim_p75 = p75,
         sim_p95 = p95, sim_duration = duration)

biosim_pupae <- readRDS(file.path("data", "processed", "biosim_phenology.rds")) %>%
  ## The pupation-date distribution (R/03 `Pupation_mf`), not the median of
  ## the pupal stock, which falls about five days later (middle of the pupal
  ## period). Corrected 2026-10-01; earlier versions used `Pupae_mf`.
  filter(series == "Pupation_mf") %>%
  select(Location, Year, biosim_pupation_p50 = p50)

bayes_pupation <- readRDS(file.path("data", "processed", "pupal_phenology.rds")) %>%
  select(Location, Year, colony,
         bayes_pupation_p05 = p05, bayes_pupation_p50 = p50,
         bayes_pupation_p95 = p95)


## ---- 4. Join and take the offsets -------------------------------------------
## Offsets are model minus observed, in days: positive means the model is late.

model_vs_observed <- trap_records %>%
  select(flightID, Location = locality, Year = year, Latitude = latitude,
         Longitude = longitude, season_total, n_nights, coverage,
         usable_phenology) %>%
  left_join(observed, by = c("Location", "Year")) %>%
  left_join(select(biosim_flight, -Latitude, -Longitude), by = c("Location", "Year")) %>%
  left_join(biosim_pupae, by = c("Location", "Year")) %>%
  left_join(bayes_pupation, by = c("Location", "Year")) %>%
  mutate(
    comparable   = usable_phenology & flight_complete,
    d_p05        = sim_p05 - obs_p05,
    d_p50        = sim_p50 - obs_p50,
    d_p95        = sim_p95 - obs_p95,
    d_peak       = sim_peak_doy - obs_peak_doy,
    d_duration   = sim_duration - obs_duration,
    ## Days from predicted pupation to observed median catch. Pupa to adult
    ## plus pre-flight maturation, as the traps see it.
    lag_bayes    = obs_p50 - bayes_pupation_p50,
    lag_biosim   = obs_p50 - biosim_pupation_p50
  ) %>%
  arrange(Year, Location)

stopifnot(nrow(model_vs_observed) == nrow(trap_records),
          !any(is.na(model_vs_observed$obs_p50)),
          !any(is.na(model_vs_observed$sim_p50)))


## ---- 4b. bayessbw on the same axis as the catches ---------------------------
## bayessbw predicts pupation, the trap sees flight, so the two are not
## directly comparable. Adding the median observed pupation-to-catch interval
## puts them on one axis. The shift is CALIBRATED ON THESE DATA, over the
## comparable rows only: it removes the mean offset by construction, so
## `d_p50_bayes` tests whether bayessbw ranks sites and years like the traps
## do, not whether it predicts the right date.

lag_calibrated <- median(model_vs_observed$lag_bayes[model_vs_observed$comparable])

model_vs_observed <- model_vs_observed %>%
  mutate(lag_calibrated = lag_calibrated,
         bayes_pred_p50 = bayes_pupation_p50 + lag_calibrated,
         d_p50_bayes    = bayes_pred_p50 - obs_p50)


## ---- 5. Write and report ----------------------------------------------------

dir.create(file.path("data", "processed"), showWarnings = FALSE, recursive = TRUE)
saveRDS(model_vs_observed, out_file)

cmp <- filter(model_vs_observed, comparable)

message(sprintf(
  paste0("Compared %d locality-years, %d of them comparable ",
         "(usable trapping season and a complete modelled cohort).\n",
         "  BioSIM flight minus observed catch, median of the record:\n",
         "    p05 %+.1f d | p50 %+.1f d | p95 %+.1f d | peak %+.1f d | duration %+.1f d\n",
         "  p50 offset: mean %+.1f d, sd %.1f, range %+.1f to %+.1f\n",
         "  correlation of modelled and observed p50: r = %.2f (n = %d)\n",
         "  observed median catch follows predicted pupation by %.1f d ",
         "(bayessbw) and %.1f d (BioSIM)\n",
         "  bayessbw shifted by that lag: offset sd %.1f d, r = %.2f, and its\n",
         "    offset correlates with BioSIM's at %.2f\n",
         "  written to %s"),
  nrow(model_vs_observed), nrow(cmp),
  median(cmp$d_p05), median(cmp$d_p50), median(cmp$d_p95),
  median(cmp$d_peak), median(cmp$d_duration),
  mean(cmp$d_p50), sd(cmp$d_p50), min(cmp$d_p50), max(cmp$d_p50),
  cor(cmp$sim_p50, cmp$obs_p50), nrow(cmp),
  median(cmp$lag_bayes), median(cmp$lag_biosim),
  sd(cmp$d_p50_bayes), cor(cmp$bayes_pred_p50, cmp$obs_p50),
  cor(cmp$d_p50, cmp$d_p50_bayes),
  out_file
))
