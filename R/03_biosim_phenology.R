# =============================================================================
# 03_biosim_phenology.R
# -----------------------------------------------------------------------------
# Simulate spruce budworm seasonal biology with BioSIM at every trap
# locality-year, and reduce it to flight-phenology percentiles comparable with
# the ones computed from the trap catches.
#
# HOW IT WORKS
#   The BioSIM Web API is called through the BioSIM R package (which drives a
#   Java client through J4R). For each year, one request covers every locality
#   sampled that year; BioSIM generates the weather for those coordinates from
#   its own historical station network and runs the "Spruce_Budworm_Biology"
#   model on it, returning the daily number of individuals in each stage.
#
#   The model starts from a fixed cohort of second-instar larvae and reports
#   every stage as a PERCENTAGE of that cohort, quantised in steps of 0.25%.
#   The R package cannot change the cohort size, so the simulation is repeated
#   `n_runs` times per locality-year and the runs are averaged, which smooths
#   that quantisation -- the earlier light-trap analysis (Old LighttrapRhainds)
#   summed them instead, which differs only by the factor `n_runs` and so
#   changes no percentile. Averaging keeps the values on the percentage scale.
#   Each run draws a different seed.
#
#   Natural mortality is switched OFF (ApplyMortality = 0, ApplyAdultMortality
#   = 0). The purpose here is phenology -- WHEN the stages occur -- not
#   survival, and with mortality on, only a few percent of the cohort reaches
#   the adult stage, so the flight curve is estimated from a handful of
#   individuals. Checked against the Web API for one Maine locality-year
#   (46.608 N, -69.516 W, 1978): with mortality off the season sums are
#   MaleEmergence + FemaleEmergence = 100% and MaleFlight + FemaleFlight =
#   97%, against 10% and 6% under the defaults.
#
#   Male and female are summed for the three quantities of interest:
#     Pupae_mf  = MalePupae  + FemalePupae    (stock: pupae present that day)
#     Adult_mf  = MaleAdult  + FemaleAdult    (stock: adults present that day)
#     Flight_mf = MaleFlight + FemaleFlight   (activity: flight that day)
#   Only Flight is a daily flux, so its cumulative curve is a true
#   distribution; the two stocks are kept for comparison (see section 5).
#   MaleEmergence / FemaleEmergence are also in the output and are a flux too,
#   should an emergence-based phenology be wanted instead.
#
# INPUT   data/processed/trap_records.rds  (locality-years and coordinates)
# OUTPUT  data/interim/biosim_raw.rds       raw daily output, one row per run
#         data/processed/biosim_daily.rds   runs averaged, daily percentages
#         data/processed/biosim_phenology.rds  three rows per locality-year
#                                           (Pupae_mf, Adult_mf, Flight_mf)
#                                           with percentiles, peak, duration
#                                           and the `flight_complete` flag
#         data/processed/biosim_flight.rds  the Flight_mf rows only, one per
#                                           locality-year, for joining to the
#                                           trap percentiles, with the
#                                           `flight_complete` flag
#
# REQUIREMENTS
#   * Java 8 or later on the machine, and the J4R + BioSIM packages:
#       install.packages("J4R", repos = c("https://repicea.r-universe.dev",
#                                         "https://cloud.r-project.org"))
#       install.packages("BioSIM", repos = c("https://repicea.r-universe.dev",
#                                            "https://cloud.r-project.org"))
#     (see https://github.com/RNCan/BioSimClient_R/wiki)
#   * network access to the BioSIM Web API
#   * dplyr, tidyr, purrr
#
# RUNTIME  22 years x `n_runs` requests to a remote service. Expect minutes to
#          tens of minutes, and set `reuse_raw = TRUE` to avoid repeating the
#          calls once they have succeeded.
#
# USAGE  source("R/03_biosim_phenology.R")   # from the project root
# =============================================================================

library(dplyr)
library(tidyr)


## ---- 1. Configuration -------------------------------------------------------

raw_file   <- file.path("data", "interim",   "biosim_raw.rds")
daily_file <- file.path("data", "processed", "biosim_daily.rds")
pheno_file <- file.path("data", "processed", "biosim_phenology.rds")
flight_file <- file.path("data", "processed", "biosim_flight.rds")

## BioSIM model. getModelList() shows what else is available and
## getModelHelp("Spruce_Budworm_Biology") documents this one.
biosim_model <- "Spruce_Budworm_Biology"

## Replicate simulations per locality-year (see the note in the header).
n_runs <- 10

## Season kept for the phenology summary, as day of year. DOY 146-243 is
## 26 May - 31 August, the window used in the earlier analysis; it brackets
## every observed flight period in the trap data.
keep_doy <- 146:243

## Percentiles of the modelled flight curve, matching the trap-side summary.
probs <- c(0.05, 0.25, 0.50, 0.75, 0.95)

## Minimum share of the cohort that must reach flight for a locality-year's
## percentiles to describe the whole population rather than the subset that
## completed development. Below it, `flight_complete` is FALSE -- nothing is
## dropped, the analysis decides. With mortality off a complete locality-year
## reaches 95-99%; the shortfall is the part of the cohort that never finishes
## development (those locality-years end the year with every stage at zero and
## DeathAdult equal to their flight total).
##
## Which way this biases the percentiles is not clear a priori and is not
## resolved by the data: in the 1968-1989 run the 20 flagged locality-years
## have a LATER mean p50 (194.0) than the 207 complete ones (190.1), which is
## consistent with incompleteness arising in cold site-years where development
## runs late, rather than with the surviving fraction being the fast developers.
min_flight_pct <- 90

## Set to FALSE to force new calls to BioSIM even if the raw output exists.
reuse_raw <- TRUE

## Model parameters, passed to BioSIM through `additionalParms`. These are the
## server defaults except for the two mortality switches, which are off (see
## the header). getModelDefaultParameters("Spruce_Budworm_Biology") lists them;
## LastSampleDate is left out because its default is empty and it only applies
## when the simulation is anchored on a field observation of average instar.
biosim_parms <- c(
  "ApplyMortality"      = 0,   # default 1
  "ApplyAdultMortality" = 0,   # default 1
  "CreateEggs"          = 0,   # one generation only
  "TreeKind"            = 0,   # host tree
  "ObservedAI"          = 0,   # no observed average instar to start from
  "Defoliation"         = 0,
  "AdultLongevity"      = 0
)

## Columns used from the model output. Everything is a percentage of the
## simulated cohort. The sex-specific columns are the ones summed below;
## `Pupae` and `Adults` (both sexes together already) are kept as a cross-check.
sex_pairs <- list(
  Pupae_mf  = c("MalePupae",  "FemalePupae"),
  Adult_mf  = c("MaleAdult",  "FemaleAdult"),
  Flight_mf = c("MaleFlight", "FemaleFlight")
)
needed_cols <- c("Year", "Month", "Day", unlist(sex_pairs), "Pupae", "Adults")

## Series summarised into phenology percentiles, in the order they are reported.
pheno_series <- names(sex_pairs)


## ---- 2. Locality-years ------------------------------------------------------
## BioSIM returns the `id` it was given as `KeyID`. Locality names carry
## commas, dots and parentheses, so a short synthetic key is sent instead and
## the mapping is kept for the join back.

trap_records <- readRDS(file.path("data", "processed", "trap_records.rds"))

locality_years <- trap_records %>%
  distinct(Location = locality, Year = year,
           Latitude = latitude, Longitude = longitude) %>%
  arrange(Location, Year)

locs <- locality_years %>%
  distinct(Location, Latitude, Longitude) %>%
  mutate(KeyID = sprintf("L%02d", row_number()))

locality_years <- left_join(locality_years, select(locs, Location, KeyID),
                            by = "Location")

message(sprintf("BioSIM will simulate %d locality-years (%d localities, %d-%d) x %d runs",
                nrow(locality_years), nrow(locs),
                min(locality_years$Year), max(locality_years$Year), n_runs))


## ---- 3. Call BioSIM ---------------------------------------------------------
## One request per year and run: BioSIM is asked for the localities sampled in
## that year only, so nothing is simulated that the trap data does not cover.
## Elevation is left unspecified, which makes BioSIM take it from its digital
## elevation model.

#' Run the model for one year and one replicate
#'
#' @param year calendar year to simulate.
#' @param run  replicate index, kept in the output.
#' @return data.frame as returned by BioSIM, plus a `run` column, or NULL.
biosim_year <- function(year, run) {
  sel <- filter(locality_years, Year == year)

  out <- tryCatch(
    BioSIM::generateWeather(
      modelNames = biosim_model,
      fromYr     = year,
      toYr       = year,
      id         = sel$KeyID,
      latDeg     = sel$Latitude,
      longDeg    = sel$Longitude,
      ## One parameter map, applied to the single model requested.
      additionalParms = list(biosim_parms)
    )[[biosim_model]],
    error = function(e) {
      warning(sprintf("BioSIM failed for %d (run %d): %s", year, run,
                      conditionMessage(e)), call. = FALSE)
      NULL
    }
  )

  if (is.null(out)) return(NULL)
  mutate(out, run = run)
}

if (reuse_raw && file.exists(raw_file)) {
  biosim_raw <- readRDS(raw_file)
  message("Re-using ", raw_file, " (", nrow(biosim_raw), " rows)")
} else {
  years <- sort(unique(locality_years$Year))
  biosim_raw <- list()
  for (run in seq_len(n_runs)) {
    for (yr in years) {
      biosim_raw[[paste(yr, run)]] <- biosim_year(yr, run)
    }
    message(sprintf("  run %d of %d done", run, n_runs))
  }
  biosim_raw <- bind_rows(biosim_raw)

  ## The Java client keeps a connection open; close it when the calls are done.
  try(BioSIM::shutdownClient(), silent = TRUE)

  dir.create(file.path("data", "interim"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(biosim_raw, raw_file)
}

## What came back must contain the columns used below.
missing_cols <- setdiff(c("KeyID", needed_cols), names(biosim_raw))
if (length(missing_cols) > 0) {
  stop("BioSIM output is missing: ", paste(missing_cols, collapse = ", "),
       "\nCheck getModelHelp(\"", biosim_model, "\") -- the model output may have changed.",
       call. = FALSE)
}

## `Rep` indexes the climate replicates; this script asks for one, so more than
## one would mean the runs are being summed over something unintended.
if ("Rep" %in% names(biosim_raw) && n_distinct(biosim_raw$Rep) > 1) {
  stop("BioSIM returned several climate replicates (Rep); the summation in ",
       "section 4 assumes one.", call. = FALSE)
}

## BioSIM marks whether the weather behind a row was observed or generated.
if ("DataType" %in% names(biosim_raw)) {
  generated <- biosim_raw %>% filter(DataType != "Real_Data")
  if (nrow(generated) > 0) {
    warning(sprintf("%.1f%% of the simulated days rest on generated rather than observed weather (%s).",
                    100 * nrow(generated) / nrow(biosim_raw),
                    paste(unique(generated$DataType), collapse = ", ")), call. = FALSE)
  }
}

## Every requested locality-year must have come back, for every run.
returned <- biosim_raw %>% distinct(KeyID, Year, run) %>% count(KeyID, Year, name = "runs")
incomplete <- locality_years %>%
  left_join(returned, by = c("KeyID", "Year")) %>%
  mutate(runs = coalesce(runs, 0L)) %>%
  filter(runs < n_runs)
if (nrow(incomplete) > 0) {
  warning(nrow(incomplete), " locality-years have fewer than ", n_runs,
          " runs; their percentiles rest on a smaller cohort.", call. = FALSE)
}


## ---- 4. Average the runs and build the daily table --------------------------
## Every column is a percentage of the simulated cohort, quantised at 0.25%.
## Averaging `n_runs` independent runs smooths that step without changing the
## scale: the result is still "percent of the cohort on this day".

biosim_daily <- biosim_raw %>%
  mutate(date = as.Date(sprintf("%d-%02d-%02d", Year, Month, Day)),
         doy  = as.integer(format(date, "%j")),
         ## Male + female, for the three quantities of interest.
         Pupae_mf  = MalePupae  + FemalePupae,
         Adult_mf  = MaleAdult  + FemaleAdult,
         Flight_mf = MaleFlight + FemaleFlight) %>%
  group_by(KeyID, Year, doy, date) %>%
  summarise(across(all_of(c(pheno_series, "Pupae", "Adults", unlist(sex_pairs))),
                   ~ mean(.x, na.rm = TRUE)),
            n_runs = n_distinct(run), .groups = "drop") %>%
  left_join(select(locality_years, KeyID, Year, Location, Latitude, Longitude),
            by = c("KeyID", "Year")) %>%
  relocate(Location, Year, doy, date) %>%
  arrange(Location, Year, doy)

## The sex-specific columns should reproduce the both-sexes ones the model also
## reports; a persistent gap would mean the columns do not mean what is assumed.
sex_sum_check <- biosim_daily %>%
  summarise(pupae_gap = max(abs(Pupae_mf - Pupae)),
            adult_gap = max(abs(Adult_mf - Adults)))
if (with(sex_sum_check, pupae_gap > 1 || adult_gap > 1)) {
  warning(sprintf(paste("Male + female does not match the both-sexes column:",
                        "largest gap %.2f%% (pupae), %.2f%% (adults)."),
                  sex_sum_check$pupae_gap, sex_sum_check$adult_gap), call. = FALSE)
}


## ---- 5. Phenology percentiles -----------------------------------------------
## Percentiles are taken from the cumulative curve of each series, the same way
## the trap percentiles are computed, so the two are directly comparable.
##
## The three series are not the same kind of quantity:
##   Flight_mf is a daily activity rate -- a flux. Its cumulative curve is a
##     genuine distribution of flight over the season (it sums to ~97% of the
##     cohort with mortality off), and it is the series to compare with trap
##     catches.
##   Pupae_mf and Adult_mf are stocks: the percentage present on a day. Their
##     cumulative curves are weighted by how long individuals stay in the
##     stage, so with adult mortality off -- adults then persist for weeks --
##     the adult curve lags flight. They are reported for comparison, not as
##     an emergence date.

## curve_percentiles() is shared with the trap and comparison scripts.
source(file.path("R", "functions", "phenology.R"))

## Names of the percentile columns, e.g. p05 ... p95.
p_names <- sprintf("p%02d", round(probs * 100))

biosim_phenology <- biosim_daily %>%
  filter(doy %in% keep_doy) %>%
  select(Location, Year, KeyID, Latitude, Longitude, doy, n_runs,
         all_of(pheno_series)) %>%
  pivot_longer(all_of(pheno_series), names_to = "series", values_to = "value") %>%
  group_by(Location, Year, KeyID, Latitude, Longitude, series) %>%
  summarise(
    n_runs     = max(n_runs),
    season_pct = sum(value, na.rm = TRUE),
    peak_doy   = if (all(is.na(value))) NA_integer_ else doy[which.max(value)],
    peak_value = suppressWarnings(max(value, na.rm = TRUE)),
    tibble::as_tibble(as.list(curve_percentiles(doy, value, probs))),
    .groups    = "drop"
  ) %>%
  mutate(duration = .data[[p_names[length(p_names)]]] - .data[[p_names[1]]]) %>%
  relocate(series, .after = Year) %>%
  arrange(Location, Year, match(series, pheno_series))

## Wide view of the flight series alone, for joining to the trap percentiles.
## `flight_complete` marks the locality-years whose curve covers essentially
## the whole cohort; where it is FALSE the percentiles describe only the
## fraction that completed development (see `min_flight_pct` above for what
## that does, and does not, imply about their timing).
## Completeness is a property of the locality-year, not of one series, so the
## flag is attached to all three rows -- a pupal or adult curve from a cohort
## that never finished flying is just as partial.
flight_completeness <- biosim_phenology %>%
  filter(series == "Flight_mf") %>%
  transmute(Location, Year, flight_complete = season_pct >= min_flight_pct)

biosim_phenology <- biosim_phenology %>%
  left_join(flight_completeness, by = c("Location", "Year")) %>%
  relocate(flight_complete, .after = series)

biosim_flight <- biosim_phenology %>%
  filter(series == "Flight_mf") %>%
  select(-series)


## ---- 6. Write and report ----------------------------------------------------

dir.create(file.path("data", "processed"), showWarnings = FALSE, recursive = TRUE)
saveRDS(biosim_daily,     daily_file)
saveRDS(biosim_phenology, pheno_file)
saveRDS(biosim_flight,    flight_file)

flight <- filter(biosim_phenology, series == "Flight_mf")
message(sprintf(
  paste0("BioSIM phenology: %d locality-years, %d localities, %d-%d, %d runs each.\n",
         "  mortality off (ApplyMortality = %d, ApplyAdultMortality = %d)\n",
         "  modelled flight (Male+Female): p50 DOY %.0f to %.0f, mean %.1f; ",
         "season total %.0f%% to %.0f%% of the cohort\n",
         "  flight_complete (>= %d%% of the cohort): %d of %d locality-years\n",
         "  written to %s, %s and %s"),
  n_distinct(paste(flight$Location, flight$Year)), n_distinct(flight$Location),
  min(flight$Year), max(flight$Year), max(flight$n_runs),
  biosim_parms[["ApplyMortality"]], biosim_parms[["ApplyAdultMortality"]],
  min(flight$p50, na.rm = TRUE), max(flight$p50, na.rm = TRUE),
  mean(flight$p50, na.rm = TRUE),
  min(flight$season_pct, na.rm = TRUE), max(flight$season_pct, na.rm = TRUE),
  min_flight_pct, sum(biosim_flight$flight_complete), nrow(biosim_flight),
  daily_file, pheno_file, flight_file
))
