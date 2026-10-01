# =============================================================================
# 04_pupal_phenology.R
# -----------------------------------------------------------------------------
# Predict the date of pupation at every trap locality-year with the Bayesian
# development model of Studens et al. (https://github.com/kdis19/bayessbw),
# driven by the hourly CaSR temperatures extracted in R/02_read_weather_data.R.
#
# WHAT THE MODEL IS
#   bayessbw fits a Sharpe-DeMichele / Schoolfield development-rate curve to
#   laboratory rearing data for the five larval stages L2-L6, in a Bayesian
#   hierarchical form: each posterior draw carries a rate curve per stage plus
#   a lognormal individual-variation term (`s_eps`). The fitted curves are
#   distributed as the R package sbwFieldPheno, whose dev_days() accumulates
#   development hour by hour through a temperature series and returns the day
#   of year on which each simulated individual reaches a given stage.
#
#   Pupation is the end of L6, so `stage = "Pupa"` is the last transition the
#   larval model describes; the pupal stage itself is not modelled in bayessbw
#   (its laboratory data stop at L6). "Pupal phenology" here therefore means
#   the timing of pupation, not the duration of the pupal stage.
#
#   Two sources of variation are simulated and kept apart:
#     individuals  within-population variation (the `s_eps` multiplier)
#     posterior    parameter uncertainty (one rate curve set per draw)
#   With n_post = 100 every posterior draw in the package is used, so the
#   posterior spread is not itself a sample.
#
# THE START OF THE SERIES IS A MODELLING CHOICE
#   dev_days() starts accumulating development on the first row of the weather
#   it is given, i.e. the first row acts as the date of post-diapause L2
#   emergence. The window below (1 March to 31 August) is the one used in
#   bayessbw's own weather simulations (five_pops/code/new_weather.R).
#
# INPUT   data/processed/weather_hourly.rds  (hourly Temp per locality-year)
# OUTPUT  data/interim/pupal_raw.rds          individuals x draws per locality-year
#         data/processed/pupal_phenology.rds  one row per locality-year
#         data/processed/pupal_population.rds the median simulated population,
#                                             one row per individual rank
#
# REQUIREMENTS  dplyr, tidyr, tibble, lubridate, parallel, and
#   remotes::install_github("studensk/sbwFieldPheno")
#
# RUNTIME  about 8 s per locality-year on one core (measured, 100 x 100), so
#          roughly 30 min for the 227 locality-years single-threaded; the jobs
#          are run in parallel over `n_cores`.
#
# USAGE  source("daily_dynamics/R/04_pupal_phenology.R")   # from the project root
# =============================================================================

library(dplyr)
library(lubridate)
## Attached, not used through `::`: dev_days() looks up the fitted curves
## (`curves`, `curve_info`) unqualified, and they are lazy-loaded package data
## that are only on the search path once the package is attached.
library(sbwFieldPheno)


## ---- 1. Configuration -------------------------------------------------------

raw_file  <- file.path("data", "interim",   "pupal_raw.rds")
pheno_file <- file.path("data", "processed", "pupal_phenology.rds")
pop_file   <- file.path("data", "processed", "pupal_population.rds")

## Source colony whose fitted curves are used. sbwFieldPheno carries six:
## AB, IPQL, NB, NWT, ON, QC. New Brunswick is the closest to Maine of the
## field colonies (IPQL is a long-standing laboratory colony). Add more to the
## vector to compare them; the colony is carried in every output table.
colonies <- "NB"

## Development window. The first hour is the assumed start of post-diapause
## L2 development, so this choice shifts every predicted date.
window_months <- 3:8

## Simulation size. n_post is capped at 100 by the package (and 100 is all the
## draws it holds); individuals is the simulated population per draw.
n_post      <- 100
individuals <- 100

## Percentiles of the simulated population, matching the trap and BioSIM sides.
probs <- c(0.05, 0.25, 0.50, 0.75, 0.95)

seed <- 1978

## detectCores() is documented to return NA on some platforms, and
## max(1, NA) is NA, which mclapply() rejects -- so guard it.
n_cores <- getOption("mc.cores", local({
  d <- parallel::detectCores()
  if (is.na(d)) 1L else max(1L, as.integer(d) - 2L)
}))

## Set to FALSE to re-simulate even if the raw output exists.
reuse_raw <- TRUE



## ---- 2. Weather, one frame per locality-year --------------------------------
## dev_days() wants Temp (degC), Date and Hour. Local standard time is used so
## that a "day" is a Maine day; the hourly series are complete, which matters
## because dev_days() warns on gaps and accumulates over whatever it is given.

weather_hourly <- readRDS(file.path("data", "processed", "weather_hourly.rds"))

weather_for <- function(nested) {
  nested %>%
    transmute(Temp = Temp,
              Date = as.Date(datetime_est),
              Hour = hour(datetime_est)) %>%
    filter(month(Date) %in% window_months) %>%
    arrange(Date, Hour)
}

jobs <- weather_hourly %>%
  select(Location, Year, Latitude, Longitude, Weather) %>%
  rowwise() %>%
  mutate(weather = list(weather_for(Weather))) %>%
  ungroup() %>%
  select(-Weather) %>%
  tidyr::expand_grid(colony = colonies)

expected_hours <- vapply(jobs$weather, nrow, integer(1))
if (any(expected_hours < 24 * 150)) {
  warning(sum(expected_hours < 24 * 150), " locality-years have a short weather window.",
          call. = FALSE)
}

message(sprintf("Simulating pupation for %d locality-years x %d colonies (%d draws x %d individuals)",
                n_distinct(paste(jobs$Location, jobs$Year)), length(colonies),
                n_post, individuals))


## ---- 3. Simulate ------------------------------------------------------------
## One dev_days() call per locality-year and colony. The call returns a matrix
## of Julian dates with `individuals` rows and `n_post` columns.

simulate_one <- function(i) {
  dev_days(
    weather     = jobs$weather[[i]],
    sbwcolony   = jobs$colony[i],
    period      = "hour",
    stage       = "Pupa",
    ecdf        = FALSE,
    n.post      = n_post,
    individuals = individuals
  )
}

if (reuse_raw && file.exists(raw_file)) {
  pupal_raw <- readRDS(raw_file)
  message("Re-using ", raw_file)
} else {
  ## L'Ecuyer streams so that the parallel runs are reproducible.
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  pupal_raw <- parallel::mclapply(seq_len(nrow(jobs)), simulate_one,
                                  mc.cores = n_cores, mc.set.seed = TRUE)

  failed <- !vapply(pupal_raw, is.matrix, logical(1))
  if (any(failed)) {
    stop(sum(failed), " simulations failed, e.g.: ",
         as.character(pupal_raw[[which(failed)[1]]]), call. = FALSE)
  }

  dir.create(file.path("data", "interim"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(list(jobs = select(jobs, -weather), matrices = pupal_raw), raw_file)
}

if (is.list(pupal_raw) && !is.null(pupal_raw$matrices)) {
  jobs_done <- pupal_raw$jobs
  pupal_mat <- pupal_raw$matrices
} else {
  jobs_done <- select(jobs, -weather)
  pupal_mat <- pupal_raw
}
stopifnot(nrow(jobs_done) == length(pupal_mat))


## ---- 4. The median simulated population ------------------------------------
## Each column of a matrix is one posterior draw: a simulated population of
## `individuals` pupation dates. Sorting a column turns it into that draw's
## population distribution; taking the median across draws, rank by rank, gives
## the median population distribution -- the summary used in bayessbw.

median_population <- function(m) apply(apply(m, 2, sort), 1, median)

pupal_population <- lapply(seq_along(pupal_mat), function(i) {
  tibble::tibble(
    Location = jobs_done$Location[i],
    Year     = jobs_done$Year[i],
    colony   = jobs_done$colony[i],
    rank     = seq_len(nrow(pupal_mat[[i]])),
    pupation_jday = median_population(pupal_mat[[i]])
  )
}) %>%
  bind_rows()


## ---- 5. Percentiles and uncertainty -----------------------------------------
## Population percentiles come from the median population above. The posterior
## interval is computed the other way round: the median pupation date of each
## draw, then its 2.5th and 97.5th percentiles across draws -- so the two
## columns separate individual variation from parameter uncertainty.

pupal_phenology <- lapply(seq_along(pupal_mat), function(i) {
  m   <- pupal_mat[[i]]
  pop <- median_population(m)
  draw_medians <- apply(m, 2, median)
  pct <- stats::quantile(pop, probs = probs, names = FALSE, type = 7)

  tibble::tibble(
    Location    = jobs_done$Location[i],
    Year        = jobs_done$Year[i],
    colony      = jobs_done$colony[i],
    n_post      = ncol(m),
    individuals = nrow(m),
    !!!setNames(as.list(pct), sprintf("p%02d", round(probs * 100))),
    spread      = pct[length(pct)] - pct[1],
    post_median = median(draw_medians),
    post_lo     = stats::quantile(draw_medians, 0.025, names = FALSE),
    post_hi     = stats::quantile(draw_medians, 0.975, names = FALSE)
  )
}) %>%
  bind_rows() %>%
  left_join(select(weather_hourly, Location, Year, Latitude, Longitude),
            by = c("Location", "Year")) %>%
  relocate(Latitude, Longitude, .after = colony) %>%
  arrange(colony, Location, Year)


## ---- 6. Write and report ----------------------------------------------------

dir.create(file.path("data", "processed"), showWarnings = FALSE, recursive = TRUE)
saveRDS(pupal_phenology,  pheno_file)
saveRDS(pupal_population, pop_file)

message(sprintf(
  paste0("Pupation predicted for %d locality-years (colony %s, window %s, ",
         "%d draws x %d individuals).\n",
         "  population median (p50): DOY %.0f to %.0f, mean %.1f\n",
         "  within-population spread (p05-p95): %.1f d on average\n",
         "  posterior interval on the median: %.1f d wide on average\n",
         "  written to %s and %s"),
  nrow(pupal_phenology), paste(colonies, collapse = "/"),
  paste(range(window_months), collapse = "-"), n_post, individuals,
  min(pupal_phenology$p50), max(pupal_phenology$p50), mean(pupal_phenology$p50),
  mean(pupal_phenology$spread),
  mean(pupal_phenology$post_hi - pupal_phenology$post_lo),
  pheno_file, pop_file
))
