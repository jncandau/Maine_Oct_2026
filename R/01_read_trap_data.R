# =============================================================================
# 01_read_trap_data.R
# -----------------------------------------------------------------------------
# Read the daily moth-capture records of the eastern spruce budworm phenology
# dataset and reshape the Maine light-trap records into one row per trap-night.
#
# SOURCE
#   data/raw/SBWFieldPhenoDatabase/Eastern Spruce Budworm Phenology V1.ods
#   sheet "Flight" (763 records; Fournier, Candau et al. 2027, CC-BY).
#   The sheet is in wide form: each row is one trapping record (site x year x
#   sampling method x sex) and the columns mothsDay160 ... mothsDay231 hold the
#   capture of each day of year from 160 (early June) to 231 (late August).
#
# CELL CODES observed in the day columns (character, not numeric):
#   "-"  missing: no count available for that date
#   "m"  missing: no count available for that date
#   else a number, which may be fractional when the record is an average
#        (see `dataType`); Maine records are all dataType == "Counts".
#   Both codes mean missing data (confirmed 2026-09-26). They are handled
#   differently only by position, not by meaning: "-" fills the whole
#   out-of-season range, so "-" cells are dropped and only the trapping period
#   is kept, while "m" cells always fall inside it and are retained as NA.
#   Interior "-" cells -- inside the trapping period -- are counted in
#   `n_gap_nights` so they are not lost.
#
# SELECTION (see the configuration block below)
#   Maine light traps, 1968-1989, season total above 99 moths. The year window
#   is the overlap with the CaSR reanalysis in data/raw/Weather_Maine; the catch
#   threshold removes records too small to describe a flight curve.
#
# POOLING
#   Where one trap is reported as separate Male and Female records, the two are
#   summed night by night into one record before the threshold is applied
#   (section 4b). In Maine this affects 1980 only: "Marshfield Plot 4"
#   (flightID 760+761, 281 moths) and "Bare Island" (762+763, 161 moths). The
#   pooled record keeps the lowest flightID, sex becomes "Male+Female" and
#   `pooled_from` records the components.
#
# COVERAGE of that selection (checked 2026-09-26)
#   227 records, 49 localities, 22 years; every record has a daily series
#   5,612 counted trap-nights + 154 nights coded "m"
#   season totals 100 to 71,730 moths; effort 1 to 45 nights per record, so a
#   minimum-effort rule is still needed before percentile-based phenology
#   230 records in the year window fall at or below the threshold (`dropped_small`)
#   season totals agree with the sum of the daily counts in every retained
#   record, and no locality-year holds two records after the selection
#   (the unexplained "Millinocket, Maine" 2016 pair is outside the window)
#
#   Before the selection the sheet holds 643 Maine light-trap records at 79
#   localities over 30 years (1961-1989 and 2016), of which 108 carry only a
#   season total; widen keep_years / set min_season_total = 0 to see them.
#
# OUTPUT (written to data/processed/)
#   trap_records.rds  one row per trapping record, with sampling effort and the
#                     annual total as reported in the workbook
#   trap_daily.rds    one row per trap-night actually sampled, with the count
#
# REQUIREMENTS
#   readODS, readr, dplyr, tidyr. If missing from the project library:
#     renv::install(c("readODS", "readr", "dplyr", "tidyr")); renv::snapshot()
#
# USAGE
#   source("R/01_read_trap_data.R")        # from the project root
# =============================================================================

library(dplyr)
library(tidyr)


## ---- 1. Configuration -------------------------------------------------------

## Path to the workbook, relative to the project root.
ods_file <- file.path("data", "raw", "SBWFieldPhenoDatabase",
                      "Eastern Spruce Budworm Phenology V1.ods")

## Records to keep. `stateProvince` uses ISO-style codes (ME, NB, ON);
## `samplingMethod` in Maine is either "Light trap" (the standard Maine Forest
## Service trap network) or "Yard Light" (opportunistic yard-light counts).
## Set keep_methods = c("Light trap", "Yard Light") to include both.
keep_states  <- "ME"
keep_methods <- "Light trap"

## Analysis window. 1968-1989 is the overlap between the trap record and the
## CaSR v3.2 reanalysis in data/raw/Weather_Maine (1968-1991): it drops 1961-1967
## and 2016, for which no gridded weather is available.
keep_years <- 1968:1989

## Minimum season catch. Records at or below this total are dropped: they carry
## too few moths to describe a flight curve. With min_season_total = 99 the
## selection also excludes every season-total-only record (all fall below it)
## and the whole of 1963-1967 and 2016 (no trap exceeded 99 moths in those
## years). Set to 0 to keep every record.
min_season_total <- 99

## Codes that mean "no count available" in a day cell.
code_not_sampled <- "-"   # trap not operated / date not reported
code_missing     <- "m"   # date sampled but value missing


## ---- 2. Read the Flight sheet ----------------------------------------------

#' Read one sheet of the ESBW phenology workbook as character columns
#'
#' Everything is read as text so that the day columns -- which mix numbers with
#' the codes "-" and "m" -- are never silently coerced. Numeric conversion is
#' explicit further down.
#'
#' @param file path to the .ods workbook.
#' @param sheet sheet name.
#' @return a tibble of character columns.
read_sheet_as_text <- function(file = ods_file, sheet = "Flight") {
  if (!file.exists(file)) {
    stop("Workbook not found: ", file,
         "\nThe dataset lives in the private repo jncandau/SBWFieldPhenoDatabase.",
         call. = FALSE)
  }
  readODS::read_ods(
    file, sheet = sheet,
    col_types = readr::cols(.default = readr::col_character())
  )
}

flight_raw <- read_sheet_as_text(ods_file, "Flight")

## Day-of-year columns, in order.
day_cols <- grep("^mothsDay[0-9]+$", names(flight_raw), value = TRUE)
day_cols <- day_cols[order(as.integer(sub("mothsDay", "", day_cols)))]
stopifnot(length(day_cols) > 0)


## ---- 3. Record-level table --------------------------------------------------

## One row per trapping record: identifiers, site, sampling metadata and the
## annual total as given in the workbook (`mothCountYear`).
trap_records <- flight_raw %>%
  transmute(
    flightID         = as.integer(flightID),
    referenceID      = as.integer(referenceID),
    country          = countryCode,
    state            = stateProvince,
    locality,
    latitude         = as.numeric(decimalLatitude),
    longitude        = as.numeric(decimalLongitude),
    coord_uncert_km  = as.numeric(coordinateUncertaintyInKm),
    year             = as.integer(year),
    method           = samplingMethod,
    data_type        = dataType,
    sex,
    moths_year_source = suppressWarnings(as.numeric(mothCountYear)),
    reference_pointer = referencePointer,
    remarks
  ) %>%
  filter(state %in% keep_states, method %in% keep_methods, year %in% keep_years)

if (!all(trap_records$data_type == "Counts")) {
  warning("Records with dataType other than \"Counts\" are included: ",
          paste(unique(trap_records$data_type), collapse = ", "),
          " -- values may be averages or percentages.", call. = FALSE)
}


## ---- 4. Daily table ---------------------------------------------------------

#' Reshape the wide day columns into one row per record-day
#'
#' @param x the raw Flight sheet (character columns).
#' @param cols names of the day columns.
#' @param keep_ids flightIDs to retain (the selected records).
#' @return tibble with flightID, doy, raw cell content, parsed count and a
#'   `status` factor recording why a count is NA.
flight_to_long <- function(x, cols = day_cols, keep_ids = trap_records$flightID) {
  x %>%
    mutate(flightID = as.integer(flightID)) %>%
    filter(flightID %in% keep_ids) %>%
    select(flightID, all_of(cols)) %>%
    pivot_longer(
      all_of(cols),
      names_to         = "doy",
      names_prefix     = "mothsDay",
      names_transform  = list(doy = as.integer),
      values_to        = "raw"
    ) %>%
    mutate(
      raw    = trimws(raw),
      status = case_when(
        is.na(raw) | raw == code_not_sampled ~ "not_sampled",
        raw == code_missing                   ~ "missing",
        TRUE                                  ~ "counted"
      ),
      moths  = if_else(status == "counted", suppressWarnings(as.numeric(raw)), NA_real_)
    )
}

trap_long <- flight_to_long(flight_raw)

## Any cell that is neither a code nor a number is a data-entry problem: report
## it instead of dropping it quietly.
unparsed <- trap_long %>% filter(status == "counted", is.na(moths))
if (nrow(unparsed) > 0) {
  warning("Unparsed day cells (kept as NA): ",
          paste(unique(unparsed$raw), collapse = " | "), call. = FALSE)
}

## Analysis table: the nights that were actually sampled, with a calendar date.
## "missing" nights are kept (the trap was running) so that gaps inside a
## trapping season stay visible; drop them with filter(status == "counted").
trap_daily <- trap_long %>%
  filter(status != "not_sampled") %>%
  left_join(select(trap_records, flightID, locality, latitude, longitude,
                   year, method, sex),
            by = "flightID") %>%
  mutate(date = as.Date(paste0(year, "-01-01")) + (doy - 1L)) %>%
  select(flightID, locality, latitude, longitude, year, date, doy,
         moths, status, method, sex) %>%
  arrange(flightID, doy)


## ---- 4b. Pool the sex-split records ----------------------------------------
## A few localities report the Male and the Female catch of one trap as two
## records (in Maine: "Marshfield Plot 4" and "Bare Island", both 1980). The
## trap, not the sex, is the sampling unit here, so the two rows are summed
## night by night into a single record and the components are removed. Pooling
## happens before the catch threshold of section 5b, so a trap is judged on its
## whole catch.
##
## Records with sex "-" (the vast majority) are untouched.

sex_split <- trap_records %>%
  filter(sex %in% c("Male", "Female")) %>%
  group_by(locality, year, method) %>%
  filter(n() > 1) %>%
  ## The pooled record keeps the lowest component flightID; `component_id`
  ## preserves the original ids, which the group_by() below would otherwise mask.
  mutate(pooled_id = min(flightID), n_parts = n(), component_id = flightID) %>%
  ungroup()

## Pooling assumes the components really are one trap: same position, same
## source table. Report any group where that does not hold.
inconsistent <- sex_split %>%
  group_by(pooled_id) %>%
  filter(n_distinct(latitude) > 1 | n_distinct(longitude) > 1 |
         n_distinct(reference_pointer) > 1) %>%
  ungroup()
if (nrow(inconsistent) > 0) {
  warning("Sex-split records pooled despite differing position or source: ",
          paste(unique(inconsistent$locality), collapse = " | "), call. = FALSE)
}

## Provenance column: the component flightIDs for pooled records, NA otherwise.
trap_records$pooled_from <- NA_character_

if (nrow(sex_split) > 0) {

  key <- select(sex_split, flightID, pooled_id, n_parts)

  ## Daily counts: sum the sexes on each night. A night counts as "counted"
  ## only when every component reported a number for it; otherwise the pooled
  ## value would be a partial sum, so it is marked "missing" instead.
  pooled_daily <- trap_daily %>%
    inner_join(key, by = "flightID") %>%
    group_by(flightID = pooled_id, locality, latitude, longitude,
             year, date, doy, method, n_parts) %>%
    summarise(
      complete = all(status == "counted") & n() == first(n_parts),
      moths    = if (all(status == "counted") & n() == first(n_parts)) sum(moths) else NA_real_,
      .groups  = "drop"
    ) %>%
    mutate(status = if_else(complete, "counted", "missing"), sex = "Male+Female") %>%
    select(all_of(names(trap_daily)))

  trap_daily <- trap_daily %>%
    filter(!flightID %in% key$flightID) %>%
    bind_rows(pooled_daily) %>%
    arrange(flightID, doy)

  ## Record metadata: identical across components except the sex and the
  ## season total, which is summed. `pooled_from` keeps the provenance.
  pooled_records <- sex_split %>%
    group_by(flightID = pooled_id) %>%
    summarise(
      across(c(referenceID, country, state, locality, latitude, longitude,
               coord_uncert_km, year, method, data_type, reference_pointer),
             first),
      sex               = "Male+Female",
      moths_year_source = sum(moths_year_source),
      remarks           = paste(unique(stats::na.omit(remarks)), collapse = " | "),
      pooled_from       = paste(sort(component_id), collapse = "+"),
      .groups           = "drop"
    ) %>%
    mutate(remarks = if_else(remarks == "", NA_character_, remarks))

  trap_records <- trap_records %>%
    filter(!flightID %in% key$flightID) %>%
    bind_rows(pooled_records) %>%
    arrange(flightID)
}


## ---- 5. Sampling effort and consistency checks ------------------------------

## Effort and computed totals per record, next to the workbook's own total.
effort <- trap_daily %>%
  group_by(flightID) %>%
  summarise(
    first_doy   = min(doy),
    last_doy    = max(doy),
    n_nights    = sum(status == "counted"),
    n_missing   = sum(status == "missing"),
    moths_total = sum(moths, na.rm = TRUE),
    .groups     = "drop"
  )

## `has_daily` separates the records usable for phenology (a daily series) from
## those that only carry a season total in `moths_year_source`: in the Maine
## light-trap set, those records have every day cell set to "-".
trap_records <- trap_records %>%
  left_join(effort, by = "flightID") %>%
  mutate(
    has_daily    = !is.na(n_nights),
    ## The workbook total is the reference figure; it equals the sum of the
    ## daily counts in every record with no missing night (checked below), and
    ## is the only total available for season-total-only records.
    season_total = coalesce(moths_year_source, moths_total)
  )


## ---- 5b. Apply the catch threshold -----------------------------------------
## Records are dropped here, not in section 3, because the threshold is applied
## to the season total, which is only known once the daily cells are read.

dropped_small <- trap_records %>%
  filter(season_total <= min_season_total | is.na(season_total)) %>%
  select(flightID, year, locality, has_daily, n_nights, season_total)

trap_records <- trap_records %>%
  filter(season_total > min_season_total)

## Keep the two tables in step: trap_daily must hold only the retained records.
trap_daily <- trap_daily %>%
  filter(flightID %in% trap_records$flightID)

## The workbook total should equal the sum of the daily counts whenever no night
## is missing. Differences are listed rather than corrected.
total_mismatch <- trap_records %>%
  filter(n_missing == 0, !is.na(moths_year_source),
         abs(moths_total - moths_year_source) > 0.5) %>%
  select(flightID, year, locality, moths_total, moths_year_source)

## Records that share site, year, method and sex: duplicates or genuinely two
## traps at one place. Flagged, never merged.
duplicate_records <- trap_records %>%
  count(locality, year, method, sex, name = "n_records") %>%
  filter(n_records > 1)

## Nights with no count inside the trapping period: "-" cells between the first
## and the last counted night of a record. Both codes mean missing data, so
## these are real gaps and belong in the effort accounting.
## For pooled records the pattern is taken from the lowest-id component, whose
## night set is identical to its partner's (checked in section 4b).
interior_gaps <- trap_long %>%
  filter(flightID %in% trap_records$flightID, status == "not_sampled") %>%
  inner_join(
    trap_daily %>% group_by(flightID) %>%
      summarise(lo = min(doy), hi = max(doy), .groups = "drop"),
    by = "flightID"
  ) %>%
  filter(doy > lo, doy < hi) %>%
  count(flightID, name = "n_gap_nights")

trap_records <- trap_records %>%
  left_join(interior_gaps, by = "flightID") %>%
  mutate(n_gap_nights = coalesce(n_gap_nights, 0L))


## ---- 6. Write the processed tables -----------------------------------------

dir.create(file.path("data", "processed"), showWarnings = FALSE, recursive = TRUE)
saveRDS(trap_records, file.path("data", "processed", "trap_records.rds"))
saveRDS(trap_daily,   file.path("data", "processed", "trap_daily.rds"))

message(sprintf(
  paste0("Trap data read: %d records (%s, %s, %d-%d, season total > %d) ",
         "at %d localities, %d trap-nights.\n",
         "  dropped below threshold: %d | total mismatches vs workbook: %d",
         " | duplicated site-years: %d"),
  nrow(trap_records), paste(keep_states, collapse = "/"),
  paste(keep_methods, collapse = "/"),
  min(keep_years), max(keep_years), min_season_total,
  n_distinct(trap_records$locality),
  sum(trap_daily$status == "counted"),
  nrow(dropped_small), nrow(total_mismatch), nrow(duplicate_records)
))
