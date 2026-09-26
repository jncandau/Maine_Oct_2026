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
#   "-"  no data for that date (trap not operated / not reported)
#   "m"  date within the trapping period but the value is missing
#   else a number, which may be fractional when the record is an average
#        (see `dataType`); Maine records are all dataType == "Counts".
#   These two codes were read off the data; confirm their exact definition in
#   data/raw/SBWFieldPhenoDatabase/esbwMetadata_v1.odt before publishing.
#
# COVERAGE of the Maine light-trap selection (checked 2026-09-26)
#   643 records, 79 localities, 30 years (1961-1989 and 2016)
#   535 records carry a daily series -> 10,230 counted trap-nights + 191 "m"
#   108 records carry only a season total (all day cells "-"): `has_daily` FALSE
#   season totals agree with the sum of the daily counts in every record with
#   no missing night; "Millinocket, Maine" 2016 has two records at the same
#   coordinates and reference (flightID 741 and 754), left unmerged
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
  filter(state %in% keep_states, method %in% keep_methods)

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
## those that only carry a season total in `moths_year_source`: 108 of the 643
## Maine light-trap records have every day cell set to "-".
trap_records <- trap_records %>%
  left_join(effort, by = "flightID") %>%
  mutate(has_daily = !is.na(n_nights))

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


## ---- 6. Write the processed tables -----------------------------------------

dir.create(file.path("data", "processed"), showWarnings = FALSE, recursive = TRUE)
saveRDS(trap_records, file.path("data", "processed", "trap_records.rds"))
saveRDS(trap_daily,   file.path("data", "processed", "trap_daily.rds"))

message(sprintf(
  paste0("Trap data read: %d records (%s, %s), %d trap-nights, ",
         "%d-%d.\n  total mismatches vs workbook: %d | duplicated site-years: %d"),
  nrow(trap_records), paste(keep_states, collapse = "/"),
  paste(keep_methods, collapse = "/"),
  sum(trap_daily$status == "counted"),
  min(trap_records$year), max(trap_records$year),
  nrow(total_mismatch), nrow(duplicate_records)
))
