# =============================================================================
# 02_read_weather_data.R
# -----------------------------------------------------------------------------
# Build an hourly weather series for every trap locality-year of the working
# dataset, from the ECCC Canadian Surface Reanalysis (CaSR v3.2).
#
# Adapted from 2026_Maine_Weather.R (LightTrapsRhainds project). What changed:
#   * the locality-years come from data/processed/trap_records.rds, so the
#     weather and the captures cannot drift apart;
#   * no setwd() and no download block -- the tiles are already in data/raw/
#     (the download code is kept, commented, at the end for reference);
#   * each locality is assigned to the grid tile that contains it, and the
#     tiles are extracted separately, instead of merging two 35,064-layer
#     rasters for every variable and period;
#   * the raw 3-hourly air temperature is kept next to the interpolated
#     hourly series, so the interpolation can always be audited.
#
# SOURCE  data/raw/Weather_Maine/CaSR_v3.2_<var>_<tile>_<period>.nc
#         5 variables x 2 tiles x 6 four-year periods covering 1968-1991 on a
#         ~10 km rotated-pole grid, hourly except air temperature:
#           A_PR0_SFC  analysis  quantity of precipitation  m
#           A_TT_1.5m  analysis  air temperature at 1.5 m   degC  (3-hourly)
#           P_FB_SFC   forecast  downward solar flux        W/m2
#           P_P0_SFC   forecast  surface pressure           mb
#           P_UVC_10m  forecast  wind modulus at 10 m       kts
#         Downloaded from
#         https://hpfx.collab.science.gc.ca/~scar700/rcas-casr/download_CaSR_regions_var_period_fr.html
#         see data/raw/Weather_Maine/How_I_got_the_data.png.
#
# OUTPUT  data/processed/weather_hourly.rds -- one row per locality-year with a
#         nested `Weather` tibble of 8,760 (8,784 in leap years) hourly rows:
#           datetime     UTC, as stored in the CaSR files
#           datetime_est local standard time (UTC-5, no daylight saving)
#           Temp         air temperature, interpolated to hourly (degC)
#           TT           air temperature as stored, 3-hourly, NA elsewhere
#           PR0, FB, P0, UVC   the four hourly variables
#
# REQUIREMENTS  dplyr, tidyr, tibble, lubridate, terra, chillR
#   renv::install(c("terra", "chillR")); renv::snapshot()
#
# RUNTIME  reads ~4 GB of NetCDF; expect tens of minutes.
#
# USAGE  source("R/02_read_weather_data.R")   # from the project root
# =============================================================================

library(dplyr)
library(tidyr)
library(lubridate)


## ---- 1. Configuration -------------------------------------------------------

weather_dir <- file.path("data", "raw", "Weather_Maine")
out_file    <- file.path("data", "processed", "weather_hourly.rds")

## Local standard time offset used for `datetime_est`. Maine is UTC-5 in
## winter and UTC-4 in summer; a constant offset is used on purpose, so that
## "hours after sunset" means the same thing in every record. Set to 0 to keep
## UTC only.
tz_offset_hours <- -5

## The two tiles of the CaSR grid that cover Maine.
tiles <- c("rlon561-595_rlat351-385", "rlon561-595_rlat386-420")

## Each file holds four calendar years.
periods <- c("1968-1971", "1972-1975", "1976-1979",
             "1980-1983", "1984-1987", "1988-1991")

var_info <- tibble::tribble(
  ~variable, ~filecode,   ~units,
  "PR0",     "A_PR0_SFC", "m",
  "TT",      "A_TT_1.5m", "degC",
  "FB",      "P_FB_SFC",  "W/m2",
  "P0",      "P_P0_SFC",  "mb",
  "UVC",     "P_UVC_10m", "kts"
)

#' Period file that contains a given calendar year
period_for_year <- function(year) {
  idx <- pmin(pmax((year - 1968) %/% 4, 0), length(periods) - 1)
  if_else(year < 1968 | year > 1991, NA_character_, periods[idx + 1])
}

#' Path of one CaSR file
nc_path <- function(variable, tile, period) {
  file.path(weather_dir,
            paste0("CaSR_v3.2_", var_info$filecode[var_info$variable == variable],
                   "_", tile, "_", period, ".nc"))
}

#' Files that together hold a whole calendar year
#'
#' A period file does not start and end at midnight: it runs from 31 December
#' 13:00 UTC of the year before the period to 31 December 12:00 UTC of its last
#' year. The last 11 hours of a period's final year therefore sit in the next
#' period file, and a year that ends a period needs both files.
files_for_year <- function(variable, tile, year) {
  pers  <- unique(c(period_for_year(year), period_for_year(year + 1L)))
  pers  <- pers[!is.na(pers)]
  paths <- vapply(pers, nc_path, character(1), variable = variable, tile = tile)
  unname(paths)
}


## ---- 2. Locality-years to extract ------------------------------------------
## One row per trapping record; the weather is the same for two records at the
## same place and year, so the extraction runs on distinct locality-years.

trap_records <- readRDS(file.path("data", "processed", "trap_records.rds"))

locality_years <- trap_records %>%
  distinct(Location = locality, Year = year,
           Latitude = latitude, Longitude = longitude) %>%
  mutate(period = period_for_year(Year)) %>%
  arrange(Location, Year)

## A locality must have one position, or the point extraction is ambiguous.
stopifnot(!anyDuplicated(distinct(locality_years, Location, Latitude, Longitude)$Location))

outside <- filter(locality_years, is.na(period))
if (nrow(outside) > 0) {
  warning(nrow(outside), " locality-years fall outside 1968-1991 and are skipped: ",
          paste(unique(outside$Year), collapse = ", "), call. = FALSE)
  locality_years <- filter(locality_years, !is.na(period))
}

## Every file that will be opened must exist before the long loop starts.
needed <- expand_grid(variable = var_info$variable, tile = tiles,
                      Year = sort(unique(locality_years$Year))) %>%
  rowwise() %>%
  mutate(path = list(files_for_year(variable, tile, Year))) %>%
  ungroup() %>%
  tidyr::unnest(path)
absent <- unique(needed$path[!file.exists(needed$path)])
if (length(absent) > 0) {
  stop("Missing CaSR file(s):\n  ", paste(basename(absent), collapse = "\n  "),
       call. = FALSE)
}


## ---- 3. Locate the trap sites on the grid ----------------------------------
## The grid is a rotated pole, so the coordinates are projected once. Each
## locality is then assigned to the tile that actually has data at its
## position; the tiles overlap slightly, and the first match is used.

locs <- distinct(locality_years, Location, Latitude, Longitude)

ref_rast <- terra::rast(nc_path("PR0", tiles[1], locality_years$period[1]), lyrs = 1)
pts <- terra::vect(as.data.frame(locs), geom = c("Longitude", "Latitude"),
                   crs = "EPSG:4326") %>%
  terra::project(terra::crs(ref_rast))

in_tile <- vapply(tiles, function(tl) {
  r <- terra::rast(nc_path("PR0", tl, locality_years$period[1]), lyrs = 1)
  !is.na(terra::extract(r, pts, ID = FALSE)[[1]])
}, logical(nrow(locs)))

locs$tile <- tiles[apply(in_tile, 1, function(x) which(x)[1])]
if (anyNA(locs$tile)) {
  stop("Locality outside both CaSR tiles: ",
       paste(locs$Location[is.na(locs$tile)], collapse = " | "), call. = FALSE)
}
message(sprintf("%d localities on the grid (%s)",
                nrow(locs),
                paste(sprintf("%s: %d", basename(tiles), table(locs$tile)[tiles]),
                      collapse = ", ")))


## ---- 4. Extract the hourly series ------------------------------------------
## One file is opened at a time, and only the layers of the years that are
## actually needed are read, for the localities that belong to that tile.

#' Extract one variable, one tile, one calendar year
#'
#' @return long tibble: Location, datetime (UTC), value
extract_year <- function(path, year, pts_sub) {
  r    <- terra::rast(path)
  ts   <- terra::time(r)
  keep <- which(year(ts) == year)
  if (length(keep) == 0) return(NULL)
  m <- as.matrix(terra::extract(r[[keep]], pts_sub, ID = FALSE))  # points x layers
  tibble::tibble(
    Location = rep(pts_sub$Location, times = length(keep)),
    datetime = rep(ts[keep], each = nrow(m)),
    value    = as.vector(m)
  )
}

weather_list <- list()

for (yr in sort(unique(locality_years$Year))) {
  for (v in var_info$variable) {
    for (tl in tiles) {
      wanted <- locality_years %>%
        filter(Year == yr, Location %in% locs$Location[locs$tile == tl]) %>%
        pull(Location)
      if (length(wanted) == 0) next

      pts_sub <- pts[pts$Location %in% wanted]
      ## One or two files, depending on whether `yr` ends a period.
      out <- lapply(files_for_year(v, tl, yr), extract_year, year = yr, pts_sub = pts_sub) %>%
        bind_rows() %>%
        distinct(Location, datetime, .keep_all = TRUE)

      if (nrow(out) == 0) next
      weather_list[[paste(yr, v, tl)]] <- mutate(out, variable = v, Year = yr)
    }
    gc()
  }
  message(sprintf("  %d extracted (%d locality-years)", yr,
                  sum(locality_years$Year == yr)))
}

weather_long <- bind_rows(weather_list)
rm(weather_list); gc()


## ---- 5. One row per locality-year, with the hourly table nested ------------

weather_wide <- weather_long %>%
  pivot_wider(id_cols = c(Location, Year, datetime),
              names_from = variable, values_from = value) %>%
  mutate(datetime_est = datetime + tz_offset_hours * 3600) %>%
  arrange(Location, Year, datetime)

## Sanity: an hourly series of a full calendar year, and temperature on every
## third hour only.
hours_check <- weather_wide %>%
  group_by(Location, Year) %>%
  summarise(n_hours = n(), n_tt = sum(!is.na(TT)), .groups = "drop") %>%
  mutate(expected = if_else(leap_year(Year), 8784L, 8760L))
if (any(hours_check$n_hours != hours_check$expected)) {
  warning("Locality-years with an incomplete hourly series: ",
          sum(hours_check$n_hours != hours_check$expected), call. = FALSE)
}


## ---- 6. Interpolate air temperature to hourly ------------------------------
## Air temperature is stored every three hours. chillR::interpolate_gaps_hourly
## fits the idealised daily temperature curve of Linvill (1990) to the daily
## minima and maxima and interpolates the residuals, which is why it needs the
## latitude. The filled series is joined back on the calendar fields rather
## than by position, so a change in the number of rows returned cannot go
## unnoticed.

interpolate_temperature <- function(w, latitude) {
  h <- w %>%
    transmute(Year = year(datetime), Month = month(datetime),
              Day = day(datetime), Hour = hour(datetime), Temp = TT)
  filled <- chillR::interpolate_gaps_hourly(h, latitude = latitude)$weather %>%
    transmute(Year, Month, Day, Hour, Temp)
  w %>%
    mutate(Year = year(datetime), Month = month(datetime),
           Day = day(datetime), Hour = hour(datetime)) %>%
    left_join(filled, by = c("Year", "Month", "Day", "Hour")) %>%
    select(datetime, datetime_est, Temp, TT, PR0, FB, P0, UVC)
}

weather_hourly <- weather_wide %>%
  nest(.by = c(Location, Year), .key = "Weather") %>%
  left_join(locs, by = "Location") %>%
  rowwise() %>%
  mutate(Weather = list(interpolate_temperature(Weather, Latitude))) %>%
  ungroup() %>%
  select(Location, Year, Latitude, Longitude, tile, Weather)


## ---- 7. Write and report ----------------------------------------------------

dir.create(file.path("data", "processed"), showWarnings = FALSE, recursive = TRUE)
saveRDS(weather_hourly, out_file)

na_report <- weather_hourly %>%
  rowwise() %>%
  mutate(n_hours = nrow(Weather), n_na_temp = sum(is.na(Weather$Temp))) %>%
  ungroup() %>%
  summarise(locality_years = n(), hours = sum(n_hours),
            hours_no_temp = sum(n_na_temp))

message(sprintf(
  "Weather written to %s: %d locality-years, %s hourly rows, %s without temperature.",
  out_file, na_report$locality_years,
  format(na_report$hours, big.mark = ","),
  format(na_report$hours_no_temp, big.mark = ",")
))


## ---- Reference: how the tiles were downloaded -------------------------------
## The file list was generated with the CaSR download page (link at the top of
## this script) and saved as data/raw/Weather_Maine/RCaSv3p2_liste-fichier.txt.
##
## files <- trimws(readLines(file.path(weather_dir, "RCaSv3p2_liste-fichier.txt")))
## for (f in files) {
##   curl::curl_download(
##     url      = paste0("https://hpfx.collab.science.gc.ca/~scar700/rcas-casr/",
##                       "data/CaSRv3.2/netcdf_tile/", f),
##     destfile = file.path(weather_dir, sub("^.*/", "", f)),
##     quiet    = FALSE, mode = "wb"
##   )
## }
