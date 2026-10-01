# =============================================================================
# 08_elevation_check.R
# -----------------------------------------------------------------------------
# Is BioSIM's spatial error an elevation error?
#
# WHY THIS TEST
#   R/07 shows that both development models impose a latitudinal gradient of
#   about four days per degree on flight timing while the catches carry almost
#   none. Elevation is the obvious suspect: R/03 leaves it unspecified, so
#   BioSIM takes it from its own digital elevation model, and elevation is
#   itself a latitude proxy in this network -- the southern localities are
#   coastal lowlands, the northern ones are 300 m and above.
#
# WHAT IS DONE
#   1. an elevation for every locality from the USGS 3DEP point-query service;
#   2. a sensitivity scan: one locality-year run over a range of elevations,
#      to establish how much elevation can move the modelled flight date;
#   3. the test itself: modelled flight with BioSIM's own elevation model
#      against modelled flight with the USGS elevation, for a sample of
#      locality-years spanning the latitude range;
#   4. a check on how precisely the location has to be right: during this work
#      a coordinate 1.5 km from the trap changed the modelled date by several
#      days at a coastal locality, so the displacement is scanned here.
#
#   The Web API is called directly (http://repicea.dynu.net/BioSIM/), not
#   through the R client, so this script needs no Java. It runs a sample
#   rather than the whole dataset: a full re-run with elevation is one switch
#   away in R/03 (`use_elevation`) and the result below is what decides
#   whether it is worth the 2,270 requests.
#
# INPUT   data/processed/trap_records.rds, model_vs_observed.rds
# OUTPUT  data/processed/locality_elevation.rds  Location, elev_m
#         data/processed/elevation_check.rds     list(preview, scan, coords)
#
# USAGE  source("daily_dynamics/R/08_elevation_check.R")   # from the project root
# =============================================================================

library(dplyr)


## ---- 1. Configuration -------------------------------------------------------

elev_file  <- file.path("data", "processed", "locality_elevation.rds")
check_file <- file.path("data", "processed", "elevation_check.rds")

## Localities sampled for the NaN-against-USGS comparison, spread evenly over
## the latitude range; each is taken in its largest-catch comparable year.
n_sample <- 20

## Elevations used for the sensitivity scan, in metres.
scan_elev <- c(0, 50, 100, 150, 200, 250, 300, 350, 400, 500, 600, 800, 1000)

## Localities for the scan: one southern, one northern.
scan_sites <- c("Blue Hill Township, Maine", "Clayton Lake, Maine")

biosim_url <- "http://repicea.dynu.net:80/BioSIM/BioSimWeather"
biosim_parms <- "ApplyMortality:0*ApplyAdultMortality:0"

reuse <- TRUE


## ---- 2. Elevation from the USGS 3DEP point service --------------------------

fetch_elevation <- function(lat, lon, tries = 4, pause = 2) {
  u <- sprintf(paste0("https://epqs.nationalmap.gov/v1/json?x=%.6f&y=%.6f",
                      "&units=Meters&wkid=4326&includeDate=false"), lon, lat)
  for (k in seq_len(tries)) {
    r <- try(suppressWarnings(paste(readLines(u, warn = FALSE), collapse = "")),
             silent = TRUE)
    if (!inherits(r, "try-error") && grepl('"value"', r)) {
      v <- suppressWarnings(as.numeric(sub('.*"value":"?([-0-9.]+)"?.*', "\\1", r)))
      if (!is.na(v)) return(v)
    }
    Sys.sleep(pause)
  }
  NA_real_
}

trap_records <- readRDS(file.path("data", "processed", "trap_records.rds"))

locs <- trap_records %>%
  distinct(Location = locality, Latitude = latitude, Longitude = longitude) %>%
  arrange(Location)

if (reuse && file.exists(elev_file)) {
  elevation <- readRDS(elev_file)
} else {
  elevation <- locs %>%
    rowwise() %>%
    mutate(elev_m = fetch_elevation(Latitude, Longitude)) %>%
    ungroup() %>%
    select(Location, elev_m)
  dir.create(file.path("data", "processed"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(elevation, elev_file)
}

locs <- left_join(locs, elevation, by = "Location")
stopifnot(!any(is.na(locs$elev_m)))


## ---- 3. One BioSIM call ------------------------------------------------------
## `elev` is passed through verbatim, so "NaN" asks BioSIM for its own
## elevation model and a number overrides it.

curve_p50 <- function(doy, value) {
  keep <- value > 0
  if (!any(keep)) return(NA_real_)
  cum <- cumsum(value[keep]) / sum(value[keep])
  stats::approx(c(0, cum), c(min(doy[keep]) - 1, doy[keep]),
                xout = 0.5, ties = "ordered")$y
}

biosim_flight_p50 <- function(lat, lon, year, elev) {
  u <- sprintf(paste0("%s?lat=%f&long=%f&elev=%s&from=%d&to=%d",
                      "&model=Spruce_Budworm_Biology&Parameters=%s"),
               biosim_url, lat, lon, elev, year, year, biosim_parms)
  x <- try(utils::read.csv(url(u), skip = 1), silent = TRUE)
  if (inherits(x, "try-error")) return(c(p50 = NA_real_, total = NA_real_))
  doy <- as.integer(format(as.Date(sprintf("%d-%02d-%02d", x$Year, x$Month, x$Day)),
                           "%j"))
  flight <- x$MaleFlight + x$FemaleFlight
  c(p50 = curve_p50(doy, flight), total = sum(flight))
}


## ---- 4. Sensitivity: how much can elevation move the date? ------------------

cmp <- readRDS(file.path("data", "processed", "model_vs_observed.rds")) %>%
  filter(comparable)

scan_years <- cmp %>%
  filter(Location %in% scan_sites) %>%
  group_by(Location) %>%
  slice_max(season_total, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Location, Year, Latitude, Longitude)

scan <- scan_years %>%
  tidyr::expand_grid(elev_m = scan_elev) %>%
  rowwise() %>%
  mutate(p50 = biosim_flight_p50(Latitude, Longitude, Year,
                                 sprintf("%.0f", elev_m))[["p50"]]) %>%
  ungroup()


## ---- 5. The test: BioSIM's elevation model against the USGS elevation ------

sample_ly <- cmp %>%
  group_by(Location) %>%
  slice_max(season_total, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(elevation, by = "Location") %>%
  arrange(Latitude) %>%
  slice(round(seq(1, n(), length.out = n_sample)))

preview <- sample_ly %>%
  rowwise() %>%
  mutate(p50_dem  = biosim_flight_p50(Latitude, Longitude, Year, "NaN")[["p50"]],
         p50_usgs = biosim_flight_p50(Latitude, Longitude, Year,
                                      sprintf("%.0f", elev_m))[["p50"]]) %>%
  ungroup() %>%
  transmute(Location, Year, Latitude, elev_m, obs_p50,
            p50_dem, p50_usgs,
            shift    = p50_usgs - p50_dem,
            off_dem  = p50_dem  - obs_p50,
            off_usgs = p50_usgs - obs_p50)


## ---- 6. How precisely does the location have to be right? -------------------
## The coastal locality of the scan, displaced east by a series of steps, with
## BioSIM's own elevation model each time. Rounding a coordinate to two
## decimals is a displacement of under a kilometre at this latitude; the
## question is at what displacement the modelled date starts to move.

coord_site <- filter(scan_years, Location == scan_sites[1])

coord_steps <- c(0, 0.005, 0.01, 0.02, 0.05, 0.1)   # degrees of longitude

coords <- tibble::tibble(d_lon = coord_steps) %>%
  rowwise() %>%
  mutate(
    km  = d_lon * 111.32 * cos(coord_site$Latitude * pi / 180),
    p50 = biosim_flight_p50(coord_site$Latitude, coord_site$Longitude + d_lon,
                            coord_site$Year, "NaN")[["p50"]]
  ) %>%
  ungroup() %>%
  mutate(shift = p50 - p50[d_lon == 0])


## ---- 7. Write and report ----------------------------------------------------

summary_list <- list(
  n_localities   = nrow(elevation),
  elev_range     = range(elevation$elev_m),
  elev_median    = median(elevation$elev_m),
  elev_lat_r     = cor(locs$Latitude, locs$elev_m),
  elev_south     = median(locs$elev_m[locs$Latitude <  45.5]),
  elev_north     = median(locs$elev_m[locs$Latitude >= 45.5]),
  n_preview      = sum(!is.na(preview$shift)),
  mean_shift     = mean(preview$shift, na.rm = TRUE),
  shift_south    = mean(preview$shift[preview$Latitude <  45.5], na.rm = TRUE),
  shift_north    = mean(preview$shift[preview$Latitude >= 45.5], na.rm = TRUE),
  abs_off_dem    = mean(abs(preview$off_dem),  na.rm = TRUE),
  abs_off_usgs   = mean(abs(preview$off_usgs), na.rm = TRUE),
  grad_dem       = coef(lm(p50_dem  ~ Latitude, preview))[[2]],
  grad_usgs      = coef(lm(p50_usgs ~ Latitude, preview))[[2]],
  grad_obs       = coef(lm(obs_p50  ~ Latitude, preview))[[2]],
  ## days per 100 m over the plausible range, from the scan
  sens_per_100m  = scan %>%
    filter(elev_m >= 100, elev_m <= 500) %>%
    group_by(Location) %>%
    summarise(d = 100 * coef(lm(p50 ~ elev_m))[[2]], .groups = "drop")
)

dir.create(file.path("data", "processed"), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(preview = preview, scan = scan, coords = coords,
             elevation = elevation, summary = summary_list),
        check_file)

message(sprintf(
  paste0("Elevation check.\n",
         "  %d localities: %.0f-%.0f m (median %.0f); elevation vs latitude r = %+.2f\n",
         "    southern median %.0f m, northern %.0f m\n",
         "  sensitivity between 100 and 500 m: %s\n",
         "  %d locality-years, BioSIM's elevation model against the USGS value:\n",
         "    mean shift %+.1f d (south %+.1f, north %+.1f)\n",
         "    mean absolute offset %.1f d -> %.1f d\n",
         "    latitude gradient %+.2f -> %+.2f d per degree (observed %+.2f)\n",
         "  written to %s"),
  summary_list$n_localities, summary_list$elev_range[1], summary_list$elev_range[2],
  summary_list$elev_median, summary_list$elev_lat_r,
  summary_list$elev_south, summary_list$elev_north,
  paste(sprintf("%s %+.1f d/100 m", sub(",.*", "", summary_list$sens_per_100m$Location),
                summary_list$sens_per_100m$d), collapse = "; "),
  summary_list$n_preview, summary_list$mean_shift,
  summary_list$shift_south, summary_list$shift_north,
  summary_list$abs_off_dem, summary_list$abs_off_usgs,
  summary_list$grad_dem, summary_list$grad_usgs, summary_list$grad_obs,
  check_file
))
