# =============================================================================
# 09_regional_weather_test.R
# -----------------------------------------------------------------------------
# Does site-specific weather buy anything over one regional series?
#
# WHY THIS TEST
#   R/07 shows the catches carry almost no latitudinal gradient while both
#   development models impose about four days per degree. If flight really is
#   that synchronous across Maine, then a model driven by a single regional
#   weather series -- one prediction per year, the same for every locality --
#   should predict the observed catch dates about as well as the site-specific
#   runs. If it predicts them better, the site-specific weather is not
#   informing the models so much as misleading them.
#
# WHAT IS COMPARED
#   Four predictors of the observed median catch date, each with one fitted
#   intercept so that a constant bias is not charged against any of them:
#     A  BioSIM flight, site-specific            (from R/03, varies by site)
#     B  BioSIM flight, one regional point       (network centroid, by year)
#     C  bayessbw pupation, site-specific        (from R/04)
#     D  bayessbw pupation, one regional series  (mean hourly temperature)
#   B and D are constant within a year, so they carry no site information at
#   all. The intercept is fitted leave-one-year-out, so no predictor is
#   scored on a bias estimated from the year it is predicting.
#
#   A fifth line gives the floor: the mean observed date of the other
#   localities in the same year, which is what a perfect regional predictor
#   would look like.
#
# INPUT   data/processed/model_vs_observed.rds, weather_hourly.rds
# OUTPUT  data/processed/regional_weather_test.rds
#           list(regional_biosim, regional_bayes, scores, summary)
#
# REQUIREMENTS  dplyr, lubridate, sbwFieldPheno (see R/04), and network access
#               to the BioSIM Web API. No Java: the API is called directly.
#
# RUNTIME  a few minutes -- one dev_days() call and `n_reps` BioSIM requests
#          per year. `reuse` keeps both.
#
# USAGE  source("R/09_regional_weather_test.R")   # from the project root
# =============================================================================

library(dplyr)
library(lubridate)
library(sbwFieldPheno)


## ---- 1. Configuration -------------------------------------------------------

out_file   <- file.path("data", "processed", "regional_weather_test.rds")

## bayessbw settings, matching R/04 so the regional run differs from the
## site-specific one only in its weather.
colony        <- "NB"
window_months <- 3:8
n_post        <- 100
individuals   <- 100
seed          <- 1978

## BioSIM: replicate requests per year at the regional point. Run-to-run
## variation is about a tenth of a day, so a few replicates suffice.
n_reps       <- 5
biosim_url   <- "http://repicea.dynu.net:80/BioSIM/BioSimWeather"
biosim_parms <- "ApplyMortality:0*ApplyAdultMortality:0"

reuse <- TRUE


## ---- 2. Data ----------------------------------------------------------------

cmp <- readRDS(file.path("data", "processed", "model_vs_observed.rds")) %>%
  filter(comparable)

weather_hourly <- readRDS(file.path("data", "processed", "weather_hourly.rds"))

## The regional point: the centroid of every locality in the trap network. It
## does not depend on which localities were sampled in a given year.
centroid <- weather_hourly %>%
  distinct(Location, Latitude, Longitude) %>%
  summarise(Latitude = mean(Latitude), Longitude = mean(Longitude))

years <- sort(unique(cmp$Year))


## ---- 3. One regional temperature series per year ----------------------------
## The hourly series are complete calendar years on a common clock, so the
## regional series is the mean over localities at each hour. All localities
## present in the weather file are used, not only those trapped that year, so
## that the series describes the region rather than the sample.

regional_weather <- function(yr) {
  weather_hourly %>%
    filter(Year == yr) %>%
    select(Weather) %>%
    tidyr::unnest(Weather) %>%
    mutate(Date = as.Date(datetime_est), Hour = hour(datetime_est)) %>%
    filter(month(Date) %in% window_months) %>%
    group_by(Date, Hour) %>%
    summarise(Temp = mean(Temp), .groups = "drop") %>%
    arrange(Date, Hour)
}


## ---- 4. bayessbw on the regional series -------------------------------------

median_population <- function(m) apply(apply(m, 2, sort), 1, median)

if (reuse && file.exists(out_file)) {
  previous <- readRDS(out_file)
} else {
  previous <- NULL
}

if (!is.null(previous)) {
  regional_bayes <- previous$regional_bayes
} else {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  regional_bayes <- bind_rows(lapply(years, function(yr) {
    wx <- regional_weather(yr)
    m  <- dev_days(weather = wx, sbwcolony = colony, period = "hour",
                   stage = "Pupa", ecdf = FALSE,
                   n.post = n_post, individuals = individuals)
    tibble::tibble(Year = yr, hours = nrow(wx),
                   bayes_regional_p50 = stats::median(median_population(m)))
  }))
}


## ---- 5. BioSIM at the regional point ----------------------------------------

curve_p50 <- function(doy, value) {
  keep <- value > 0
  if (!any(keep)) return(NA_real_)
  cum <- cumsum(value[keep]) / sum(value[keep])
  stats::approx(c(0, cum), c(min(doy[keep]) - 1, doy[keep]),
                xout = 0.5, ties = "ordered")$y
}

biosim_flight_p50 <- function(lat, lon, year) {
  u <- sprintf(paste0("%s?lat=%f&long=%f&elev=NaN&from=%d&to=%d",
                      "&model=Spruce_Budworm_Biology&Parameters=%s"),
               biosim_url, lat, lon, year, year, biosim_parms)
  x <- try(utils::read.csv(url(u), skip = 1), silent = TRUE)
  if (inherits(x, "try-error")) return(NA_real_)
  doy <- as.integer(format(as.Date(sprintf("%d-%02d-%02d", x$Year, x$Month, x$Day)),
                           "%j"))
  curve_p50(doy, x$MaleFlight + x$FemaleFlight)
}

if (!is.null(previous)) {
  regional_biosim <- previous$regional_biosim
} else {
  regional_biosim <- bind_rows(lapply(years, function(yr) {
    reps <- vapply(seq_len(n_reps),
                   function(i) biosim_flight_p50(centroid$Latitude,
                                                 centroid$Longitude, yr),
                   numeric(1))
    tibble::tibble(Year = yr, n_reps = sum(!is.na(reps)),
                   biosim_regional_p50 = mean(reps, na.rm = TRUE))
  }))
}


## ---- 6. Score the predictors ------------------------------------------------
## Each predictor gets one intercept, fitted on the other years. The predictor
## itself is never refitted -- only its constant offset.

d <- cmp %>%
  select(Location, Year, Latitude, obs_p50,
         biosim_site = sim_p50, bayes_site = bayes_pupation_p50) %>%
  left_join(regional_biosim, by = "Year") %>%
  left_join(regional_bayes,  by = "Year") %>%
  rename(biosim_regional = biosim_regional_p50,
         bayes_regional  = bayes_regional_p50)

## The floor: the mean observed date of the OTHER localities in the same year.
d <- d %>%
  group_by(Year) %>%
  mutate(peers = (sum(obs_p50) - obs_p50) / (n() - 1)) %>%
  ungroup()

predictors <- c(`BioSIM, site-specific`      = "biosim_site",
                `BioSIM, one regional point` = "biosim_regional",
                `bayessbw, site-specific`    = "bayes_site",
                `bayessbw, one regional run` = "bayes_regional",
                `Other localities that year` = "peers")

loyo_pred <- function(v) {
  out <- rep(NA_real_, nrow(d))
  for (y in unique(d$Year)) {
    te <- d$Year == y
    tr <- !te & !is.na(v)
    if (!any(tr) || all(is.na(v[te]))) next
    out[te] <- v[te] + mean(d$obs_p50[tr] - v[tr])   # intercept from other years
  }
  out
}

scores <- bind_rows(lapply(names(predictors), function(nm) {
  v <- d[[predictors[[nm]]]]
  p <- loyo_pred(v)
  ok <- !is.na(p)
  tibble::tibble(
    predictor      = nm,
    varies_by_site = predictors[[nm]] %in% c("biosim_site", "bayes_site", "peers"),
    n              = sum(ok),
    rmse           = sqrt(mean((d$obs_p50[ok] - p[ok])^2)),
    mae            = mean(abs(d$obs_p50[ok] - p[ok])),
    r              = if (sd(v[ok]) > 0) cor(v[ok], d$obs_p50[ok]) else NA_real_
  )
}))


## ---- 7. Write and report ----------------------------------------------------

summary_list <- list(
  n_rows     = nrow(d),
  n_years    = length(years),
  centroid   = centroid,
  scores     = scores,
  ## Change in RMSE when the site-specific run is replaced by the regional
  ## one. Negative means the regional series predicts the catches better.
  rmse_change  = scores$rmse[scores$predictor == "BioSIM, one regional point"] -
                 scores$rmse[scores$predictor == "BioSIM, site-specific"],
  bayes_change = scores$rmse[scores$predictor == "bayessbw, one regional run"] -
                 scores$rmse[scores$predictor == "bayessbw, site-specific"]
)

dir.create(file.path("data", "processed"), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(regional_biosim = regional_biosim, regional_bayes = regional_bayes,
             rows = d, scores = scores, summary = summary_list),
        out_file)

message(sprintf(
  paste0("Regional against site-specific weather, %d locality-years, %d years.\n",
         "  regional point: %.2f N, %.2f W\n%s\n",
         "  replacing the site-specific run with the regional one changes RMSE by\n",
         "    %+.2f d (BioSIM) and %+.2f d (bayessbw); negative means the single\n",
         "    regional series predicts the catches BETTER than site-specific weather.\n",
         "  written to %s"),
  nrow(d), length(years), centroid$Latitude, abs(centroid$Longitude),
  paste(sprintf("    %-28s RMSE %5.2f d | MAE %5.2f d | r %s",
                scores$predictor, scores$rmse, scores$mae,
                ifelse(is.na(scores$r), "  --", sprintf("%.2f", scores$r))),
        collapse = "\n"),
  summary_list$rmse_change, summary_list$bayes_change, out_file
))
