# =============================================================================
# 14_censoring.R
# -----------------------------------------------------------------------------
# Quantify what the missing ("m") nights cost, now that they are read as
# censoring rather than as nights to be ignored.
#
# THE RULE (applied in R/01_read_trap_data.R)
#   A missing night does not carry an unknown catch drawn from the whole
#   season: it censors one tail of the flight curve. Nights before 15 July
#   bound the rise of the season (left-censored), nights after it bound the
#   decline (right-censored). `trap_daily$censoring` records the side and
#   `trap_records$n_cens_left` / `n_cens_right` count them per record.
#
# WHAT THIS SCRIPT ADDS
#   The side alone does not say whether a gap matters. A block of twelve
#   censored nights in mid-June, before anything is flying, costs nothing,
#   while two censored nights at the peak can hide a third of the season. So
#   for every affected record the fitted predictor (R/13) is used to ask how
#   much of the season is expected inside the censored nights, and the
#   percentiles are then re-estimated with those nights completed at their
#   expected level. The difference between the two estimates is the bias that
#   treating a censored night as absent introduces.
#
#   This is a diagnostic, not a correction: nothing downstream is changed. It
#   identifies the records whose percentiles are not safe to use at face
#   value, and shows that the published aggregates do not move.
#
# INPUT   data/processed/trap_daily.rds, trap_records.rds
#         data/processed/flight_timing.rds      (regional emergence curves)
#         data/processed/weather_hourly.rds     (hourly temperature)
#         data/processed/model_vs_observed.rds  (flags and observed curves)
#         data/processed/flight_predictor.rds   (via predict_flight())
# OUTPUT  data/processed/censoring_check.rds -- list(rows, headline, summary)
#
# REQUIREMENTS  dplyr, splines, glmmTMB (loaded by predict_flight)
#
# RUNTIME  a few seconds.
# =============================================================================

library(dplyr)
library(splines)

source(file.path("R", "functions", "phenology.R"))
source(file.path("R", "functions", "predict_flight.R"))


## ---- 1. Inputs ---------------------------------------------------------------

trap_daily   <- readRDS(file.path("data", "processed", "trap_daily.rds"))
trap_records <- readRDS(file.path("data", "processed", "trap_records.rds"))
timing       <- readRDS(file.path("data", "processed", "flight_timing.rds"))
weather      <- readRDS(file.path("data", "processed", "weather_hourly.rds"))
cmp          <- readRDS(file.path("data", "processed", "model_vs_observed.rds"))

## The driver used for the completion. Emergence is the series available at
## prediction time and scores within 0.15 d of the flight series (R/12).
driver   <- "emergence"
kappa    <- "width"
probs    <- c(0.05, 0.50, 0.95)
material <- 0.01   # censored share above which a record is called material


## ---- 2. Which records are affected ------------------------------------------

affected <- trap_records %>%
  filter(n_cens_left + n_cens_right > 0) %>%
  transmute(flightID, locality, year, season_total,
            n_cens_left, n_cens_right,
            side = case_when(n_cens_right == 0 ~ "left",
                             n_cens_left  == 0 ~ "right",
                             TRUE              ~ "left+right"))


## ---- 3. Expected mass inside the censored nights -----------------------------
## The predictor is given the record's own night set -- counted AND censored --
## so the shares it returns are conditional on the same window the trap
## covered. `cens_share` is then the expected fraction of the season that fell
## on nights with no count.

assess_one <- function(id) {
  nights <- trap_daily %>% filter(flightID == id) %>% arrange(doy)
  loc <- nights$locality[1]
  yr  <- nights$year[1]

  curve <- timing$regional_emergence %>% filter(Year == yr) %>%
    transmute(doy, pct = emergence)
  wx <- weather %>% filter(Location == loc, Year == yr) %>% pull(Weather)

  ## A year with no comparable locality has no regional curve, and a record
  ## outside the weather set has no temperature: report, do not guess.
  if (nrow(curve) == 0 || length(wx) == 0) {
    return(tibble(flightID = id, assessed = FALSE))
  }

  wx <- wx[[1]] %>% transmute(datetime = datetime_est, Temp)
  pred <- predict_flight(curve, wx, sum(nights$moths, na.rm = TRUE),
                         driver = driver, nights_doy = nights$doy,
                         kappa = kappa)

  nights <- left_join(nights, select(pred, doy, share), by = "doy")
  cens   <- nights$status == "missing"

  ## Complete the censored nights at their expected level: the observed total
  ## is rescaled onto the shares of the counted nights, and the censored
  ## nights get their own share at that rate.
  obs_total <- sum(nights$moths, na.rm = TRUE)
  rate      <- obs_total / sum(nights$share[!cens])
  completed <- ifelse(cens, nights$share * rate, nights$moths)

  p_obs  <- curve_percentiles(nights$doy[!cens], nights$moths[!cens], probs)
  p_full <- curve_percentiles(nights$doy, completed, probs)

  tibble(
    flightID      = id,
    assessed      = TRUE,
    cens_share    = sum(nights$share[cens]),
    implied_extra = sum(completed[cens]),
    cens_doy_lo   = min(nights$doy[cens]),
    cens_doy_hi   = max(nights$doy[cens]),
    d_p05         = p_full[1] - p_obs[1],
    d_p50         = p_full[2] - p_obs[2],
    d_p95         = p_full[3] - p_obs[3]
  )
}

rows <- affected %>%
  left_join(bind_rows(lapply(affected$flightID, assess_one)), by = "flightID") %>%
  left_join(select(cmp, flightID, usable_phenology, flight_complete, comparable),
            by = "flightID") %>%
  mutate(is_material = assessed & cens_share > material) %>%
  arrange(desc(cens_share))


## ---- 4. Do the published aggregates move? ------------------------------------
## The two headline spatial numbers are recomputed with the completed
## percentiles substituted for the affected records.

head_rows <- cmp %>%
  filter(comparable) %>%
  left_join(select(filter(rows, assessed), flightID, d_p50c = d_p50), by = "flightID") %>%
  mutate(d_p50c = coalesce(d_p50c, 0),
         obs_p50_completed = obs_p50 + d_p50c,
         ## d_p50 is sim - obs, so completing the observed curve moves the
         ## offset the other way.
         offset_completed  = d_p50 - d_p50c)

grad <- function(y) unname(coef(lm(y ~ Latitude + factor(Year), data = head_rows))["Latitude"])

headline <- tibble(
  quantity  = c("median offset, d", "observed gradient, d per degree N"),
  published = c(median(head_rows$d_p50), grad(head_rows$obs_p50)),
  completed = c(median(head_rows$offset_completed), grad(head_rows$obs_p50_completed))
)


## ---- 5. Summary and output ---------------------------------------------------

summary_list <- list(
  n_cens_nights   = sum(trap_daily$status == "missing"),
  n_left          = sum(trap_daily$censoring == "left",  na.rm = TRUE),
  n_right         = sum(trap_daily$censoring == "right", na.rm = TRUE),
  n_records       = nrow(rows),
  n_assessed      = sum(rows$assessed),
  n_comparable    = sum(rows$comparable, na.rm = TRUE),
  n_material      = sum(rows$is_material, na.rm = TRUE),
  median_share    = median(rows$cens_share[rows$assessed]),
  max_share       = max(rows$cens_share[rows$assessed]),
  max_shift_p50   = max(abs(rows$d_p50[rows$assessed])),
  max_shift_p95   = max(abs(rows$d_p95[rows$assessed])),
  driver          = driver,
  kappa           = kappa,
  material_cut    = material,
  headline        = headline
)

saveRDS(list(rows = rows, headline = headline, summary = summary_list),
        file.path("data", "processed", "censoring_check.rds"))

cat(sprintf(
  paste0("Censoring: %d missing nights (%d left, %d right) in %d records, ",
         "%d assessed, %d material (>%.0f%% of the season censored).\n"),
  summary_list$n_cens_nights, summary_list$n_left, summary_list$n_right,
  summary_list$n_records, summary_list$n_assessed, summary_list$n_material,
  100 * material))
cat(sprintf("  censored share: median %.2f%%, max %.1f%% | largest shift: p50 %.1f d, p95 %.1f d\n",
            100 * summary_list$median_share, 100 * summary_list$max_share,
            summary_list$max_shift_p50, summary_list$max_shift_p95))
print(headline)
cat("  written to data/processed/censoring_check.rds\n")
