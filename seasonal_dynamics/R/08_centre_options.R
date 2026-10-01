# =============================================================================
# seasonal_dynamics/R/08_centre_options.R
# -----------------------------------------------------------------------------
# Alternatives for the centre of the seasonal predictor (07), whose error is
# almost all of the predictor's error. Every option uses the split form of 03
# -- the year mean of the driver over the network plus the locality's
# departure from it, each with its own slope -- and is scored on the same
# locality-years, both on the median catch date alone and inside the full
# predictor (width from the season total, logistic curve averaged over the
# centre's uncertainty).
#
#   OPTIONS
#     BioSIM pupation        the current centre (07)
#     bayessbw pupation      the current alternative
#     BioSIM emergence       median adult emergence date (R/03 Emergence_mf)
#     BioSIM flight          median of the modelled flight series (R/03
#                            Flight_mf)
#     bayessbw pupation + T  bayessbw pupation plus the mean temperature over
#                            the `t_window` days after the locality's
#                            predicted median pupation, split like the driver
#                            into a network-mean and a local-departure term.
#                            The temperature of the pupal and pre-flight
#                            period is the obvious thing the pupation date
#                            alone cannot know.
#
# TEMPERATURE   daily mean of the hourly CaSR series (R/02), averaged over
#               days pupation + 0 ... pupation + t_window - 1, at the trap.
#
# VALIDATION    leave-one-year-out and leave-one-locality-out, everything
#               refitted per fold; coefficients from mixed models with a
#               random year intercept.
#
# INPUT   data/processed/seasonal/median_vs_pupation.rds   (analysis set, 03)
#         data/processed/seasonal/width_vs_catch.rds       (04)
#         data/processed/seasonal/pupation_dates.rds       (02)
#         data/processed/biosim_phenology.rds               (R/03)
#         data/processed/weather_hourly.rds                 (R/02)
#         data/processed/trap_daily.rds                     (R/01)
#
# OUTPUT  data/processed/seasonal/centre_options.rds
#         list(data, scores, coefs, summary)
#
# REQUIREMENTS  dplyr, lme4
#
# USAGE  source("seasonal_dynamics/R/08_centre_options.R")   # from the project root
# =============================================================================

library(dplyr)
library(lme4)

source(file.path("seasonal_dynamics", "R", "functions", "seasonal_curve.R"))

out_dir   <- file.path("data", "processed", "seasonal")
t_windows <- c(14, 21)
eps       <- 0.001


## ---- 1. The analysis set and the candidate drivers ---------------------------------

base <- readRDS(file.path(out_dir, "median_vs_pupation.rds"))$data %>%
  distinct(flightID, Location, Year, p50, season_total) %>%
  inner_join(readRDS(file.path(out_dir, "width_vs_catch.rds"))$data %>%
               select(flightID, width90_adj, log_total), by = "flightID")

pup <- readRDS(file.path(out_dir, "pupation_dates.rds")) %>%
  select(Location, Year, biosim_pup = biosim_pup_p50, bayes_pup = bayes_pup_p50)

pheno <- readRDS(file.path("data", "processed", "biosim_phenology.rds")) %>%
  filter(series %in% c("Emergence_mf", "Flight_mf")) %>%
  select(Location, Year, series, p50) %>%
  tidyr::pivot_wider(names_from = series, values_from = p50) %>%
  rename(biosim_emerg = Emergence_mf, biosim_flight = Flight_mf)

base <- base %>%
  left_join(pup,   by = c("Location", "Year")) %>%
  left_join(pheno, by = c("Location", "Year"))


## ---- 2. Temperature after predicted pupation ----------------------------------------

weather <- readRDS(file.path("data", "processed", "weather_hourly.rds"))

daily_temp <- weather %>%
  semi_join(base, by = c("Location", "Year")) %>%
  select(Location, Year, Weather) %>%
  tidyr::unnest(Weather) %>%
  mutate(doy = as.integer(format(datetime_est, "%j"))) %>%
  group_by(Location, Year, doy) %>%
  summarise(tmean = mean(Temp, na.rm = TRUE), .groups = "drop")

post_temp <- function(loc, yr, start, len) {
  d <- daily_temp[daily_temp$Location == loc & daily_temp$Year == yr, ]
  days <- round(start) + 0:(len - 1)
  v <- d$tmean[match(days, d$doy)]
  if (any(is.na(v))) NA_real_ else mean(v)
}

for (L in t_windows) {
  base[[paste0("T", L)]] <- mapply(post_temp, base$Location, base$Year,
                                   base$bayes_pup, L)
}
stopifnot(!any(is.na(base[paste0("T", t_windows)])))


## ---- 3. Split terms -----------------------------------------------------------------
## For each driver x: xY = mean of x over the analysis set's localities that
## year; xW = x - xY. pupY/pupW for bayessbw are reused by the temperature
## specifications.

split_terms <- function(d, var) {
  d %>% group_by(Year) %>%
    mutate("{var}_Y" := mean(.data[[var]]), "{var}_W" := .data[[var]] - mean(.data[[var]])) %>%
    ungroup()
}
for (v in c("biosim_pup", "bayes_pup", "biosim_emerg", "biosim_flight",
            paste0("T", t_windows))) base <- split_terms(base, v)

specs <- list(
  `BioSIM pupation`         = c("biosim_pup_Y", "biosim_pup_W"),
  `bayessbw pupation`       = c("bayes_pup_Y", "bayes_pup_W"),
  `BioSIM emergence`        = c("biosim_emerg_Y", "biosim_emerg_W"),
  `BioSIM flight`           = c("biosim_flight_Y", "biosim_flight_W"),
  `bayessbw pupation + T14` = c("bayes_pup_Y", "bayes_pup_W", "T14_Y", "T14_W"),
  `bayessbw pupation + T21` = c("bayes_pup_Y", "bayes_pup_W", "T21_Y", "T21_W")
)


## ---- 4. Cross-validation ---------------------------------------------------------------

nights <- readRDS(file.path("data", "processed", "trap_daily.rds")) %>%
  filter(status == "counted", flightID %in% base$flightID) %>%
  select(flightID, doy, moths) %>%
  arrange(flightID, doy)

run_spec <- function(spec, scheme) {
  terms <- specs[[spec]]
  f_c <- stats::as.formula(paste("p50 ~", paste(terms, collapse = " + ")))
  fold_var <- if (scheme == "leave-one-year-out") "Year" else "Location"
  pred <- bind_rows(lapply(unique(base[[fold_var]]), function(g) {
    train <- base[base[[fold_var]] != g, ]
    test  <- base[base[[fold_var]] == g, ]
    cf <- stats::lm(f_c, data = train)
    wf <- stats::lm(width90_adj ~ log_total, data = train)
    tibble::tibble(flightID = test$flightID, obs = test$p50,
                   season_total = test$season_total,
                   mu = stats::predict(cf, newdata = test),
                   w = pmax(stats::predict(wf, newdata = test), 1),
                   sigma = stats::sigma(cf))
  }))
  sn <- pred %>%
    inner_join(nights, by = "flightID", relationship = "one-to-many") %>%
    group_by(flightID) %>%
    mutate(share = night_shares(doy, first(mu), first(w), "logistic",
                                sigma = first(sigma), eps = eps)) %>%
    ungroup()
  tibble::tibble(
    spec = spec, scheme = scheme, n = nrow(pred),
    mae = mean(abs(pred$mu - pred$obs)),
    rmse = sqrt(mean((pred$mu - pred$obs)^2)),
    mae_moth = stats::weighted.mean(abs(pred$mu - pred$obs), pred$season_total),
    r = stats::cor(pred$mu, pred$obs),
    sigma = mean(pred$sigma),
    log_score = sum(sn$moths * log(sn$share)) / sum(sn$moths)
  )
}

scores <- bind_rows(lapply(names(specs), function(sp)
  bind_rows(lapply(c("leave-one-year-out", "leave-one-locality-out"),
                   function(sc) run_spec(sp, sc))))) %>%
  arrange(scheme, desc(log_score))


## ---- 5. Coefficients (mixed model, random year intercept) ---------------------------------

coefs <- bind_rows(lapply(names(specs), function(sp) {
  f <- stats::as.formula(paste("p50 ~", paste(specs[[sp]], collapse = " + "),
                               "+ (1 | Year)"))
  m  <- lmer(f, data = base, REML = TRUE)
  fe <- lme4::fixef(m)[-1]
  se <- sqrt(diag(as.matrix(stats::vcov(m))))[-1]
  tibble::tibble(spec = sp, term = names(fe), estimate = unname(fe),
                 lo = unname(fe - 1.96 * se), hi = unname(fe + 1.96 * se))
}))


## ---- 6. Write and report --------------------------------------------------------------

summary_list <- list(n = nrow(base), n_years = n_distinct(base$Year),
                     t_windows = t_windows,
                     temp_range = sapply(paste0("T", t_windows),
                                         function(v) range(base[[v]])))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(data = base, scores = scores, coefs = coefs, summary = summary_list),
        file.path(out_dir, "centre_options.rds"))

cat(sprintf("Centre options: %d locality-years, %d years.\n", nrow(base),
            n_distinct(base$Year)))
for (sc in unique(scores$scheme)) {
  cat("  ", sc, "\n", sep = "")
  s <- filter(scores, scheme == sc)
  for (i in seq_len(nrow(s))) cat(sprintf(
    "    %-25s log score %7.3f | MAE %4.2f d | RMSE %4.2f | MAE per moth %4.2f | r %4.2f | sigma %3.1f\n",
    s$spec[i], s$log_score[i], s$mae[i], s$rmse[i], s$mae_moth[i], s$r[i], s$sigma[i]))
}
cat("  coefficients:\n")
for (i in seq_len(nrow(coefs))) cat(sprintf(
  "    %-25s %-16s %6.3f [%6.3f, %6.3f]\n", coefs$spec[i], coefs$term[i],
  coefs$estimate[i], coefs$lo[i], coefs$hi[i]))
cat("  written to data/processed/seasonal/centre_options.rds\n")
