# =============================================================================
# seasonal_dynamics/R/03_median_vs_pupation.R
# -----------------------------------------------------------------------------
# Model the observed median catch date of each locality-year as a function of
# the median pupation date predicted by BioSIM or by bayessbw.
#
# INPUT   data/processed/seasonal/season_percentiles.rds  (01)
#         data/processed/seasonal/pupation_dates.rds      (02)
#         data/processed/biosim_phenology.rds              (R/03: cohort flag)
#
# ANALYSIS SET
#   Locality-years usable for phenology (R/01), with a complete BioSIM cohort
#   (R/03 `flight_complete`) and a BioSIM cohort that fully pupated (02). The
#   same set is used for both drivers so that they are compared on identical
#   data.
#
# MODELS, for each driver (pup = predicted median pupation date)
#   lag      p50 = pup + c                         a fixed pupal-plus-pre-flight
#                                                  interval; slope fixed at 1
#   linear   p50 = a + b pup                       slope estimated
#   split    p50 = a + bY pupY + bW (pup - pupY)   pupY = mean of pup over the
#                                                  localities of that year.
#            bY is the response to a year being early or late, bW the response
#            to a locality being early or late relative to the others that
#            year. The daily analysis found the models' spatial gradient far
#            steeper than the observed one; this split measures it directly.
#            pupY uses predicted pupation only, so it is available at
#            prediction time.
#   Baselines, for scale:
#     climatology   p50 = the mean observed p50 of the training data
#     same-year     p50 = the mean observed p50 of the OTHER localities that
#                   year -- not a forecast (it uses observations of the year
#                   being predicted) but the floor any year-level model faces
#
# VALIDATION
#   leave-one-year-out (LOYO): every locality-year of a year is predicted
#   from a fit to the other years -- the test for a new year.
#   leave-one-locality-out (LOLO): every year of a locality is predicted from
#   a fit to the other localities -- the test for a new place.
#
# INFERENCE
#   Coefficients of the full-data fits are reported from mixed models with a
#   random year intercept (lme4), because locality-years of one year share
#   that year's weather and are not independent; ordinary least squares would
#   overstate the precision of the year-level slope.
#
# OUTPUT  data/processed/seasonal/median_vs_pupation.rds
#         list(data, cv, scores, coefs, summary)
#
# REQUIREMENTS  dplyr, lme4
#
# USAGE  source("seasonal_dynamics/R/03_median_vs_pupation.R")   # from the project root
# =============================================================================

library(dplyr)
library(lme4)

out_dir <- file.path("data", "processed", "seasonal")


## ---- 1. Analysis set -------------------------------------------------------------

season   <- readRDS(file.path(out_dir, "season_percentiles.rds"))
pupation <- readRDS(file.path(out_dir, "pupation_dates.rds"))
cohort   <- readRDS(file.path("data", "processed", "biosim_phenology.rds")) %>%
  distinct(Location, Year, flight_complete)

dat <- season %>%
  select(flightID, Location = locality, Year = year, latitude, longitude,
         p50, season_total, usable_phenology) %>%
  inner_join(pupation, by = c("Location", "Year")) %>%
  inner_join(cohort,   by = c("Location", "Year")) %>%
  filter(usable_phenology, flight_complete, biosim_pup_complete) %>%
  ## Long form: one row per locality-year and driver.
  tidyr::pivot_longer(c(biosim_pup_p50, bayes_pup_p50),
                      names_to = "driver", values_to = "pup") %>%
  mutate(driver = recode(driver, biosim_pup_p50 = "BioSIM",
                         bayes_pup_p50 = "bayessbw")) %>%
  group_by(driver, Year) %>%
  mutate(pupY = mean(pup), pupW = pup - pupY) %>%
  ungroup() %>%
  select(flightID, Location, Year, latitude, longitude, season_total,
         driver, p50, pup, pupY, pupW)

stopifnot(!any(is.na(dat$pup)), !any(is.na(dat$p50)))


## ---- 2. Model definitions ---------------------------------------------------------
## Each returns predictions for `test` from a fit to `train`.

fit_predict <- list(
  lag = function(train, test) {
    test$pup + mean(train$p50 - train$pup)
  },
  linear = function(train, test) {
    stats::predict(stats::lm(p50 ~ pup, data = train), newdata = test)
  },
  split = function(train, test) {
    stats::predict(stats::lm(p50 ~ pupY + pupW, data = train), newdata = test)
  },
  climatology = function(train, test) {
    rep(mean(train$p50), nrow(test))
  }
)

## The same-year floor is not a fitted model: computed once, per driver set.
same_year <- function(d) {
  d %>%
    group_by(Year) %>%
    mutate(pred = if (n() > 1) (sum(p50) - p50) / (n() - 1) else NA_real_) %>%
    ungroup() %>%
    pull(pred)
}


## ---- 3. Cross-validation ------------------------------------------------------------

cross_validate <- function(d, fold_var) {
  folds <- unique(d[[fold_var]])
  bind_rows(lapply(folds, function(f) {
    test  <- d[d[[fold_var]] == f, ]
    train <- d[d[[fold_var]] != f, ]
    ## pupY of the held-out data is recomputed from the held-out year itself,
    ## which is legitimate: it is a mean of PREDICTED pupation dates.
    bind_rows(lapply(names(fit_predict), function(m) {
      tibble::tibble(flightID = test$flightID, Year = test$Year,
                     Location = test$Location, model = m,
                     obs = test$p50, pred = fit_predict[[m]](train, test))
    }))
  }))
}

cv <- bind_rows(lapply(c("BioSIM", "bayessbw"), function(dr) {
  d <- filter(dat, driver == dr)
  bind_rows(
    mutate(cross_validate(d, "Year"),     scheme = "leave-one-year-out"),
    mutate(cross_validate(d, "Location"), scheme = "leave-one-locality-out"),
    tibble::tibble(flightID = d$flightID, Year = d$Year, Location = d$Location,
                   model = "same-year", obs = d$p50, pred = same_year(d),
                   scheme = "leave-one-year-out"),
    tibble::tibble(flightID = d$flightID, Year = d$Year, Location = d$Location,
                   model = "same-year", obs = d$p50, pred = same_year(d),
                   scheme = "leave-one-locality-out")
  ) %>% mutate(driver = dr)
}))

scores <- cv %>%
  filter(!is.na(pred)) %>%
  mutate(err = pred - obs) %>%
  group_by(scheme, driver, model) %>%
  summarise(n = n(), rmse = sqrt(mean(err^2)), mae = mean(abs(err)),
            bias = mean(err), r = stats::cor(pred, obs), .groups = "drop") %>%
  arrange(scheme, driver, rmse)


## ---- 4. Coefficients from the full data -------------------------------------------

coef_table <- function(fit, term) {
  est <- lme4::fixef(fit)[term]
  se  <- sqrt(diag(as.matrix(stats::vcov(fit))))[term]
  tibble::tibble(term = term, estimate = unname(est),
                 lo = unname(est - 1.96 * se), hi = unname(est + 1.96 * se))
}

coefs <- bind_rows(lapply(c("BioSIM", "bayessbw"), function(dr) {
  d <- filter(dat, driver == dr)
  lin   <- lmer(p50 ~ pup + (1 | Year), data = d, REML = TRUE)
  split <- lmer(p50 ~ pupY + pupW + (1 | Year), data = d, REML = TRUE)
  bind_rows(
    mutate(coef_table(lin, "pup"), model = "linear"),
    mutate(coef_table(split, c("pupY", "pupW")), model = "split"),
    tibble::tibble(term = "lag", model = "lag",
                   estimate = mean(d$p50 - d$pup),
                   lo = NA_real_, hi = NA_real_)
  ) %>% mutate(driver = dr)
}))


## ---- 5. Write and report ---------------------------------------------------------

summary_list <- list(
  n_rows      = n_distinct(dat$flightID),
  n_years     = n_distinct(dat$Year),
  n_locations = n_distinct(dat$Location),
  scores      = scores,
  coefs       = coefs
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(data = dat, cv = cv, scores = scores, coefs = coefs,
             summary = summary_list),
        file.path(out_dir, "median_vs_pupation.rds"))

cat(sprintf("Median catch date vs predicted pupation: %d locality-years, %d years, %d localities.\n",
            summary_list$n_rows, summary_list$n_years, summary_list$n_locations))
for (sc in unique(scores$scheme)) {
  cat("  ", sc, "\n", sep = "")
  s <- filter(scores, scheme == sc)
  for (i in seq_len(nrow(s))) cat(sprintf(
    "    %-9s %-12s RMSE %5.2f d | MAE %5.2f d | bias %+5.2f | r %5.2f\n",
    s$driver[i], s$model[i], s$rmse[i], s$mae[i], s$bias[i], s$r[i]))
}
cat("  coefficients (mixed model, random year intercept):\n")
for (i in seq_len(nrow(coefs))) cat(sprintf(
  "    %-9s %-7s %-5s %6.2f  [%6.2f, %6.2f]\n", coefs$driver[i], coefs$model[i],
  coefs$term[i], coefs$estimate[i], coefs$lo[i], coefs$hi[i]))
cat("  written to data/processed/seasonal/median_vs_pupation.rds\n")
