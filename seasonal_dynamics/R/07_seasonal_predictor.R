# =============================================================================
# seasonal_dynamics/R/07_seasonal_predictor.R
# -----------------------------------------------------------------------------
# Assemble the seasonal results into a predictor of the nightly distribution
# of a season's catch, and score it out of sample.
#
#   centre  the median catch date, from the split form of 03: the year mean of
#           a development-model date over the network plus the locality's
#           departure from it, each with its own slope. Three drivers are
#           carried: BioSIM flight (the adopted centre, 08), BioSIM pupation
#           and bayessbw pupation (the drivers the predictor was developed
#           with, kept so that the log's development section reproduces)
#   width   the central 90% of the season, from the season total (04), on the
#           widths with their sampling bias removed, so that the curve
#           describes the season and the count noise is left to the nights.
#           Variant: the predicted spread of the driver's own distribution
#           (p95 - p05) as a second covariate.
#   shape   symmetric about the median (decision in the log): normal, logistic
#           or Laplace, with its scale set so that its central 90% equals the
#           predicted width
#   centre uncertainty   the predicted centre is itself uncertain by several
#           days. The predictive curve is the season curve averaged over that
#           uncertainty: a normal error with the residual SD of the centre
#           model in the training data, integrated on a fixed grid. Variant
#           `integrate = TRUE`; without it the curve is placed as if its centre
#           were known.
#
# NIGHTLY SHARES
#   With F the cumulative distribution of the curve, the share of night d is
#   F(d) - F(d - 1): the observed percentiles treat the catch of night d as
#   accumulated between d - 1 and d (R/functions/phenology.R), so this places
#   the curve's median on the same scale as the observed median. Shares are
#   mixed with a uniform floor of `eps` (a fraction of the season spread evenly
#   over its nights) so that no night has a share of zero, then normalised over
#   the record's counted nights.
#
# VALIDATION
#   leave-one-year-out and leave-one-locality-out; every component (centre,
#   width) is refitted on the training data of each fold.
#   Scores: log score per moth (sum of moths x log share over all nights,
#   divided by the moths -- the multinomial log-likelihood of the catch given
#   its total), mean absolute error of the median,
#   and mean error of the central 90%.
#   Oracle variants, leave-one-year-out only, replace the predicted centre
#   and/or width by the observed ones: they show what each component costs.
#
# INPUT   data/processed/seasonal/median_vs_pupation.rds (analysis set),
#         width_vs_catch.rds, pupation_dates.rds
#         data/processed/biosim_phenology.rds (R/03: the flight series)
#         data/processed/trap_daily.rds
#
# OUTPUT  data/processed/seasonal/seasonal_predictor.rds
#         list(scores, oracle, nights, fits, summary)
#
# REQUIREMENTS  dplyr
#
# USAGE  source("seasonal_dynamics/R/07_seasonal_predictor.R")   # from the project root
# =============================================================================

library(dplyr)

source(file.path("seasonal_dynamics", "R", "functions", "seasonal_curve.R"))

out_dir <- file.path("data", "processed", "seasonal")
eps     <- 0.001


## ---- 1. Inputs -------------------------------------------------------------------

## The analysis set of 03 (usable for phenology, complete BioSIM cohort).
base <- readRDS(file.path(out_dir, "median_vs_pupation.rds"))$data %>%
  distinct(flightID, Location, Year, p50)
widths <- readRDS(file.path(out_dir, "width_vs_catch.rds"))$data %>%
  select(flightID, width90, width90_adj, log_total)

## Each driver: its median date (x) and the spread of its distribution.
pupation <- readRDS(file.path(out_dir, "pupation_dates.rds"))
flight   <- readRDS(file.path("data", "processed", "biosim_phenology.rds")) %>%
  filter(series == "Flight_mf")
drivers <- bind_rows(
  transmute(flight, Location, Year, driver = "BioSIM flight",
            x = p50, spread = p95 - p05),
  transmute(pupation, Location, Year, driver = "BioSIM pupation",
            x = biosim_pup_p50, spread = biosim_pup_p95 - biosim_pup_p05),
  transmute(pupation, Location, Year, driver = "bayessbw pupation",
            x = bayes_pup_p50, spread = bayes_pup_p95 - bayes_pup_p05)
)
driver_names <- c("BioSIM flight", "BioSIM pupation", "bayessbw pupation")

## Split terms over the analysis set, per driver and year.
dat <- base %>%
  inner_join(widths, by = "flightID") %>%
  inner_join(drivers, by = c("Location", "Year"), relationship = "one-to-many") %>%
  group_by(driver, Year) %>%
  mutate(xY = mean(x), xW = x - mean(x)) %>%
  ungroup()

nights <- readRDS(file.path("data", "processed", "trap_daily.rds")) %>%
  filter(status == "counted", flightID %in% dat$flightID) %>%
  select(flightID, doy, moths) %>%
  arrange(flightID, doy)

stopifnot(!any(is.na(dat$width90_adj)), !any(is.na(dat$x)), !any(is.na(dat$spread)),
          all(table(dat$driver) == nrow(base)))


## ---- 2. Curves --------------------------------------------------------------------
## Shapes and nightly shares: seasonal_dynamics/R/functions/seasonal_curve.R.

shapes <- season_shapes


## ---- 3. Component models ------------------------------------------------------------

centre_fit <- function(train) stats::lm(p50 ~ xY + xW, data = train)
width_fits <- list(
  catch          = function(train) stats::lm(width90_adj ~ log_total, data = train),
  `catch + driver spread` =
    function(train) stats::lm(width90_adj ~ log_total + spread, data = train)
)

## Predicted (centre, width) for the test records of one fold. `oracle`
## substitutes the observed value for one or both components.
predict_fold <- function(train, test, width_model, oracle = "none") {
  cf <- centre_fit(train)
  mu <- stats::predict(cf, newdata = test)
  w  <- stats::predict(width_fits[[width_model]](train), newdata = test)
  sigma <- stats::sigma(cf)            # residual SD of the centre, training data
  if (oracle %in% c("centre", "both")) { mu <- test$p50; sigma <- 0 }
  if (oracle %in% c("width", "both"))  w  <- test$width90_adj
  tibble::tibble(flightID = test$flightID, mu = mu, w = pmax(w, 1), sigma = sigma)
}


## ---- 4. Scoring -----------------------------------------------------------------------

score_nights <- function(pred, shape, integrate) {
  pred %>%
    inner_join(nights, by = "flightID", relationship = "one-to-many") %>%
    group_by(flightID) %>%
    mutate(share = night_shares(doy, first(mu), first(w), shape,
                                sigma = if (integrate) first(sigma) else 0,
                                eps = eps)) %>%
    ungroup()
}

run <- function(driver, scheme, width_model, shape, oracle = "none",
                integrate = FALSE, keep_nights = FALSE) {
  d <- dat[dat$driver == driver, ]
  fold_var <- if (scheme == "leave-one-year-out") "Year" else "Location"
  pred <- bind_rows(lapply(unique(d[[fold_var]]), function(f) {
    predict_fold(d[d[[fold_var]] != f, ], d[d[[fold_var]] == f, ],
                 width_model, oracle)
  }))
  sn <- score_nights(pred, shape, integrate)
  rec <- pred %>% left_join(select(d, flightID, p50, width90, width90_adj),
                            by = "flightID")
  out <- tibble::tibble(
    driver = driver, scheme = scheme, width_model = width_model,
    shape = shape, oracle = oracle, integrate = integrate,
    sigma = mean(pred$sigma), n = nrow(rec), n_nights = nrow(sn),
    log_score = sum(sn$moths * log(sn$share)) / sum(sn$moths),
    mae_median = mean(abs(rec$mu - rec$p50)),
    width_err_adj = mean(rec$w - rec$width90_adj),
    width_err_obs = mean(rec$w - rec$width90),
    mae_width_obs = mean(abs(rec$w - rec$width90))
  )
  if (keep_nights) attr(out, "nights") <- sn %>% mutate(driver = driver, scheme = scheme)
  out
}

grid <- expand.grid(driver = driver_names,
                    scheme = c("leave-one-year-out", "leave-one-locality-out"),
                    width_model = names(width_fits),
                    shape = names(shapes), integrate = c(FALSE, TRUE),
                    stringsAsFactors = FALSE)

scores <- bind_rows(lapply(seq_len(nrow(grid)), function(i)
  run(grid$driver[i], grid$scheme[i], grid$width_model[i], grid$shape[i],
      integrate = grid$integrate[i]))) %>%
  arrange(scheme, desc(log_score))

## Oracles: leave-one-year-out, logistic, width from catch.
oracle <- bind_rows(lapply(driver_names, function(dr)
  bind_rows(lapply(c("none", "width", "centre", "both"), function(o)
    run(dr, "leave-one-year-out", "catch", "logistic", oracle = o)))))


## ---- 5. Full-data fits (the predictor to use) ------------------------------------------

## What seasonal_dynamics/R/functions/predict_season.R needs: the centre
## coefficients and its residual SD (the centre uncertainty), and the width
## coefficients on the season total.
fits <- lapply(setNames(driver_names, driver_names), function(dr) {
  d  <- dat[dat$driver == dr, ]
  cf <- centre_fit(d)
  list(centre = stats::coef(cf), sigma = stats::sigma(cf),
       width  = stats::coef(width_fits[["catch"]](d)),
       shape = "logistic", eps = eps, n = nrow(d))
})


## ---- 6. Write and report ------------------------------------------------------------------

best <- scores %>% filter(scheme == "leave-one-year-out") %>% slice(1)
summary_list <- list(n = n_distinct(dat$flightID), n_nights = nrow(nights),
                     adopted = "BioSIM flight",
                     eps = eps, best = best, fits = fits)

nights_best <- attr(run(best$driver, best$scheme, best$width_model, best$shape,
                        integrate = best$integrate, keep_nights = TRUE), "nights")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(scores = scores, oracle = oracle,
             nights = nights_best, fits = fits, summary = summary_list),
        file.path(out_dir, "seasonal_predictor.rds"))

cat(sprintf("Seasonal predictor: %d locality-years, %d nights.\n",
            summary_list$n, summary_list$n_nights))
for (sc in unique(scores$scheme)) {
  cat("  ", sc, "\n", sep = "")
  s <- filter(scores, scheme == sc)
  for (i in seq_len(nrow(s))) cat(sprintf(
    "    %-17s %-22s %-9s %-10s log score %7.3f | median MAE %4.2f d | width err %+5.2f d (adj) %+5.2f d (obs)\n",
    s$driver[i], s$width_model[i], s$shape[i],
    if (s$integrate[i]) sprintf("sigma %.1f", s$sigma[i]) else "fixed", s$log_score[i], s$mae_median[i],
    s$width_err_adj[i], s$width_err_obs[i]))
}
cat("  oracles (leave-one-year-out, logistic, width from catch):\n")
for (i in seq_len(nrow(oracle))) cat(sprintf(
  "    %-17s observed %-7s log score %7.3f | median MAE %4.2f d\n",
  oracle$driver[i], oracle$oracle[i], oracle$log_score[i], oracle$mae_median[i]))
cat("  written to data/processed/seasonal/seasonal_predictor.rds\n")
