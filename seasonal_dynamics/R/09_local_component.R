# =============================================================================
# seasonal_dynamics/R/09_local_component.R
# -----------------------------------------------------------------------------
# The local component of the median catch date: how a locality departs from
# the other localities of its year. The adopted centre (BioSIM flight, 07/08)
# follows the year signal with a slope near 0.9 but the local departure of
# BioSIM's flight date with a slope near 0.1, so almost all of the
# within-year variation is left in its residuals. This script asks what that
# residual is made of and whether anything known at prediction time explains
# it.
#
#   1. DECOMPOSITION  residuals of the adopted centre model split into a year
#      part (what the driver misses about the year), a persistent locality
#      part (habitual lateness), and a locality-year remainder, from which the
#      sampling variance of the observed median (multinomial resampling of
#      each record's own catch) is subtracted to leave the transient part.
#   2. SPATIAL STRUCTURE  correlation of within-year residuals between pairs of
#      localities trapped in the same year, by distance; and between years at
#      the same locality. A smooth spatial field would show as positive
#      correlation at short distance; pure trap-level noise as none.
#   2b. YEAR-SPECIFIC GRADIENTS  in years with at least `min_grad_n`
#      localities, a plane (latitude + longitude) fitted to the residuals of
#      each year: the share of within-year residual variance it explains,
#      against a null from permuting residuals among the localities of a year.
#      The same planes fitted to the observed median and to BioSIM's flight
#      date give each year's observed and modelled gradients, which are
#      compared (mean, year-to-year SD net of the slopes' own sampling error,
#      correlation).
#   3. CANDIDATE LOCAL TERMS  each added to the centre model p50 ~ xY + xW
#      and scored like 08 (leave-one-year-out, leave-one-locality-out, log
#      score per moth of the full predictor). All covariates are entered as
#      departures from their year mean over the analysis set, so that they
#      can only act on the local component:
#        geography    latitude, longitude, elevation
#        temperature  local anomaly of the mean CaSR temperature over two
#                     windows fixed in advance: 1 May - 30 June (larval
#                     development) and 15 June - 15 July (pupal and pre-flight
#                     period)
#        catch        local anomaly of log10 season total (an input of the
#                     predictor already)
#        habitual lateness  a random locality intercept: the locality's
#                     shrunken mean residual in its training years. Known only
#                     for localities trapped before, so it can help a new
#                     year but not a new locality (where it falls back to 0
#                     and the centre uncertainty widens by the locality SD).
#
# INPUT   data/processed/seasonal/median_vs_pupation.rds   (analysis set, 03)
#         data/processed/seasonal/width_vs_catch.rds       (04)
#         data/processed/biosim_phenology.rds               (R/03)
#         data/processed/weather_hourly.rds                 (R/02: coordinates, temperature)
#         data/processed/locality_elevation.rds             (R/03)
#         data/processed/trap_daily.rds                     (R/01)
#
# OUTPUT  data/processed/seasonal/local_component.rds
#         list(data, variance, correlogram, gradients, scores, coefs,
#              lateness, summary)
#
# REQUIREMENTS  dplyr, tidyr, lme4
#
# USAGE  source("seasonal_dynamics/R/09_local_component.R")   # from the project root
# =============================================================================

library(dplyr)
library(lme4)

source(file.path("R", "functions", "phenology.R"))
source(file.path("seasonal_dynamics", "R", "functions", "seasonal_curve.R"))

out_dir <- file.path("data", "processed", "seasonal")
eps     <- 0.001
n_boot  <- 200            # resamples per record for the sampling variance
n_yboot <- 500            # resamples of years for correlogram intervals
t_windows <- list(MayJun = c(121, 181), JunJul = c(166, 196))   # DOY, non-leap
dist_breaks <- c(0, 25, 50, 100, 200, Inf)                       # km
min_grad_n  <- 6          # localities a year needs for its own gradient
n_perm      <- 500        # permutations for the gradient null
set.seed(20261001)


## ---- 1. The analysis set ------------------------------------------------------------

base <- readRDS(file.path(out_dir, "median_vs_pupation.rds"))$data %>%
  distinct(flightID, Location, Year, p50, season_total) %>%
  inner_join(readRDS(file.path(out_dir, "width_vs_catch.rds"))$data %>%
               select(flightID, width90_adj, log_total), by = "flightID")

flight <- readRDS(file.path("data", "processed", "biosim_phenology.rds")) %>%
  filter(series == "Flight_mf") %>%
  select(Location, Year, x = p50)

weather <- readRDS(file.path("data", "processed", "weather_hourly.rds"))
coords  <- weather %>% distinct(Location, lat = Latitude, lon = Longitude) %>%
  group_by(Location) %>% slice(1) %>% ungroup()
elev    <- readRDS(file.path("data", "processed", "locality_elevation.rds"))

base <- base %>%
  left_join(flight, by = c("Location", "Year")) %>%
  left_join(coords, by = "Location") %>%
  left_join(elev,   by = "Location")


## ---- 2. Local temperature ------------------------------------------------------------

daily_temp <- weather %>%
  semi_join(base, by = c("Location", "Year")) %>%
  select(Location, Year, Weather) %>%
  tidyr::unnest(Weather) %>%
  mutate(doy = as.integer(format(datetime_est, "%j"))) %>%
  group_by(Location, Year, doy) %>%
  summarise(tmean = mean(Temp, na.rm = TRUE), .groups = "drop")

for (nm in names(t_windows)) {
  win <- t_windows[[nm]]
  tw <- daily_temp %>%
    filter(doy >= win[1], doy <= win[2]) %>%
    group_by(Location, Year) %>%
    summarise("T_{nm}" := mean(tmean), n_days = n(), .groups = "drop")
  stopifnot(all(tw$n_days == diff(win) + 1))
  base <- left_join(base, select(tw, -n_days), by = c("Location", "Year"))
}

stopifnot(!anyNA(base[c("x", "lat", "lon", "elev_m", "T_MayJun", "T_JunJul")]))


## ---- 3. Split terms ------------------------------------------------------------------
## xY/xW for the driver as in 07; every local covariate as a departure from its
## year mean (suffix _W). Elevation in hundreds of metres.

base <- base %>%
  mutate(elev100 = elev_m / 100) %>%
  group_by(Year) %>%
  mutate(xY = mean(x), xW = x - mean(x),
         across(c(lat, lon, elev100, T_MayJun, T_JunJul, log_total),
                ~ .x - mean(.x), .names = "{.col}_W"),
         n_year = n()) %>%
  ungroup()


## ---- 4. Decomposition of the centre residual ---------------------------------------------

nights <- readRDS(file.path("data", "processed", "trap_daily.rds")) %>%
  filter(status == "counted", flightID %in% base$flightID) %>%
  select(flightID, doy, moths) %>%
  arrange(flightID, doy)

## Sampling variance of each observed median: the record's own catch
## resampled with its own total.
boot_p50 <- function(id, total) {
  n_id <- nights[nights$flightID == id, ]
  draws <- stats::rmultinom(n_boot, round(total), n_id$moths / sum(n_id$moths))
  p <- apply(draws, 2, function(v) curve_percentiles(n_id$doy, v, 0.5))
  tibble::tibble(flightID = id, p50_se = stats::sd(p))
}
base <- left_join(base, bind_rows(Map(boot_p50, base$flightID, base$season_total)),
                  by = "flightID")

centre_lm <- stats::lm(p50 ~ xY + xW, data = base)
base$resid <- stats::residuals(centre_lm)

vc_of <- function(m) {
  v <- as.data.frame(lme4::VarCorr(m))
  setNames(v$vcov, v$grp)
}
m_res <- lmer(resid ~ 1 + (1 | Year) + (1 | Location), data = base, REML = TRUE)
m_obs <- lmer(p50   ~ 1 + (1 | Year) + (1 | Location), data = base, REML = TRUE)
samp_var <- mean(base$p50_se^2)

## Parametric bootstrap of the residual decomposition for intervals.
vc_boot <- lme4::bootMer(m_res, function(m) vc_of(m), nsim = 200, seed = 1)$t

variance <- bind_rows(
  tibble::tibble(model = "observed median", component = names(vc_of(m_obs)),
                 var = vc_of(m_obs), lo = NA_real_, hi = NA_real_),
  tibble::tibble(model = "centre residual", component = names(vc_of(m_res)),
                 var = vc_of(m_res),
                 lo = apply(vc_boot, 2, stats::quantile, 0.025),
                 hi = apply(vc_boot, 2, stats::quantile, 0.975))
) %>%
  mutate(component = recode(component, Year = "year", Location = "locality",
                            Residual = "locality-year"))
## Split the locality-year remainder into sampling and transient parts.
variance <- bind_rows(
  variance,
  tibble::tibble(model = c("observed median", "centre residual"),
                 component = "sampling", var = samp_var, lo = NA, hi = NA)
) %>%
  group_by(model) %>%
  mutate(var = ifelse(component == "locality-year", var - samp_var, var),
         lo = ifelse(component == "locality-year", lo - samp_var, lo),
         hi = ifelse(component == "locality-year", hi - samp_var, hi),
         component = ifelse(component == "locality-year", "transient", component),
         share = var / sum(var), sd = sqrt(pmax(var, 0))) %>%
  ungroup()


## ---- 5. Spatial structure of the within-year residual ------------------------------------
## Within-year residual: the centre residual minus its year mean (the year
## part is not local). Centring within a year makes the expected correlation
## between two records of the same year -1 / (n_year - 1) under no structure;
## that reference is reported with the observed values.

haversine <- function(lat1, lon1, lat2, lon2) {
  r <- pi / 180
  a <- sin((lat2 - lat1) * r / 2)^2 +
       cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  2 * 6371 * asin(sqrt(a))
}

base <- base %>% group_by(Year) %>% mutate(resid_W = resid - mean(resid)) %>% ungroup()
s2 <- sum(base$resid_W^2) / (nrow(base) - n_distinct(base$Year))

pairs <- base %>%
  select(Year, i = flightID, Li = Location, lat_i = lat, lon_i = lon, e_i = resid_W, n_year) %>%
  inner_join(base %>% select(Year, j = flightID, Lj = Location, lat_j = lat, lon_j = lon, e_j = resid_W),
             by = "Year", relationship = "many-to-many") %>%
  filter(i < j) %>%
  mutate(km = haversine(lat_i, lon_i, lat_j, lon_j),
         bin = cut(km, dist_breaks, right = FALSE, dig.lab = 4),
         prod = e_i * e_j, null = -1 / (n_year - 1))

correlogram_of <- function(p) {
  p %>% group_by(bin) %>%
    summarise(n_pairs = n(), km = mean(km), r = mean(prod) / s2,
              null = mean(null), .groups = "drop")
}
corr_obs <- correlogram_of(pairs)
yrs <- unique(base$Year)
corr_boot <- bind_rows(lapply(seq_len(n_yboot), function(b) {
  ys <- sample(yrs, replace = TRUE)
  bind_rows(lapply(ys, function(y) pairs[pairs$Year == y, ])) %>% correlogram_of()
}))
correlogram <- corr_obs %>%
  left_join(corr_boot %>% group_by(bin) %>%
              summarise(lo = stats::quantile(r, 0.025), hi = stats::quantile(r, 0.975)),
            by = "bin")

## Same locality, different years: correlation of within-year residuals.
same_loc <- base %>%
  select(Location, i = flightID, Yi = Year, e_i = resid_W) %>%
  inner_join(base %>% select(Location, j = flightID, Yj = Year, e_j = resid_W),
             by = "Location", relationship = "many-to-many") %>%
  filter(i < j)
same_loc_r <- mean(same_loc$e_i * same_loc$e_j) / s2


## ---- 5b. Year-specific gradients ------------------------------------------------------

grad_years <- base %>% count(Year) %>% filter(n >= min_grad_n) %>% pull(Year)
plane <- function(dd, v) {
  m <- stats::lm(stats::reformulate(c("lat", "lon"), v), data = dd)
  cs <- summary(m)$coefficients
  c(ss = sum(stats::residuals(m)^2), tot = sum((dd[[v]] - mean(dd[[v]]))^2),
    lat = cs["lat", 1], lat_se = cs["lat", 2], lon = cs["lon", 1], lon_se = cs["lon", 2])
}
gradients <- bind_rows(lapply(grad_years, function(y) {
  dd <- base[base$Year == y, ]
  r <- plane(dd, "resid"); o <- plane(dd, "p50"); x <- plane(dd, "x")
  tibble::tibble(Year = y, n = nrow(dd), ss = r[["ss"]], tot = r[["tot"]],
                 obs_lat = o[["lat"]], obs_lat_se = o[["lat_se"]],
                 obs_lon = o[["lon"]], obs_lon_se = o[["lon_se"]],
                 biosim_lat = x[["lat"]], biosim_lon = x[["lon"]])
}))
plane_share <- 1 - sum(gradients$ss) / sum(gradients$tot)
plane_null <- replicate(n_perm, {
  f <- sapply(grad_years, function(y) {
    dd <- base[base$Year == y, ]; dd$resid <- sample(dd$resid)
    plane(dd, "resid")[c("ss", "tot")]
  })
  1 - sum(f["ss", ]) / sum(f["tot", ])
})
## Year-to-year SD of the observed gradient net of its sampling error.
true_sd <- function(est, se) sqrt(max(stats::var(est) - mean(se^2), 0))
grad_summary <- list(
  n_years = length(grad_years), n = sum(gradients$n),
  share = plane_share, null_mean = mean(plane_null),
  null_q95 = unname(stats::quantile(plane_null, 0.95)),
  p = (1 + sum(plane_null >= plane_share)) / (1 + n_perm),
  obs_lat_mean = mean(gradients$obs_lat), obs_lat_sd = stats::sd(gradients$obs_lat),
  obs_lat_sd_true = true_sd(gradients$obs_lat, gradients$obs_lat_se),
  obs_lon_mean = mean(gradients$obs_lon), obs_lon_sd = stats::sd(gradients$obs_lon),
  obs_lon_sd_true = true_sd(gradients$obs_lon, gradients$obs_lon_se),
  biosim_lat_mean = mean(gradients$biosim_lat), biosim_lat_sd = stats::sd(gradients$biosim_lat),
  biosim_lon_mean = mean(gradients$biosim_lon), biosim_lon_sd = stats::sd(gradients$biosim_lon),
  r_lat = stats::cor(gradients$obs_lat, gradients$biosim_lat),
  r_lon = stats::cor(gradients$obs_lon, gradients$biosim_lon)
)


## ---- 6. Candidate local terms: cross-validation ---------------------------------------------

specs <- list(
  `adopted (xY + xW)`        = list(terms = character(0),                          re = FALSE),
  `+ geography`              = list(terms = c("lat_W", "lon_W", "elev100_W"),       re = FALSE),
  `+ elevation`              = list(terms = "elev100_W",                           re = FALSE),
  `+ temperature May-Jun`    = list(terms = "T_MayJun_W",                          re = FALSE),
  `+ temperature Jun-Jul`    = list(terms = "T_JunJul_W",                          re = FALSE),
  `+ catch`                  = list(terms = "log_total_W",                         re = FALSE),
  `+ habitual lateness`      = list(terms = character(0),                          re = "locality"),
  `+ habitual lateness, year-adjusted` = list(terms = character(0),               re = "locality + year")
)
for (nm in names(specs)) if (isFALSE(specs[[nm]]$re)) specs[[nm]]$re <- "none"

## One fold: fit on train, predict the centre (mu) and its uncertainty (sigma)
## for test. With a locality intercept the uncertainty of a known locality is
## the residual SD plus the posterior variance of its intercept; of a new
## locality, the residual SD plus the locality variance.
##   re = "locality"          p50 ~ xY + xW + (1 | Location): the intercept
##                            also absorbs the year errors of the few years a
##                            locality was trapped
##   re = "locality + year"   p50 ~ xY + xW + (1 | Year) + (1 | Location),
##                            predicted with the locality intercept only (a
##                            held-out year has no year intercept): habitual
##                            lateness net of the year errors; the year
##                            variance is added back to the uncertainty
predict_centre <- function(spec, train, test) {
  rhs <- paste(c("xY", "xW", spec$terms), collapse = " + ")
  if (spec$re == "none") {
    cf <- stats::lm(stats::as.formula(paste("p50 ~", rhs)), data = train)
    return(list(mu = unname(stats::predict(cf, newdata = test)),
                sigma = rep(stats::sigma(cf), nrow(test))))
  }
  yr <- spec$re == "locality + year"
  m <- lmer(stats::as.formula(paste("p50 ~", rhs, "+ (1 | Location)",
                                    if (yr) "+ (1 | Year)" else "")),
            data = train, REML = TRUE)
  mu <- unname(stats::predict(m, newdata = test, re.form = ~ (1 | Location),
                              allow.new.levels = TRUE))
  re <- lme4::ranef(m, condVar = TRUE)$Location
  pv <- setNames(as.numeric(attr(re, "postVar")), rownames(re))
  vc <- vc_of(m)
  extra <- ifelse(test$Location %in% names(pv), pv[test$Location], vc[["Location"]])
  if (yr) extra <- extra + vc[["Year"]]
  list(mu = mu, sigma = sqrt(stats::sigma(m)^2 + extra))
}

run_spec <- function(name, scheme) {
  spec <- specs[[name]]
  fold_var <- if (scheme == "leave-one-year-out") "Year" else "Location"
  pred <- bind_rows(lapply(unique(base[[fold_var]]), function(g) {
    train <- base[base[[fold_var]] != g, ]
    test  <- base[base[[fold_var]] == g, ]
    pc <- predict_centre(spec, train, test)
    wf <- stats::lm(width90_adj ~ log_total, data = train)
    tibble::tibble(flightID = test$flightID, Year = test$Year, obs = test$p50,
                   season_total = test$season_total,
                   mu = pc$mu, sigma = pc$sigma,
                   w = pmax(unname(stats::predict(wf, newdata = test)), 1))
  }))
  sn <- pred %>%
    inner_join(nights, by = "flightID", relationship = "one-to-many") %>%
    group_by(flightID) %>%
    mutate(share = night_shares(doy, first(mu), first(w), "logistic",
                                sigma = first(sigma), eps = eps)) %>%
    ungroup()
  ## Within-year skill: predicted and observed departures from their year means.
  wy <- pred %>% group_by(Year) %>% filter(n() > 1) %>%
    mutate(mu_W = mu - mean(mu), obs_W = obs - mean(obs)) %>% ungroup()
  tibble::tibble(
    spec = name, scheme = scheme, n = nrow(pred),
    log_score = sum(sn$moths * log(sn$share)) / sum(sn$moths),
    mae = mean(abs(pred$mu - pred$obs)),
    rmse = sqrt(mean((pred$mu - pred$obs)^2)),
    mae_moth = stats::weighted.mean(abs(pred$mu - pred$obs), pred$season_total),
    r_within = stats::cor(wy$mu_W, wy$obs_W),
    sigma = mean(pred$sigma)
  )
}

scores <- bind_rows(lapply(names(specs), function(sp)
  bind_rows(lapply(c("leave-one-year-out", "leave-one-locality-out"),
                   function(sc) run_spec(sp, sc))))) %>%
  group_by(scheme) %>%
  mutate(gain = log_score - log_score[spec == "adopted (xY + xW)"]) %>%
  ungroup() %>%
  arrange(scheme, desc(log_score))


## ---- 7. Coefficients of the local terms (random year intercept) ------------------------------

coefs <- bind_rows(lapply(names(specs)[-1], function(sp) {
  s <- specs[[sp]]
  if (length(s$terms) == 0) return(NULL)
  f <- stats::as.formula(paste("p50 ~", paste(c("xY", "xW", s$terms), collapse = " + "),
                               "+ (1 | Year)", if (s$re != "none") "+ (1 | Location)" else ""))
  m  <- lmer(f, data = base, REML = TRUE)
  fe <- lme4::fixef(m)
  se <- sqrt(diag(as.matrix(stats::vcov(m))))
  k  <- names(fe) %in% s$terms
  tibble::tibble(spec = sp, term = names(fe)[k], estimate = unname(fe[k]),
                 lo = unname(fe[k] - 1.96 * se[k]), hi = unname(fe[k] + 1.96 * se[k]))
}))

## Habitual lateness: locality intercepts on the full data (net of year
## errors), and how much of their variance geography explains (localities
## with >= 3 years).
m_hab <- lmer(p50 ~ xY + xW + (1 | Year) + (1 | Location), data = base, REML = TRUE)
hab <- tibble::tibble(Location = rownames(lme4::ranef(m_hab)$Location),
                      lateness = lme4::ranef(m_hab)$Location[, 1]) %>%
  left_join(base %>% count(Location, name = "n_years"), by = "Location") %>%
  left_join(coords, by = "Location") %>% left_join(elev, by = "Location")
hab3 <- filter(hab, n_years >= 3)
hab_geo_r2 <- summary(stats::lm(lateness ~ lat + lon + elev_m, data = hab3))$r.squared


## ---- 8. Write and report --------------------------------------------------------------------

summary_list <- list(
  n = nrow(base), n_years = n_distinct(base$Year), n_loc = n_distinct(base$Location),
  n_loc_multi = sum(table(base$Location) >= 2), n_boot = n_boot,
  samp_sd = sqrt(samp_var), resid_sd = stats::sigma(centre_lm),
  resid_W_sd = sqrt(s2), same_loc_r = same_loc_r, n_same_loc_pairs = nrow(same_loc),
  hab_sd = sqrt(vc_of(m_hab)[["Location"]]), hab_geo_r2 = hab_geo_r2,
  n_hab3 = nrow(hab3), t_windows = t_windows, gradients = grad_summary,
  temp_W_sd = sapply(c("T_MayJun_W", "T_JunJul_W"), function(v) stats::sd(base[[v]]))
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(data = base, variance = variance, correlogram = correlogram,
             gradients = gradients, scores = scores, coefs = coefs, lateness = hab, summary = summary_list),
        file.path(out_dir, "local_component.rds"))

cat(sprintf("Local component: %d locality-years, %d years, %d localities (%d with >= 2 years).\n",
            summary_list$n, summary_list$n_years, summary_list$n_loc, summary_list$n_loc_multi))
cat(sprintf("  centre residual SD %.2f d; within-year %.2f d; sampling SD of the median %.2f d\n",
            summary_list$resid_sd, summary_list$resid_W_sd, summary_list$samp_sd))
cat("  variance components:\n")
for (i in seq_len(nrow(variance))) cat(sprintf(
  "    %-16s %-12s var %6.2f [%6.2f, %6.2f] | SD %4.2f | share %4.2f\n",
  variance$model[i], variance$component[i], variance$var[i], variance$lo[i],
  variance$hi[i], variance$sd[i], variance$share[i]))
cat("  correlogram of within-year residuals:\n")
for (i in seq_len(nrow(correlogram))) cat(sprintf(
  "    %-10s pairs %4d | r %+5.2f [%+5.2f, %+5.2f] | null %+5.2f\n",
  as.character(correlogram$bin[i]), correlogram$n_pairs[i], correlogram$r[i],
  correlogram$lo[i], correlogram$hi[i], correlogram$null[i]))
cat(sprintf("    same locality, other years: r %+5.2f (%d pairs)\n", same_loc_r, nrow(same_loc)))
g <- grad_summary
cat(sprintf("  year-specific planes (%d years, %d records): share %.2f (null %.2f, 95%% %.2f, p %.3f)\n",
            g$n_years, g$n, g$share, g$null_mean, g$null_q95, g$p))
cat(sprintf("    latitude  d/deg: observed %+5.2f (SD %.2f, net %.2f) | BioSIM %+5.2f (SD %.2f) | r %+5.2f\n",
            g$obs_lat_mean, g$obs_lat_sd, g$obs_lat_sd_true, g$biosim_lat_mean, g$biosim_lat_sd, g$r_lat))
cat(sprintf("    longitude d/deg: observed %+5.2f (SD %.2f, net %.2f) | BioSIM %+5.2f (SD %.2f) | r %+5.2f\n",
            g$obs_lon_mean, g$obs_lon_sd, g$obs_lon_sd_true, g$biosim_lon_mean, g$biosim_lon_sd, g$r_lon))
cat("  scores:\n")
for (i in seq_len(nrow(scores))) cat(sprintf(
  "    %-24s %-36s log score %7.3f (%+6.3f) | MAE %4.2f | per moth %4.2f | RMSE %4.2f | r within %+5.2f | sigma %4.2f\n",
  scores$scheme[i], scores$spec[i], scores$log_score[i], scores$gain[i], scores$mae[i], scores$mae_moth[i],
  scores$rmse[i], scores$r_within[i], scores$sigma[i]))
cat("  coefficients (random year intercept):\n")
for (i in seq_len(nrow(coefs))) cat(sprintf(
  "    %-28s %-12s %+6.3f [%+6.3f, %+6.3f]\n", coefs$spec[i], coefs$term[i],
  coefs$estimate[i], coefs$lo[i], coefs$hi[i]))
cat(sprintf("  habitual lateness SD %.2f d; R2 on lat + lon + elevation %.2f (%d localities with >= 3 years)\n",
            summary_list$hab_sd, hab_geo_r2, nrow(hab3)))
cat("  written to data/processed/seasonal/local_component.rds\n")
