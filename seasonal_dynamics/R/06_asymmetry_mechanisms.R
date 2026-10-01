# =============================================================================
# seasonal_dynamics/R/06_asymmetry_mechanisms.R
# -----------------------------------------------------------------------------
# Why does a locality that flies later than the others in its year have a
# left-skewed season (05)? Two explanations, plus the null, make different
# predictions about WHICH part of the season is shared among the localities of
# a year.
#
#   pure shift      each locality's whole season moves with its own timing:
#                   every percentile varies alike among localities, and every
#                   percentile follows the locality's lateness alike.
#   common end      flight stops at about the same date everywhere in a year
#                   (a regional cool spell, say): the late percentiles
#                   (p90, p95) vary LESS among localities than the early ones
#                   and follow the locality's lateness less.
#   regional mixing part of every catch is moths flying across the region, so
#                   every locality's curve is pulled toward the regional curve
#                   in BOTH tails: early and late percentiles vary less than
#                   the median and follow the locality's lateness less, the
#                   centre most.
#
# TWO PROFILES ACROSS PERCENTILES q = p05 ... p95
#   1. Within-year spread. For every year with at least `min_per_year`
#      localities, each percentile is expressed as a departure from that
#      year's mean; the pooled within-year SD is computed per percentile and
#      corrected for sampling error (the mean resampling variance of that
#      percentile is subtracted, see below). Spread relative to the median.
#   2. Slope on lateness. Lateness of a locality is its mean within-year
#      departure of p50 over its OTHER years (leave-one-year-out), so it
#      carries no sampling error of the season being explained. The
#      within-year departure of each percentile is regressed on it.
#   Both profiles get 95% intervals from a bootstrap over YEARS (years are the
#   independent units: localities of one year share its weather).
#
# SAMPLING CORRECTION
#   Each record's catch is redrawn `n_boot` times at its own size from its own
#   nightly shares (multinomial) and its percentiles recomputed; the variance
#   of those draws is the record's sampling variance for each percentile.
#   Tail percentiles rest on fewer moths and are noisier, which would inflate
#   their spread and bias the test against a common end.
#
# ANALYSIS SET   as in 05: usable for phenology, at least 3 positive nights.
#
# OUTPUT  data/processed/seasonal/asymmetry_mechanisms.rds
#         list(data, spread, slopes, summary)
#
# REQUIREMENTS  dplyr
#
# RUNTIME  about two minutes.
#
# USAGE  source("seasonal_dynamics/R/06_asymmetry_mechanisms.R")   # from the project root
# =============================================================================

library(dplyr)

source(file.path("R", "functions", "phenology.R"))

out_dir      <- file.path("data", "processed", "seasonal")
probs        <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)
qn           <- sprintf("p%02d", round(probs * 100))
min_per_year <- 3
n_boot       <- 200          # resampling draws per record
n_yboot      <- 1000         # bootstrap over years
set.seed(1978)


## ---- 1. Data -------------------------------------------------------------------

season <- readRDS(file.path(out_dir, "season_percentiles.rds")) %>%
  filter(usable_phenology, n_positive >= 3) %>%
  group_by(year) %>%
  filter(n() >= min_per_year) %>%
  ungroup()

nights <- readRDS(file.path("data", "processed", "trap_daily.rds")) %>%
  filter(status == "counted", flightID %in% season$flightID) %>%
  arrange(flightID, doy)


## ---- 2. Sampling variance of every percentile, per record --------------------------

samp_var <- bind_rows(lapply(seq_len(nrow(season)), function(i) {
  n_id  <- nights[nights$flightID == season$flightID[i], ]
  draws <- stats::rmultinom(n_boot, round(season$season_total[i]),
                            n_id$moths / sum(n_id$moths))
  ## Not named `q`: inside tibble() a column called q would shadow it.
  pct <- apply(draws, 2, function(x) unname(curve_percentiles(n_id$doy, x, probs)))
  tibble::tibble(flightID = season$flightID[i], q = qn,
                 samp_var = apply(pct, 1, stats::var))
}))


## ---- 3. Within-year departures and leave-one-year-out lateness --------------------

long <- season %>%
  select(flightID, locality, year, all_of(qn)) %>%
  tidyr::pivot_longer(all_of(qn), names_to = "q", values_to = "doy") %>%
  group_by(year, q) %>%
  mutate(dep = doy - mean(doy)) %>%
  ungroup() %>%
  left_join(samp_var, by = c("flightID", "q"))

## Lateness: the locality's mean p50 departure in its other years.
p50_dep <- long %>% filter(q == "p50") %>% select(flightID, locality, year, dep)
lateness <- p50_dep %>%
  group_by(locality) %>%
  mutate(n_loc = n(),
         lateness = if_else(n_loc > 1, (sum(dep) - dep) / (n_loc - 1), NA_real_)) %>%
  ungroup() %>%
  select(flightID, lateness)
long <- left_join(long, lateness, by = "flightID")


## ---- 4. The two profiles, and their bootstrap over years --------------------------

## Pooled within-year variance: departures are from the year mean, so each
## year loses one degree of freedom.
profiles <- function(d) {
  df_within <- d %>% filter(q == "p50") %>% count(year) %>%
    summarise(df = sum(n - 1)) %>% pull(df)
  spread <- d %>%
    group_by(q) %>%
    summarise(var_obs  = sum(dep^2) / df_within,
              var_samp = mean(samp_var),
              .groups  = "drop") %>%
    mutate(sd_obs  = sqrt(var_obs),
           sd_corr = sqrt(pmax(var_obs - var_samp, 0)))
  spread <- spread %>%
    mutate(rel_obs  = sd_obs  / sd_obs[q == "p50"],
           rel_corr = sd_corr / sd_corr[q == "p50"])
  slopes <- d %>%
    filter(!is.na(lateness)) %>%
    group_by(q) %>%
    summarise(slope = unname(stats::coef(stats::lm(dep ~ lateness))[2]),
              .groups = "drop") %>%
    mutate(rel_slope = slope / slope[q == "p50"])
  list(spread = spread, slopes = slopes)
}

point <- profiles(long)

years <- unique(long$year)
boot <- lapply(seq_len(n_yboot), function(b) {
  ys <- sample(years, replace = TRUE)
  ## Re-label resampled years so that a year drawn twice counts as two years.
  d <- bind_rows(lapply(seq_along(ys), function(k)
    mutate(long[long$year == ys[k], ], year = k)))
  p <- profiles(d)
  list(spread = mutate(p$spread, b = b), slopes = mutate(p$slopes, b = b))
})
ci <- function(x) stats::quantile(x, c(0.025, 0.975), na.rm = TRUE, names = FALSE)

spread <- point$spread %>%
  left_join(bind_rows(lapply(boot, `[[`, "spread")) %>%
              group_by(q) %>%
              summarise(rel_corr_lo = ci(rel_corr)[1], rel_corr_hi = ci(rel_corr)[2],
                        sd_corr_lo = ci(sd_corr)[1], sd_corr_hi = ci(sd_corr)[2],
                        .groups = "drop"),
            by = "q") %>%
  mutate(q = factor(q, qn)) %>% arrange(q)

slopes <- point$slopes %>%
  left_join(bind_rows(lapply(boot, `[[`, "slopes")) %>%
              group_by(q) %>%
              summarise(slope_lo = ci(slope)[1], slope_hi = ci(slope)[2],
                        rel_lo = ci(rel_slope)[1], rel_hi = ci(rel_slope)[2],
                        .groups = "drop"),
            by = "q") %>%
  mutate(q = factor(q, qn)) %>% arrange(q)

## The contrast that separates the hypotheses: early tail against late tail,
## each relative to the centre. Common end -> late < early; regional mixing ->
## both < 1 and similar; pure shift -> both ~ 1.
contrast <- function(tab, col) {
  v <- setNames(tab[[col]], as.character(tab$q))
  c(early = mean(v[c("p05", "p10")]), late = mean(v[c("p90", "p95")]))
}
boot_contrast <- bind_rows(lapply(boot, function(x) {
  s <- arrange(mutate(x$spread, q = factor(q, qn)), q)
  l <- arrange(mutate(x$slopes, q = factor(q, qn)), q)
  tibble::tibble(spread_early = contrast(s, "rel_corr")[["early"]],
                 spread_late  = contrast(s, "rel_corr")[["late"]],
                 slope_early  = contrast(l, "rel_slope")[["early"]],
                 slope_late   = contrast(l, "rel_slope")[["late"]])
})) %>%
  mutate(spread_diff = spread_late - spread_early,
         slope_diff  = slope_late - slope_early)

contrasts <- tibble::tibble(
  quantity = c("spread, early tail / median", "spread, late tail / median",
               "spread, late minus early",
               "slope, early tail / median", "slope, late tail / median",
               "slope, late minus early"),
  estimate = c(contrast(spread, "rel_corr"),
               diff(rev(contrast(spread, "rel_corr")))  * -1,
               contrast(slopes, "rel_slope"),
               diff(rev(contrast(slopes, "rel_slope"))) * -1),
  lo = c(ci(boot_contrast$spread_early)[1], ci(boot_contrast$spread_late)[1],
         ci(boot_contrast$spread_diff)[1],
         ci(boot_contrast$slope_early)[1],  ci(boot_contrast$slope_late)[1],
         ci(boot_contrast$slope_diff)[1]),
  hi = c(ci(boot_contrast$spread_early)[2], ci(boot_contrast$spread_late)[2],
         ci(boot_contrast$spread_diff)[2],
         ci(boot_contrast$slope_early)[2],  ci(boot_contrast$slope_late)[2],
         ci(boot_contrast$slope_diff)[2])
)


## ---- 5. Write and report -------------------------------------------------------------

summary_list <- list(
  n = nrow(season), n_years = n_distinct(season$year),
  n_lateness = sum(!is.na(lateness$lateness)),
  n_localities_lateness = n_distinct(long$locality[!is.na(long$lateness)]),
  n_boot = n_boot, n_yboot = n_yboot, contrasts = contrasts
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(data = long, spread = spread, slopes = slopes,
             contrasts = contrasts, summary = summary_list),
        file.path(out_dir, "asymmetry_mechanisms.rds"))

cat(sprintf("Asymmetry mechanisms: %d locality-years in %d years (lateness available for %d).\n",
            summary_list$n, summary_list$n_years, summary_list$n_lateness))
cat("  within-year SD by percentile (observed -> sampling-corrected, relative to p50):\n")
for (i in seq_len(nrow(spread))) cat(sprintf(
  "    %s  sd %5.2f -> %5.2f d | relative %4.2f [%4.2f, %4.2f]\n",
  spread$q[i], spread$sd_obs[i], spread$sd_corr[i], spread$rel_corr[i],
  spread$rel_corr_lo[i], spread$rel_corr_hi[i]))
cat("  slope of each percentile on the locality's lateness (other years):\n")
for (i in seq_len(nrow(slopes))) cat(sprintf(
  "    %s  slope %5.2f [%5.2f, %5.2f] | relative to p50 %4.2f [%4.2f, %4.2f]\n",
  slopes$q[i], slopes$slope[i], slopes$slope_lo[i], slopes$slope_hi[i],
  slopes$rel_slope[i], slopes$rel_lo[i], slopes$rel_hi[i]))
cat("  contrasts (bootstrap over years):\n")
for (i in seq_len(nrow(contrasts))) cat(sprintf(
  "    %-30s %+5.2f [%+5.2f, %+5.2f]\n", contrasts$quantity[i],
  contrasts$estimate[i], contrasts$lo[i], contrasts$hi[i]))
cat("  written to data/processed/seasonal/asymmetry_mechanisms.rds\n")
