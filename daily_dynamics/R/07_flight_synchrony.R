# =============================================================================
# 07_flight_synchrony.R
# -----------------------------------------------------------------------------
# Is the observed flight synchronous across Maine, and are southern catches in
# synchrony with northern ones?
#
# WHY THIS TEST
#   The comparison in R/05 leaves a spatial residual: both models run early in
#   the south, by up to a fortnight at some localities. Two explanations remain
#   after R/06 ruled out a temperature-driven pupal stage -- a site-level
#   effect neither model contains, and immigration, since a locality's catch
#   need not be its own cohort. This script tests the second directly, without
#   using weather at all, by asking whether southern catches happen on the
#   region's dates or on their own.
#
#   Two predictions separate them:
#     seasonal scale  if southern catch is the local cohort, the observed
#                     dates should carry the latitudinal gradient the models
#                     predict; immigration from the north would flatten it.
#     nightly scale   a wave of arriving moths should make distant localities
#                     spike on the same calendar nights; locally emerged
#                     moths, gated by local weather, should not.
#
# WHAT IS COMPUTED
#   1. days per degree of latitude in observed and modelled dates, with year
#      fixed effects so the shared interannual signal is removed;
#   2. the south-minus-north gap in each year, observed against modelled;
#   3. the correlation of nightly catch shares between every pair of
#      localities trapped in the same year, against their separation, with a
#      null obtained by pairing the same two localities in different years.
#
#   The null matters: nightly catch curves are spiky, and two unrelated spiky
#   series can correlate by chance. Shuffling the year keeps each locality's
#   own seasonal window and spikiness and destroys only the shared nights.
#
# INPUT   data/processed/model_vs_observed.rds, trap_daily.rds
# OUTPUT  data/processed/flight_synchrony.rds
#           list(gradients, gaps, pairs, null_draws, summary)
#
# USAGE  source("daily_dynamics/R/07_flight_synchrony.R")   # from the project root
# =============================================================================

library(dplyr)

## ---- 1. Configuration -------------------------------------------------------

## Latitude splitting the trap network into two halves of comparable size.
south_of <- 45.5

## Nightly synchrony: lags searched, and the minimum number of nights both
## localities were trapped for a pair to be used.
lag_window   <- -7:7
min_overlap  <- 10

## Draws of the shuffled-year null.
n_null_draws <- 20
seed         <- 1978

out_file <- file.path("data", "processed", "flight_synchrony.rds")


## ---- 2. Data ----------------------------------------------------------------

cmp <- readRDS(file.path("data", "processed", "model_vs_observed.rds")) %>%
  filter(comparable)

## Nightly catch as a share of the locality-year's own season total, so that a
## correlation compares shape rather than size.
nightly <- readRDS(file.path("data", "processed", "trap_daily.rds")) %>%
  filter(status == "counted") %>%
  semi_join(cmp, by = c("locality" = "Location", "year" = "Year")) %>%
  group_by(locality, year) %>%
  mutate(share = moths / sum(moths)) %>%
  ungroup()

sites <- distinct(cmp, Location, Latitude, Longitude)


## ---- 3. The latitudinal gradient --------------------------------------------
## Year fixed effects remove the interannual signal common to all localities,
## leaving the gradient across localities within years.

gradient_of <- function(response, data, label) {
  f <- lm(reformulate(c("Latitude", "factor(Year)"), response), data = data)
  ci <- confint(f)["Latitude", ]
  tibble::tibble(series = label, n = nrow(data),
                 d_per_degree = coef(f)[["Latitude"]],
                 lo = ci[[1]], hi = ci[[2]])
}

gradients <- bind_rows(
  gradient_of("obs_p50", cmp, "Observed median catch"),
  gradient_of("sim_p50", cmp, "BioSIM modelled flight"),
  gradient_of("bayes_pupation_p50", cmp, "bayessbw pupation"),
  gradient_of("obs_p50", filter(cmp, season_total >= 1000),
              "Observed, catches of 1,000 or more"),
  gradient_of("obs_p50", filter(cmp, season_total >= 5000),
              "Observed, catches of 5,000 or more"),
  gradient_of("sim_p50", filter(cmp, season_total >= 1000),
              "BioSIM, catches of 1,000 or more")
)

lat_span <- diff(range(cmp$Latitude))


## ---- 4. South against north, year by year -----------------------------------

gaps <- cmp %>%
  mutate(region = ifelse(Latitude < south_of, "south", "north")) %>%
  group_by(Year, region) %>%
  summarise(n = n(), obs = mean(obs_p50), sim = mean(sim_p50), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = region, values_from = c(n, obs, sim)) %>%
  filter(!is.na(n_north), !is.na(n_south), n_north >= 2, n_south >= 2) %>%
  mutate(obs_gap = obs_south - obs_north,
         sim_gap = sim_south - sim_north)


## ---- 5. Nightly synchrony between pairs of localities -----------------------

## Correlation of two nightly share series on the nights both were trapped,
## with `b` shifted by `lag` days.
pair_cor <- function(a, b, lag = 0) {
  bb_doy <- b$doy + lag
  cm <- intersect(a$doy, bb_doy)
  if (length(cm) < min_overlap) return(NA_real_)
  s1 <- a$share[match(cm, a$doy)]
  s2 <- b$share[match(cm, bb_doy)]
  if (sd(s1) == 0 || sd(s2) == 0) return(NA_real_)
  cor(s1, s2)
}

series_of <- function(loc, yr) {
  x <- nightly[nightly$locality == loc & nightly$year == yr, c("doy", "share")]
  x[order(x$doy), ]
}

region_of <- function(loc) ifelse(sites$Latitude[match(loc, sites$Location)] < south_of,
                                  "S", "N")

## Great-circle distance in km.
dist_km <- function(la, lo, lb, lob) {
  rad <- pi / 180
  6371 * acos(pmin(1, sin(la * rad) * sin(lb * rad) +
                     cos(la * rad) * cos(lb * rad) * cos((lob - lo) * rad)))
}

pairs_in_year <- function(y) {
  locs <- sort(unique(nightly$locality[nightly$year == y]))
  if (length(locs) < 2) return(NULL)
  grid <- t(utils::combn(locs, 2))
  out <- vector("list", nrow(grid))
  for (k in seq_len(nrow(grid))) {
    a <- series_of(grid[k, 1], y); b <- series_of(grid[k, 2], y)
    r0 <- pair_cor(a, b)
    if (is.na(r0)) next
    cc <- vapply(lag_window, function(L) pair_cor(a, b, L), numeric(1))
    ia <- match(grid[k, 1], sites$Location); ib <- match(grid[k, 2], sites$Location)
    out[[k]] <- tibble::tibble(
      Year = y, a = grid[k, 1], b = grid[k, 2], r0 = r0,
      r_best   = if (all(is.na(cc))) NA_real_ else max(cc, na.rm = TRUE),
      lag_best = if (all(is.na(cc))) NA_real_ else
        lag_window[which.max(replace(cc, is.na(cc), -Inf))],
      dlat    = abs(sites$Latitude[ia] - sites$Latitude[ib]),
      dist_km = dist_km(sites$Latitude[ia], sites$Longitude[ia],
                        sites$Latitude[ib], sites$Longitude[ib]),
      pair_type = paste(sort(c(region_of(grid[k, 1]), region_of(grid[k, 2]))),
                        collapse = "")
    )
  }
  bind_rows(out)
}

pairs_df <- bind_rows(lapply(sort(unique(nightly$year)), pairs_in_year))

## Null: the same two localities, but the second one in a different year.
set.seed(seed)
null_draws <- vapply(seq_len(n_null_draws), function(draw) {
  r <- vapply(seq_len(nrow(pairs_df)), function(i) {
    other <- setdiff(unique(nightly$year[nightly$locality == pairs_df$b[i]]),
                     pairs_df$Year[i])
    if (!length(other)) return(NA_real_)
    pair_cor(series_of(pairs_df$a[i], pairs_df$Year[i]),
             series_of(pairs_df$b[i], other[sample.int(length(other), 1)]))
  }, numeric(1))
  c(median = median(r, na.rm = TRUE), above_half = mean(r > 0.5, na.rm = TRUE))
}, numeric(2))


## ---- 6. Write and report ----------------------------------------------------

summary_list <- list(
  n_rows        = nrow(cmp),
  lat_span      = lat_span,
  obs_spread    = gradients$d_per_degree[1] * lat_span,
  sim_spread    = gradients$d_per_degree[2] * lat_span,
  bayes_spread  = gradients$d_per_degree[3] * lat_span,
  obs_gap       = mean(gaps$obs_gap),
  sim_gap       = mean(gaps$sim_gap),
  n_pairs       = nrow(pairs_df),
  r0_median     = median(pairs_df$r0),
  r0_above_half = mean(pairs_df$r0 > 0.5),
  null_median   = mean(null_draws["median", ]),
  null_above_half = mean(null_draws["above_half", ]),
  r_dist        = cor(pairs_df$r0, pairs_df$dist_km),
  r_dlat        = cor(pairs_df$r0, pairs_df$dlat)
)

dir.create(file.path("data", "processed"), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(gradients = gradients, gaps = gaps, pairs = pairs_df,
             null_draws = null_draws, summary = summary_list),
        out_file)

message(sprintf(
  paste0("Flight synchrony over %d comparable locality-years, %.1f degrees of latitude.\n",
         "  days per degree: observed %+.2f | BioSIM %+.2f | bayessbw %+.2f\n",
         "  implied spread across the network: observed %.1f d, BioSIM %.1f d, bayessbw %.1f d\n",
         "  south minus north: observed %+.1f d, modelled %+.1f d (%d years)\n",
         "  nightly synchrony, %d locality pairs: median r %.2f against a null of %.2f;\n",
         "    %.0f%% above 0.5 against %.0f%% in the null; r falls with distance (%+.2f)\n",
         "  written to %s"),
  summary_list$n_rows, lat_span,
  gradients$d_per_degree[1], gradients$d_per_degree[2], gradients$d_per_degree[3],
  summary_list$obs_spread, summary_list$sim_spread, summary_list$bayes_spread,
  summary_list$obs_gap, summary_list$sim_gap, nrow(gaps),
  summary_list$n_pairs, summary_list$r0_median, summary_list$null_median,
  100 * summary_list$r0_above_half, 100 * summary_list$null_above_half,
  summary_list$r_dist, out_file
))
