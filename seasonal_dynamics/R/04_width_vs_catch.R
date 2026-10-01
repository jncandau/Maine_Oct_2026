# =============================================================================
# seasonal_dynamics/R/04_width_vs_catch.R
# -----------------------------------------------------------------------------
# Is the width of the flight season a function of the total trap catch?
#
# INPUT   data/processed/seasonal/season_percentiles.rds  (01)
#         data/processed/trap_daily.rds                    (R/01: nightly
#                                                           counts, for the
#                                                           rarefaction null)
#
# THE QUESTION HAS TWO PARTS
#   1. Is there a relationship? width90 (p95 - p05) and width50 (p75 - p25)
#      are regressed on log10(season total), alone, with year fixed effects,
#      and split into a between-locality and a within-locality component (a
#      locality's large-catch years against its small-catch years).
#   2. Is it biology or sampling? A percentile estimated from 150 moths is
#      noisier than one estimated from 15,000, and with a skewed, spiky curve
#      that noise need not average out to zero: a small catch can be
#      systematically narrower or wider than the season it samples. The
#      rarefaction null measures that directly. Every record with at least
#      `null_min_total` moths is treated as a known season; `null_reps`
#      catches of size N are drawn from it (multinomial over its nights, with
#      the observed nightly shares as probabilities) for each N in `null_sizes`,
#      and the width of each draw is computed exactly as for the real data.
#      The change in mean width between N and the full catch is the width
#      change that sampling alone produces. If the observed relationship is
#      no steeper than this, catch size affects the precision of the width,
#      not the season.
#
# ANALYSIS SET   records usable for phenology (R/01). A robustness fit drops
#                seasons with fewer than 3 positive nights.
#
# OUTPUT  data/processed/seasonal/width_vs_catch.rds
#         list(data, fits, decomposition, null, summary)
#
# REQUIREMENTS  dplyr
#
# RUNTIME  under a minute.
#
# USAGE  source("seasonal_dynamics/R/04_width_vs_catch.R")   # from the project root
# =============================================================================

library(dplyr)

source(file.path("R", "functions", "phenology.R"))

out_dir <- file.path("data", "processed", "seasonal")

null_min_total <- 3000                          # seasons treated as "known"
null_sizes     <- c(100, 200, 500, 1000, 2000)  # rarefied catch sizes
null_reps      <- 200
set.seed(1978)


## ---- 1. Data ---------------------------------------------------------------------

season <- readRDS(file.path(out_dir, "season_percentiles.rds")) %>%
  filter(usable_phenology) %>%
  mutate(log_total = log10(season_total))


## ---- 2. The relationship ----------------------------------------------------------

slope_of <- function(fit, term = "log_total") {
  ci <- stats::confint(fit)[term, ]
  tibble::tibble(slope = unname(stats::coef(fit)[term]),
                 lo = unname(ci[1]), hi = unname(ci[2]),
                 r2 = summary(fit)$r.squared)
}

fits <- bind_rows(lapply(c("width90", "width50"), function(w) {
  f <- stats::as.formula(paste(w, "~ log_total"))
  g <- stats::as.formula(paste(w, "~ log_total + factor(year)"))
  bind_rows(
    mutate(slope_of(stats::lm(f, data = season)), model = "alone"),
    mutate(slope_of(stats::lm(g, data = season)), model = "year fixed effects"),
    mutate(slope_of(stats::lm(f, data = filter(season, n_positive >= 3))),
           model = "alone, >= 3 positive nights")
  ) %>% mutate(width = w, n = c(nrow(season), nrow(season),
                                 sum(season$n_positive >= 3)))
}))

## Between and within localities: localities with at least 3 locality-years.
multi <- season %>% group_by(locality) %>% filter(n() >= 3) %>% ungroup()
decomposition <- bind_rows(lapply(c("width90", "width50"), function(w) {
  between <- multi %>% group_by(locality) %>%
    summarise(wid = mean(.data[[w]]), lt = mean(log_total), .groups = "drop")
  within  <- multi %>% group_by(locality) %>%
    mutate(wid = .data[[w]] - mean(.data[[w]]), lt = log_total - mean(log_total)) %>%
    ungroup()
  tibble::tibble(
    width = w,
    component = c("between localities", "within localities"),
    r = c(stats::cor(between$wid, between$lt), stats::cor(within$wid, within$lt)),
    slope = c(unname(stats::coef(stats::lm(wid ~ lt, between))[2]),
              unname(stats::coef(stats::lm(wid ~ lt, within))[2])),
    n = c(nrow(between), nrow(within)),
    n_localities = n_distinct(multi$locality)
  )
}))


## ---- 3. Rarefaction null ------------------------------------------------------------

nights <- readRDS(file.path("data", "processed", "trap_daily.rds")) %>%
  filter(status == "counted", flightID %in% season$flightID) %>%
  arrange(flightID, doy)

widths <- function(doy, moths) {
  q <- curve_percentiles(doy, moths, c(0.05, 0.25, 0.75, 0.95))
  c(width90 = unname(q[4] - q[1]), width50 = unname(q[3] - q[2]))
}

known <- season %>% filter(season_total >= null_min_total)

null <- bind_rows(lapply(known$flightID, function(id) {
  n_id <- filter(nights, flightID == id)
  prob <- n_id$moths / sum(n_id$moths)
  full <- widths(n_id$doy, n_id$moths)
  bind_rows(lapply(null_sizes, function(N) {
    draws <- stats::rmultinom(null_reps, N, prob)
    w <- apply(draws, 2, function(x) widths(n_id$doy, x))
    tibble::tibble(flightID = id, N = N,
                   d_width90 = mean(w["width90", ]) - full[["width90"]],
                   d_width50 = mean(w["width50", ]) - full[["width50"]],
                   sd_width90 = stats::sd(w["width90", ]),
                   sd_width50 = stats::sd(w["width50", ]))
  }))
}))

null_summary <- null %>%
  group_by(N) %>%
  summarise(across(c(d_width90, d_width50, sd_width90, sd_width50), mean),
            .groups = "drop")

## The width change sampling alone produces per tenfold increase in catch,
## over the rarefied range, to set against the observed slope.
null_slope <- c(
  width90 = unname(stats::coef(stats::lm(d_width90 ~ log10(N), null))[2]),
  width50 = unname(stats::coef(stats::lm(d_width50 ~ log10(N), null))[2])
)


## The sampling bias expected at each record's own catch size, interpolated on
## log10(N) from the null, taken as zero from 10,000 moths up (the null's
## reference seasons are themselves catches of 3,000 or more, so the bias
## between them and an infinite catch is not measured and is small). The width
## minus that bias is the width the season would show if catch size only
## affected precision; regressing it on catch size is the direct test.
bias_at <- function(N, d) {
  stats::approx(log10(c(null_sizes, 1e4)), c(d, 0), xout = log10(N), rule = 2)$y
}
season <- season %>%
  mutate(width90_adj = width90 - bias_at(season_total, null_summary$d_width90),
         width50_adj = width50 - bias_at(season_total, null_summary$d_width50))

fits <- bind_rows(fits, bind_rows(lapply(c("width90", "width50"), function(w) {
  f <- stats::as.formula(paste0(w, "_adj ~ log_total"))
  mutate(slope_of(stats::lm(f, data = season)),
         model = "alone, sampling bias removed", width = w, n = nrow(season))
})))


## ---- 4. Write and report -----------------------------------------------------------

summary_list <- list(
  n            = nrow(season),
  n_known      = nrow(known),
  null_sizes   = null_sizes,
  null_reps    = null_reps,
  null_slope   = null_slope,
  r            = c(width90 = stats::cor(season$width90, season$log_total),
                   width50 = stats::cor(season$width50, season$log_total))
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(data = season, fits = fits, decomposition = decomposition,
             null = null, null_summary = null_summary, summary = summary_list),
        file.path(out_dir, "width_vs_catch.rds"))

cat(sprintf("Season width vs total catch: %d locality-years.\n", nrow(season)))
for (i in seq_len(nrow(fits))) cat(sprintf(
  "    %-7s %-28s n %3d | %+5.2f d per tenfold [%+5.2f, %+5.2f] | R2 %.3f\n",
  fits$width[i], fits$model[i], fits$n[i], fits$slope[i], fits$lo[i], fits$hi[i],
  fits$r2[i]))
cat(sprintf("  between/within (%d localities with >= 3 years):\n",
            decomposition$n_localities[1]))
for (i in seq_len(nrow(decomposition))) cat(sprintf(
  "    %-7s %-19s r %+5.2f | slope %+5.2f d per tenfold (n %d)\n",
  decomposition$width[i], decomposition$component[i], decomposition$r[i],
  decomposition$slope[i], decomposition$n[i]))
cat(sprintf("  rarefaction null from %d seasons of >= %d moths, %d draws per size:\n",
            nrow(known), null_min_total, null_reps))
for (i in seq_len(nrow(null_summary))) cat(sprintf(
  "    N = %5d: width90 %+5.2f d (sd %4.2f) | width50 %+5.2f d (sd %4.2f)\n",
  null_summary$N[i], null_summary$d_width90[i], null_summary$sd_width90[i],
  null_summary$d_width50[i], null_summary$sd_width50[i]))
cat(sprintf("  null slope: width90 %+5.2f, width50 %+5.2f d per tenfold\n",
            null_slope[["width90"]], null_slope[["width50"]]))
cat("  written to data/processed/seasonal/width_vs_catch.rds\n")
