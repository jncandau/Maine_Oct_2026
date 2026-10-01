# =============================================================================
# seasonal_dynamics/R/02_pupation_dates.R
# -----------------------------------------------------------------------------
# The predicted pupation date of every trap locality-year, from the two
# development models of the shared data layer, as distributions comparable
# with each other and with the observed catch percentiles.
#
# INPUT   data/processed/biosim_phenology.rds  (R/03_biosim_phenology.R:
#                                           series Pupation_mf and Pupae_mf)
#         data/interim/biosim_raw.rds      (R/03: daily output, for the
#                                           vanished-pupae check only)
#         data/processed/pupal_phenology.rds   (R/04_pupal_phenology.R)
#
# BIOSIM PUPATION IS DERIVED (in R/03), NOT A MODEL OUTPUT
#   BioSIM reports pupae as a STOCK -- the share of the cohort that is a pupa
#   on a given day -- not as a distribution of pupation dates. The median of
#   that stock curve (the `Pupae_mf` series of R/03, used as "BioSIM pupation"
#   in the daily analysis) is the middle of the pupal period, roughly the
#   pupation date plus half the pupal stage. The pupation-date distribution is
#   recovered from mass balance: with mortality switched off (as in R/03), an
#   individual that has pupated by day t is either still a pupa or has already
#   emerged as an adult, so
#
#       cumulative pupation(t) = pupae(t) + cumulative emergence(t)
#
#   with pupae = MalePupae + FemalePupae and emergence = MaleEmergence +
#   FemaleEmergence, all in percent of the cohort. Its daily increments are the
#   pupation-date distribution, and its percentiles are taken with the shared
#   helper. Runs (seeds) are averaged first, as in R/03.
#
# BAYESSBW PUPATION
#   R/04 already reports pupation-date percentiles of the median simulated
#   population (`p05` ... `p95`), which are used as they are.
#
# PUPAE THAT VANISH
#   In a few locality-years BioSIM's pupae drop to zero within a day without
#   emerging, although mortality is switched off; the cases coincide with the
#   incomplete cohorts flagged in R/03 (Clayton Lake recurs), and the
#   mechanism inside the model is not identified. On such a day the mass
#   balance goes negative. Pupation ENTRIES are unaffected on every other day,
#   so the percentiles are taken over the positive increments only, and the
#   loss is carried as `biosim_pup_lost` (percent of the cohort).
#
# COMPLETENESS
#   `biosim_pup_complete` is TRUE when the pupation entries add up to at least
#   99% of the cohort; a percentile of an incomplete cohort is the percentile
#   of the part that developed.
#
# OUTPUT  data/processed/seasonal/pupation_dates.rds -- one row per
#         locality-year: biosim_pup_p05/p50/p95, biosim_pup_complete,
#         biosim_pup_lost,
#         biosim_stock_p50 (the stock median, for comparison),
#         bayes_pup_p05/p50/p95
#
# REQUIREMENTS  dplyr
#
# USAGE  source("seasonal_dynamics/R/02_pupation_dates.R")   # from the project root
# =============================================================================

library(dplyr)

source(file.path("R", "functions", "phenology.R"))

out_dir <- file.path("data", "processed", "seasonal")


## ---- 1. BioSIM ------------------------------------------------------------------
## The pupation-date distribution is built once, in the shared data layer
## (R/03, series `Pupation_mf`, by the mass balance described above); its
## percentiles and its total are read from there, and the stock median
## (series `Pupae_mf`) is kept for comparison.

pheno <- readRDS(file.path("data", "processed", "biosim_phenology.rds"))

biosim_pup <- pheno %>%
  filter(series == "Pupation_mf") %>%
  select(Location, Year, KeyID, biosim_pup_p05 = p05, biosim_pup_p50 = p50,
         biosim_pup_p95 = p95, biosim_pup_final = season_pct) %>%
  left_join(pheno %>% filter(series == "Pupae_mf") %>%
              select(Location, Year, biosim_stock_p50 = p50),
            by = c("Location", "Year"))

## The vanished pupae need the raw series: on a day when pupae disappear
## without emerging the mass balance falls, and R/03 drops that increment.
raw <- readRDS(file.path("data", "interim", "biosim_raw.rds"))

loss <- raw %>%
  mutate(doy = as.integer(format(as.Date(sprintf("%d-%02d-%02d", Year, Month, Day)), "%j")),
         pupae = MalePupae + FemalePupae,
         emergence = MaleEmergence + FemaleEmergence) %>%
  group_by(KeyID, Year, doy) %>%                     # average over runs
  summarise(pupae = mean(pupae), emergence = mean(emergence), .groups = "drop") %>%
  arrange(KeyID, Year, doy) %>%
  group_by(KeyID, Year) %>%
  mutate(cum = pupae + cumsum(emergence), inc = c(cum[1], diff(cum))) %>%
  summarise(biosim_pup_lost = -sum(pmin(inc, 0)), .groups = "drop")
n_loss <- sum(loss$biosim_pup_lost > 0.5)

biosim_pup <- biosim_pup %>%
  left_join(loss, by = c("KeyID", "Year")) %>%
  mutate(biosim_pup_complete = biosim_pup_final >= 99) %>%
  select(Location, Year, biosim_pup_p05, biosim_pup_p50, biosim_pup_p95,
         biosim_pup_final, biosim_pup_lost, biosim_pup_complete,
         biosim_stock_p50)


## ---- 2. bayessbw -----------------------------------------------------------------

bayes_pup <- readRDS(file.path("data", "processed", "pupal_phenology.rds")) %>%
  select(Location, Year, bayes_pup_p05 = p05, bayes_pup_p50 = p50,
         bayes_pup_p95 = p95)


## ---- 3. Join, write, report ----------------------------------------------------

pupation <- full_join(biosim_pup, bayes_pup, by = c("Location", "Year")) %>%
  arrange(Year, Location)

stopifnot(!any(duplicated(pupation[c("Location", "Year")])))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(pupation, file.path(out_dir, "pupation_dates.rds"))

both <- filter(pupation, !is.na(biosim_pup_p50), !is.na(bayes_pup_p50))
cat(sprintf(
  paste0("Pupation dates: %d locality-years (%d with both models; BioSIM cohort ",
         "fully pupated in %d; pupae vanish in %d).\n",
         "  BioSIM pupation p50: median DOY %.1f | stock median: %.1f ",
         "(stock - pupation: median %+.1f d)\n",
         "  bayessbw pupation p50: median DOY %.1f | BioSIM - bayessbw: ",
         "mean %+.1f d, r = %.2f\n",
         "  written to %s\n"),
  nrow(pupation), nrow(both), sum(pupation$biosim_pup_complete, na.rm = TRUE),
  n_loss,
  median(both$biosim_pup_p50), median(both$biosim_stock_p50),
  median(both$biosim_stock_p50 - both$biosim_pup_p50),
  median(both$bayes_pup_p50),
  mean(both$biosim_pup_p50 - both$bayes_pup_p50),
  cor(both$biosim_pup_p50, both$bayes_pup_p50),
  file.path(out_dir, "pupation_dates.rds")))
