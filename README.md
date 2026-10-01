# Spruce budworm flight phenology in Maine

Daily light-trap catches of eastern spruce budworm (*Choristoneura fumiferana*)
moths in Maine, 1968-1989, analysed to characterise and predict adult flight
activity. Statistical analyses are in R and every step is documented in a
Quarto analysis log.

## Layout

The project holds two analyses of one dataset.

| Path | Contents |
|------|----------|
| `R/01_read_trap_data.R` | **Shared data layer.** Reads the Flight sheet, selects Maine light traps, applies the season-window, zero-fill and censoring rules, writes `data/processed/trap_records.rds` and `trap_daily.rds` |
| `R/02_read_weather_data.R` | **Shared data layer.** Hourly CaSR v3.2 weather for every trap locality-year, `data/processed/weather_hourly.rds` |
| `R/03_biosim_phenology.R` | **Shared data layer.** BioSIM *Spruce_Budworm_Biology* at every trap locality-year (mortality off), `data/processed/biosim_phenology.rds`, raw daily output in `data/interim/biosim_raw.rds` |
| `R/04_pupal_phenology.R` | **Shared data layer.** bayessbw pupation dates at every trap locality-year from the hourly CaSR temperatures, `data/processed/pupal_phenology.rds` |
| `R/functions/phenology.R` | Shared helper: percentiles of a seasonal curve |
| `daily_dynamics/` | **Daily dynamics** — the night-by-night distribution of catches: model-versus-observed comparison, diagnostics, the nightly negative binomial model and the `predict_flight()` predictor. Scripts `R/03`-`R/14`, its analysis log and the manuscript drafts |
| `seasonal_dynamics/` | **Seasonal dynamics** — each locality-year as one observation: the timing, spread and shape of the season as functions of temperature. Its own scripts (`R/01`, ...) and analysis log; uses the shared development-model outputs |
| `data/` | Raw, interim and processed data; not version-controlled (see `.gitignore`) |

Scripts are run from the project root, e.g.
`source("daily_dynamics/R/10_nightly_count_model.R")` or
`source("seasonal_dynamics/R/01_season_percentiles.R")`. The analysis logs are
rendered with `quarto render <analysis>/analysis/analysis_log.qmd`.

## Running order

1. Data layer: `R/01_read_trap_data.R`, then `R/02_read_weather_data.R`, then
   `R/03_biosim_phenology.R` and `R/04_pupal_phenology.R`.
2. Daily dynamics: `daily_dynamics/R/05`; then `06`-`10`;
   then `11` and `12` (both read the output of `10`); then `13`, then `14`.
3. Seasonal dynamics: `seasonal_dynamics/R/01_season_percentiles.R`, then
   later scripts in numeric order.

Any change to the data layer requires the downstream scripts of both analyses
to be re-run.

## Data

The trap data come from the Eastern Spruce Budworm Phenology database
(`jncandau/SBWFieldPhenoDatabase`, private until its data paper is published)
and are cloned into `data/raw/`. Only Maine records are used.
