# Materials and methods

## Light-trap records

Nightly captures of spruce budworm (*Choristoneura fumiferana*) adults came
from the Eastern Spruce Budworm Phenology database (Fournier, Candau et al.,
CC-BY), whose *Flight* sheet holds 763 trapping records in wide form: one row
per site × year × sampling method × sex, with one column per day of year from
160 to 231. Two cell codes, `-` and `m`, both denote a missing count; `-`
additionally pads the out-of-season range, so leading and trailing `-` cells
were dropped and only the trapping period retained, within which both codes
became `NA`.

Records were restricted to Maine light traps operated between 1968 and 1989 —
the overlap with the reanalysis described below — with a season total above 99
moths, a threshold that removes records too sparse to describe a flight curve.
Where one trap was reported as separate male and female records, the two were
summed night by night before the threshold was applied (two traps in 1980).
The working dataset comprised **227 trapping records at 49 localities over 22
years** (43.60–47.44° N, 67.33–70.85° W), reshaped to one row per trap-night:
5,612 counted nights and 154 missing.

Records were flagged rather than filtered, so that every subsequent table
could report the population it referred to. `usable_phenology` (n = 187) marks
records whose daily series is complete enough to describe a flight curve;
`flight_complete` (n = 207) marks those whose trapping period brackets the
modelled flight season; their intersection, `comparable` (n = **170**
locality-years), is the set used for all model comparisons.

## Observed flight phenology

For each locality-year the cumulative catch was formed over the trapped nights
and the 5th, 25th, 50th, 75th and 95th percentiles of catch date obtained by
linear interpolation of that cumulative curve, anchored at zero on the day
before the first trapped night, with the season duration defined as
p95 − p05. The same routine was applied to the simulated curves
(below), so that observed and modelled phenology are summarised identically;
it is implemented once in `R/functions/phenology.R` and shared by every
script.

## Weather

Hourly weather for every trap locality-year was extracted from the
Environment and Climate Change Canada Canadian Surface Reanalysis (CaSR
v3.2), on a ~10 km rotated-pole grid: five variables (precipitation, air
temperature at 1.5 m, downward solar flux, surface pressure and 10 m wind
modulus) over two tiles and six four-year periods spanning 1968–1991. Each
locality was assigned to the grid cell containing it. Air temperature is
stored 3-hourly and was interpolated to hourly with
`chillR::interpolate_gaps_hourly()`, which fits a daily temperature cycle from
the daily minima and maxima and interpolates the residuals; the raw 3-hourly
values were retained alongside so the interpolation remains auditable. The result is one
nested hourly series per locality-year — 227 series, 1,989,744 hourly rows —
timestamped in UTC and in local standard time (UTC−5, no daylight saving),
the latter used for all nightly aggregation.

## Developmental phenology models

Two independent developmental models provided the timing input.

**BioSIM.** The *Spruce_Budworm_Biology* model (Régnière) was run through the
BioSIM Web API at every trap locality-year. BioSIM generates weather for the
requested coordinates from its own historical station network and returns the
daily percentage of a fixed second-instar cohort in each stage, quantised in
steps of 0.25%; the simulation was therefore repeated 10 times per
locality-year with different seeds and the runs averaged to smooth that
quantisation. Natural mortality was switched off (`ApplyMortality = 0`,
`ApplyAdultMortality = 0`) because the object is phenology rather than
survival: under the defaults only a few percent of the cohort reaches the
adult stage and the flight curve is estimated from a handful of individuals.
Male and female series were summed. Of the stage outputs, emergence and
flight are daily fluxes whose cumulative curves are true distributions (male +
female emergence sums to 100% of the cohort), whereas pupae and adults are
stocks and were used only for comparison.

**bayessbw.** Pupation dates were predicted with the Bayesian development
model of Studens et al. (`kdis19/bayessbw`), distributed as the R package
**sbwFieldPheno**. The model fits Sharpe–DeMichele/Schoolfield development-rate
curves to laboratory rearing data for the five larval instars L2–L6 in a
hierarchical form, each posterior draw carrying a rate curve per stage plus a
lognormal individual-variation term. `dev_days()` accumulates development hour
by hour through the CaSR temperature series; pupation is the end of L6. The
New Brunswick colony was used, with 100 posterior draws × 100 simulated
individuals per locality-year, so that within-population variation and
posterior parameter uncertainty are propagated separately. Development was
accumulated from 1 March to 31 August, the window used in bayessbw's own
weather simulations, the first row of the series acting as the date of
post-diapause L2 emergence.

**Regional versus site-specific driving.** Both models were run twice: once at
each trap's own coordinates and once at a single regional point, the centroid
of the trap network (45.718° N, 69.164° W). The two were compared by their
ability to predict the observed median catch date out of sample, and the
regional run was carried forward (see Results).

## Nightly flight-hour temperature

Each trapping night was summarised by the mean air temperature over the hours
20:00–02:00 local standard time, with post-midnight hours assigned to the
preceding evening. A night required at least five of these seven hours to be
present; all nights in the analysis set met that condition. Wind and
precipitation were extracted but not used in the models reported here.

## Model of nightly catches

The catch of night *j* at locality-year *i* was modelled conditionally on the
season total, so that the model describes the *shape* of the flight season and
never its magnitude. Expected shares were formed as

  s_ij = exp(η_ij) / Σ_j exp(η_ij),  E[y_ij] = s_ij · Σ_j y_ij,

with the linear predictor

  η_ij = β_d log(π_ij) + f(d_ij) + γ₁ T_ij + γ₂ T_ij² + u_i,

where π_ij is the driver's simulated percentage on that date (floored at
10⁻³), d_ij is the day of year minus the driver's median date, f(·) a natural
cubic spline with 4 degrees of freedom, T_ij the centred flight-hour
temperature, and u_i a locality-year effect. Models were fitted by maximum
likelihood in **glmmTMB** with a type-2 negative binomial response,
`nbinom2()`, and `(1 | site_year)`.

Seven specifications were compared, crossing three developmental drivers
(BioSIM flight, BioSIM emergence, bayessbw pupation) with the presence or
absence of the temperature terms, plus a temperature-only specification in
which the driver spline was replaced by a 5-degree-of-freedom spline in day of
year. Earlier specifications reported alongside these — a seasonal climatology
fitted from other years, the raw simulated curve without reshaping, and
site-specific rather than regional drivers — were fitted in the same framework.

## Conditional distribution of a night's catch

The distributional family was selected on the 4,319 nights of the analysis set
by AIC among Poisson, type-1 and type-2 negative binomial, Poisson-lognormal
and zero-inflated negative binomial fits with the same mean structure. The
empirical variance–mean relationship was estimated by binning nights on a
fitted mean conditioned only on the locality-year effect, and the candidate
families were further compared by simulating replicate seasons and comparing
their maxima and zero frequencies with the observed ones. Fits were
cross-checked in **lme4** and **MASS**.

## Out-of-sample evaluation

All model comparison was leave-one-year-out: each of the 19 years was
predicted from a model fitted to the other 18, using the training fixed
effects only. The locality-year effect of a held-out year is never estimated
and never needed, because shares are renormalised within each locality-year
and any locality-year-constant term cancels; the season total is an input, so
the only quantity predicted is the distribution of catch over nights. Spline
basis knots were placed once on the pooled data, the design being a fixed
transformation of the driver's own median offset. Three criteria were reported: the absolute error of
the predicted median catch date, the logarithmic score per moth (Σ y log s / Σ
y, larger is better), and the coverage of nominal 90% prediction intervals.
The width of the predicted season was assessed by comparing observed p95 − p05
with the same statistic computed from 100 simulated negative binomial
realisations of each predicted season, rather than from the mean curve, since
a mean curve is necessarily broader than any realisation of it.

## Sharpening

The fitted shape is a log-share curve, so raising it to a power κ > 1 and
renormalising sharpens it around its own peak without displacing the peak. κ
was searched over 1.00–2.50 in steps of 0.05 and selected out of sample by
stacking: for each year, the κ applied to it is the value that scored best on
all other years, so that no year contributes to the κ used on itself. Two
selection criteria were carried in parallel — maximising the log score and
matching the predicted season width — because they do not agree, and the
resulting calibration is reported under both.

## Predictive implementation

The final models were refitted on all years and stored with everything needed
to predict an untrapped locality-year: coefficients, the dispersion parameter,
the spline basis, the temperature centring constant and the calibrated κ.
`predict_flight()` (`R/functions/predict_flight.R`) takes the three quantities
that will be available at prediction time — a developmental distribution
(BioSIM flight or emergence, or bayessbw pupation), an hourly temperature
series, and the season total — and returns, for each night, the expected share
of the season's catch, the expected count, and a negative binomial prediction
interval. Shares are normalised over the nights supplied, so the output is
always conditional on the window requested.

## Diagnostic analyses

Several candidate explanations for the discrepancy between modelled and
observed flight timing were tested and are reported as such.

- **Latitudinal gradient.** Days of flight-date advance per degree of latitude
  were estimated for observed and modelled curves by linear regression with
  year fixed effects, and the variance of the model–observation offset was
  partitioned between and within localities and years.
- **A temperature-dependent pupal stage.** A constant lag from predicted
  pupation to observed median flight was compared, leave-one-year-out, with
  degree-day accumulations over that interval, and with a spline in mean
  interval temperature.
- **Elevation.** Locality elevations were obtained from the USGS 3DEP
  point-query service and compared with latitude; BioSIM was then re-run with
  the measured elevation in place of its own interpolated value for a sample
  of 20 locality-years spread evenly over the elevation range, using the Web
  API directly.
- **Immigration.** Nightly catch series were correlated pairwise within years
  across localities at lags of −7 to +7 days, and the resulting distribution
  compared with a null obtained by pairing localities across shuffled years;
  correlations were also stratified by the north–south position of each pair
  and by inter-trap distance.
- **Catch size.** The season total was added to the gradient regressions as a
  second covariate on the log scale, and the spread of the model–observation
  offset compared across catch-size classes, both between and within
  localities.
- **Coordinate precision.** Trap coordinates were perturbed over 0.4–8 km and
  rounded to 0.01°, and the phenology re-simulated, to bound the sensitivity of
  BioSIM's weather generation to trap position.

## Software

Analyses were carried out in R (4.5.3) with **dplyr**, **tidyr**,
**lubridate**, **splines**, **glmmTMB**, **lme4**, **MASS**, **terra**,
**ncdf4**, **chillR**, **sbwFieldPheno** and **ggplot2**. Every step is a
numbered script in `R/`, writing intermediate objects to `data/processed/`,
and is documented in a Quarto analysis log that recomputes every quoted
figure from those objects at render time. Code and log are available at
`https://github.com/jncandau/Maine_Oct_2026`; the source data are in the
Eastern Spruce Budworm Phenology database.
