# Results

## The observed flight season

The working dataset comprised 226 light-trap records at 49 Maine localities
between 1968 and 1989, 5,590 counted trap-nights with 152 nights missing
inside the trapping periods. Of these, 169 locality-years were both complete
enough to describe a flight curve and long enough to bracket the modelled
flight season, and all model comparison is on that set.

The missing nights were not randomly placed. They occur as blocks of dates
repeated verbatim across localities hundreds of kilometres apart — a
survey-wide servicing schedule rather than individual lost nights — and were
read as censoring: a gap before 15 July bounds the rise of the season, one
after it bounds the decline. Of the 152 missing nights retained in the
dataset, 122 are left-censored and 30 right-censored. Their cost is extremely
skewed: the median affected record had 0.2% of its expected season inside
censored nights and only 3 of 14 assessable records exceeded 1%, the largest
reaching 2.9%. One further record — a trap whose only two censored nights, 8
and 9 July 1978, fell at the peak and carried an expected 33% of the season,
4,395 moths against 9,060 recorded — was removed from the dataset, since
completing those nights moved its median catch date 1.1 days and its p95 1.5
days and no percentile taken from it is identifiable. Completing every
remaining affected record left both headline spatial quantities unchanged —
median offset −1.80 days, observed gradient +0.492 against +0.493 days per
degree north — so censoring constrains which individual records support a
width or tail analysis rather than correcting the results below.

Observed median catch dates ranged from day of year 177 to 206, with a median
of 194.3. Flight was concentrated: the central 90% of a season's catch
occupied a median of 10.6 days, and a median 34% of a season's catch fell on
the single peak night. Zero-catch nights inside the
trapping period accounted for 23.2% of all nights.

## The two developmental models agree with each other

BioSIM and bayessbw were run independently, on different weather (BioSIM's own
generated series against the CaSR reanalysis) and from different biology (a
full seasonal model against a hierarchical development-rate model of the
larval instars). Their predicted pupation dates nevertheless agreed closely:
BioSIM was later by a mean of 0.7 days, with a correlation of 0.80 across the
169 locality-years. In the bayessbw simulations, within-population variation
exceeded posterior parameter uncertainty about fourfold — a population took a
median 9.0 days to pupate from its 5th to its 95th percentile, while the
posterior placed that population's median date to within 2.0 days.

Both models also sat the same distance from the catches. The constant lag
from predicted pupation to observed median catch that minimised the residual
error was 13.0 days, and the two models required lags within a day of each
other.

## Modelled flight is early, and the error is in the spatial gradient

Against the catches, modelled flight was early by a median of 1.8 days at the
median of the season, 3.8 days at p05 and 0.9 days at p95, so the modelled
season was also 3.5 days too long. The average conceals a split: year means of
modelled and observed median date correlated at 0.84, but within a year,
across localities, the correlation was only 0.14. The models place the season
correctly in time and incorrectly in space.

The cause is a gradient-slope error. With year fixed effects removed, the
observed median catch date advanced by only +0.49 days per degree of latitude,
a slope not distinguishable from zero, while both models predicted four to
five times that (Table 1). Across the 3.57° of latitude spanned by the
network, the observed flight date spread 1.8 days; BioSIM spread it 15.0 days
and bayessbw 11.7 days. The two models' offsets correlated at 0.74, so this is
one error shared by both, not two independent ones.

**Table 1. Days per degree of latitude, with year fixed effects.**

| Series | n | d per °N | 95% CI | implied spread, d |
|---|---|---|---|---|
| Observed median catch | 169 | +0.49 | −0.14 to +1.13 | 1.8 |
| BioSIM modelled flight | 169 | +4.19 | +3.37 to +5.01 | 15.0 |
| bayessbw pupation | 169 | +3.29 | +2.55 to +4.02 | 11.7 |

## Five candidate explanations were tested and set aside

**A temperature-dependent pupal stage.** Replacing the constant
pupation-to-flight lag with degree-day accumulation improved leave-one-year-out
RMSE from 6.58 to 6.49 days — 0.09 days. The interval from predicted pupation
to observed flight was uncorrelated with mean temperature over that interval
(r = −0.04), and a spline in interval temperature explained 5.5% of its
variance. Most of the variance in the interval was spatial rather than thermal:
63.9% lay between localities, 12.4% between years. At site level the sign was
backwards — the interval grew longer at warmer July temperatures (r = +0.47)
and shorter at higher latitude (r = −0.40). No development function of
temperature can absorb an error structured this way.

**Elevation.** The 49 localities ranged from 7 to 573 m (median 218 m) and
elevation was itself a latitude proxy in this network (r = +0.56; southern
localities averaged 68 m against 303 m in the north), so it was a plausible
confound for the whole gradient. Supplying measured USGS 3DEP elevations to
BioSIM in place of its own interpolated values shifted predicted flight by
+0.27 days on average and reduced the latitudinal gradient only from +6.12 to
+5.98 days per degree. A local sensitivity scan gave −1.84 d per 100 m at a
coastal locality and +1.46 d per 100 m inland, opposite in sign, which
suggests elevation partly acts on BioSIM through station selection rather than
through lapse rate alone; the mechanism was not pursued.

**A detection threshold at low catch.** If small catches simply clipped the
tails of the flight curve, they would appear narrower and later-starting. The
opposite held: the share of a season's catch falling after the modelled window
correlated −0.31 with log catch size, and season width correlated −0.23, so
small catches carried *more* late catch and slightly wider curves.

**Immigration at nightly scale.** Across 875 within-year locality pairs the
median correlation between nightly catch series was 0.134, against −0.088 for
a null in which localities were paired across shuffled years; 20.2% of real
pairs exceeded r = 0.5, against 6.5% of null pairs. Correlation decayed with
inter-trap distance (r = −0.17) and was weakest for north–south pairs
(r = −0.18), and the best lag was zero. Localities therefore share weather at
nightly scale, as expected, but there is no signature of synchronous arrivals
that would explain a systematic southern excess.

**Catch size.** Season total survived as a second covariate on the offset,
nearly independent of latitude (r = 0.14 between the two). Fitted jointly, the
offset increased by +3.33 days per degree north and +1.87 days per tenfold
increase in catch, with 0.210 of the variance unique to latitude and 0.043
unique to catch. The effect is between localities, not within them: mean
offset correlated 0.40 with mean catch across localities, but only 0.07 with
catch within a locality across years. A locality's poor years are not its
low-catch years. The spread of the offset was also wider at small catches (sd
7.4 days below 1,000 moths against 4.6 days above 5,000), so part of the
catch-size effect is estimation precision rather than biology — but the shift
in the mean is not.

## Regional weather outperforms site-specific weather

Because the error is spatial, the site-specific weather driving each
simulation was tested directly. Predicting the observed median catch date with
a leave-one-year-out intercept, BioSIM run at each trap's own coordinates gave
an RMSE of 6.50 days, while a single BioSIM run at the network centroid gave
4.30 days (Table 2). The same held for bayessbw, 6.50 against 5.15 days. The
floor for this comparison — predicting a locality-year from the mean observed
date of the other localities that year — was 3.64 days, so one regional run
came within 0.66 days of knowing the answer from the other traps.

Site-specific weather thus makes both models worse by 1.4–2.2 days of RMSE.
Whatever spatial signal the weather carries into the models is, in this
network, noise with respect to the catches.

**Table 2. Predicting the observed median catch date (169 locality-years).**

| Predictor | varies by site | RMSE, d | MAE, d | r |
|---|---|---|---|---|
| BioSIM, site-specific | yes | 6.50 | 4.95 | 0.49 |
| BioSIM, one regional point | no | 4.30 | 3.39 | 0.68 |
| bayessbw, site-specific | yes | 6.50 | 4.95 | 0.40 |
| bayessbw, one regional run | no | 5.15 | 4.11 | 0.52 |
| Other localities that year | yes | 3.64 | 2.60 | 0.77 |

## A nightly model of the shape of the season

Nine specifications of the nightly count model were compared over 4,297 nights
(Table 3). The regional driver with an empirically reshaped curve (`reg+`)
gave a mean absolute error of 3.03 days on the median catch date, with 83% of
locality-years within 5 days, against 3.76 days for the same model on
site-specific curves and 4.46 days for the raw simulated curve or a seasonal
climatology fitted from other years. Adding latitude as a covariate did not
help (3.09 days).

The empirical correction is large. The exponent on the regional curve was
0.240, so only about a quarter of the simulated curve's shape survived
fitting; the correction moved the predicted median by +1.91 days on average
(sd 4.9) and widened the predicted central 90% from 12.6 to 15.8 days, against
11.0 observed. Replacing the whole regional curve by its median date alone
(`reg-med`) cost nothing — 3.00 days — which locates the models' contribution
precisely: they supply the date, and the empirical spline supplies the shape.

**Table 3. Leave-one-year-out skill of the nightly count specifications.**

| Specification | mean \|median error\|, d | within 5 d |
|---|---|---|
| Regional curve, median only (`reg-med`) | 3.00 | 82% |
| Regional curve, reshaped (`reg+`) | 3.03 | 83% |
| Regional + latitude | 3.09 | 80% |
| Site-specific curve, reshaped | 3.76 | 70% |
| Site-specific + latitude | 3.99 | 66% |
| bayessbw pupation, site-specific | 4.01 | 69% |
| bayessbw pupation, regional | 4.21 | 64% |
| Seasonal climatology | 4.46 | 65% |
| Raw simulated curve | 4.46 | 63% |

## The conditional distribution of a night's catch

Over the same nights, the Poisson model is untenable by orders of magnitude
(ΔAIC > 2.2 × 10⁶). The Poisson-lognormal fit had the lowest AIC (38,111),
with the type-2 negative binomial 54 AIC behind (38,165), the zero-inflated
negative binomial a further 2 behind and the type-1 negative binomial 1,682
behind (Table 4). Zero inflation is therefore unnecessary: the negative
binomial already reproduces the 998 observed zero nights (1,150 simulated),
and adding a zero component changes nothing but the parameter count.

The empirical variance–mean slope was 2.07, close to the quadratic of NB2 and
away from the linear form of NB1. The fitted NB2 dispersion was θ = 0.290,
i.e. Var = μ + 3.45 μ², a multiplicative noise with a coefficient of variation
of 1.86: a night expected to yield 100 moths has a standard deviation of 186
and an 18% chance of catching none. The two leading families separate on the
extremes rather than the centre: simulating from the Poisson-lognormal fit
produced a maximum night of 221,600 moths against 27,500 observed, while NB2
produced 46,137. NB2 was therefore adopted for prediction and simulation, and
the Poisson-lognormal retained as the interpretable description; neither
should be trusted for extreme quantiles.

**Table 4. Distribution of a night's catch, same mean structure (4,297 nights).**

| Family | parameters | AIC | ΔAIC | zeros simulated | max simulated |
|---|---|---|---|---|---|
| Poisson-lognormal | 8 | 38,111 | 0 | 943 | 221,600 |
| Negative binomial (NB2) | 8 | 38,165 | 54 | 1,150 | 46,137 |
| Zero-inflated NB2 | 9 | 38,167 | 56 | 1,160 | 53,778 |
| Negative binomial (NB1) | 8 | 39,793 | 1,682 | 1,320 | 6,377 |
| Poisson | 7 | 2,234,604 | 2,196,494 | 164 | 14,169 |

Observed: 998 zero nights, maximum 27,500 moths.

## Which forecast input is needed

The three developmental distributions that will be available at prediction
time were then scored against each other in the same framework (Table 5).
BioSIM's flight and emergence curves were interchangeable for timing — 3.05
and 3.19 days of median absolute error — while bayessbw pupation cost about a
day (4.17 days) and a temperature-only model with no developmental input at
all gave 4.18 days.

Nightly flight-hour temperature was the first covariate to improve the nightly
numbers without improving the median: adding it moved the log score per moth
from −2.912 to −2.839 with the flight curve and from −3.016 to −2.885 with
emergence, while the median error stayed within 0.1 days. Nominal 90%
prediction intervals covered 93% of nights under every specification.

**Table 5. Forecast inputs, leave-one-year-out (169 locality-years).**

| Input | median \|error\|, d | width error, d | log score per moth | 90% coverage |
|---|---|---|---|---|
| BioSIM flight | 3.05 | +3.44 | −2.912 | 0.93 |
| BioSIM flight + temperature | 3.15 | +3.19 | −2.839 | 0.93 |
| BioSIM emergence + temperature | 3.17 | +3.85 | −2.885 | 0.93 |
| BioSIM emergence | 3.19 | +4.07 | −3.016 | 0.93 |
| bayessbw pupation + temperature | 3.82 | +5.39 | −2.960 | 0.93 |
| bayessbw pupation | 4.17 | +5.55 | −3.160 | 0.94 |
| Temperature only | 4.18 | +6.70 | −3.141 | 0.93 |

## The predicted season is too broad, and can be sharpened

Every specification predicted a season broader than observed, and the excess
is real but smaller than a naive comparison suggests. Against an observed
central 90% of 11.0 days, the mean predicted curve spanned 14.2 days under the
best specification, but a simulated realisation of that curve spanned 13.4
days, since any single realisation is narrower than the curve it comes from.
Only 73% of observed widths fell inside the simulated 90% band (69% for the
flight curve without temperature, 60% for pupation, 51% for temperature only),
so the over-width is of the order of 2.5 days rather than the 5 days a
curve-to-realisation comparison would imply.

Sharpening the fitted log-share curve by a power κ, chosen out of sample by
stacking, removed almost all of it — but only under one of two criteria
(Table 6). By log score, κ = 1 for every driver: the proper scoring rule
prefers the broad curve, because sharpening moves mass off the shoulders where
a few moths still arrive. By width matching, κ ranged from 1.34 to 1.73 and
reduced the width error from +3.19 to +0.19 days for the flight curve, at a
cost of 0.09 days on the median and 0.064 on the log score per moth. Both
calibrations were retained, since the choice depends on whether the quantity
of interest is the window or a particular night's count.

**Table 6. Sharpening, chosen out of sample by stacking.**

| Driver | criterion | κ | median error, d | width error, d | log score per moth |
|---|---|---|---|---|---|
| BioSIM flight | none | 1.00 | 3.15 | +3.19 | −2.839 |
| BioSIM flight | log score | 1.00 | 3.15 | +3.19 | −2.839 |
| BioSIM flight | width | 1.34 | 3.24 | +0.19 | −2.903 |
| BioSIM emergence | none | 1.00 | 3.17 | +3.85 | −2.885 |
| BioSIM emergence | log score | 1.00 | 3.17 | +3.85 | −2.885 |
| BioSIM emergence | width | 1.46 | 3.28 | +0.24 | −2.984 |
| bayessbw pupation | none | 1.00 | 3.82 | +5.39 | −2.960 |
| bayessbw pupation | log score | 1.00 | 3.82 | +5.39 | −2.960 |
| bayessbw pupation | width | 1.73 | 4.10 | +0.54 | −3.219 |

## Worked prediction

Applied to two locality-years from their emergence curve, hourly temperature
and season total alone, the calibrated predictor placed the median catch date
2.3 days early at East Allagash Township in 1983 — the largest record in the
dataset, 71,730 moths over 32 nights — and 2.6 days late at Meddybemps
Township in 1981 (10,466 moths). The temperature term is visible in the nightly trace: at
East Allagash, day 187 averaged 6.9 °C through the flight hours and caught a
single moth, against hundreds on the adjacent warm nights.
