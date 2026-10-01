# Results

## The observed flight season

The working dataset comprised 226 light-trap records at 49 Maine localities
between 1968 and 1989. Reading a `-` cell inside a year's flight season as a
zero rather than as a missing value gives 8,286 counted trap-nights, 2,720 of
them such zeros, with only 30 nights left missing. Of these records, 203
locality-years were both complete enough to describe a flight curve and long
enough to bracket the modelled flight season, and all model comparison is on
that set. Half of all counted nights caught nothing.

The nights still coded as missing were treated as censoring, a gap before 15
July bounding the rise of the season and one after it the decline. After the
season window is applied only 30 such nights remain, in four records, and all
of them are right-censored: the large blocks of missing June dates fall before
the first capture of their year and leave the dataset with the rest of the
out-of-season range. None of the four carries as much as 1% of its expected
season in those nights. One further record — a trap whose only two censored
nights, 8 and 9 July 1978, fell at its peak and carried an expected 26% of the
season — was removed from the dataset, since no percentile taken from it is
identifiable.

Observed median catch dates ranged from day of year 176 to 206, with a median
of 194.0. Flight was concentrated: the central 90% of a season's catch
occupied a median of 10.9 days, and a median 34% of a season's catch fell on
the single peak night.

## The two developmental models agree on the pattern, not the date

BioSIM and bayessbw were run independently, on different weather (BioSIM's own
generated series against the CaSR reanalysis) and from different biology (a
full seasonal model against a hierarchical development-rate model of the
larval instars). Their predicted pupation dates ranked the 203 locality-years
alike (r = 0.77), but BioSIM pupated earlier, by a mean of 4.3 days (median
4.0). In the bayessbw simulations, within-population variation
exceeded posterior parameter uncertainty about fourfold — a population took a
median 9.0 days to pupate from its 5th to its 95th percentile, while the
posterior placed that population's median date to within 2.0 days.

The constant lag from predicted pupation to observed median catch was
correspondingly different: a median 13.4 days for bayessbw and 17.5 days for
BioSIM. The four-day gap is the models' disagreement on the pupation date
itself; the spread of the interval around its median was the same for both
(sd 6.6 and 6.7 days).

## Modelled flight is early, and the error is in the spatial gradient

Against the catches, modelled flight was early by a median of 1.8 days at the
median of the season, 2.1 days at the median, 3.7 days at p05 and 1.0 days at
p95, so the modelled season was also 3.3 days too long. The average conceals a
split: year means of modelled and observed median date correlated at 0.84, but
within a year, across localities, the correlation was only 0.20. The models
place the season correctly in time and incorrectly in space.

The cause is a gradient-slope error. With year fixed effects removed, the
observed median catch date advanced by +0.66 days per degree of latitude,
against four to six times that in both models (Table 1). Across the 3.8° of
latitude spanned by the network, the observed flight date spread 2.5 days;
BioSIM spread it 16.0 days and bayessbw 11.5 days. The observed slope is small
but, on this larger set of records, no longer indistinguishable from zero. The
two models' offsets correlated at 0.71, so this is one error shared by both,
not two independent ones.

**Table 1. Days per degree of latitude, with year fixed effects.**

| Series | n | d per °N | 95% CI | implied spread, d |
|---|---|---|---|---|
| Observed median catch | 203 | +0.66 | +0.07 to +1.24 | 2.5 |
| BioSIM modelled flight | 203 | +4.17 | +3.43 to +4.90 | 16.0 |
| bayessbw pupation | 203 | +2.99 | +2.30 to +3.69 | 11.5 |

## Five candidate explanations were tested and set aside

**A temperature-dependent pupal stage.** Replacing the constant
pupation-to-flight lag with degree-day accumulation did not help at all:
leave-one-year-out RMSE was 6.75 days against 6.76 days. The interval from
predicted pupation to observed flight was uncorrelated with mean temperature
over that interval (r = −0.01), and a spline in interval temperature explained
4.7% of its variance. Most of the variance in the interval was spatial rather
than thermal: 64.5% lay between localities, 11.5% between years. At site level the sign was
backwards — the interval grew longer at warmer July temperatures (r = +0.47)
and shorter at higher latitude (r = −0.40). No development function of
temperature can absorb an error structured this way.

**Elevation.** The 49 localities ranged from 7 to 573 m (median 218 m) and
elevation was itself a latitude proxy in this network (r = +0.56; southern
localities averaged 68 m against 303 m in the north), so it was a plausible
confound for the whole gradient. Supplying measured USGS 3DEP elevations to
BioSIM in place of its own interpolated values shifted predicted flight by
−0.08 days on average and reduced the latitudinal gradient only from +4.43 to
+4.33 days per degree in the sampled locality-years. A local sensitivity scan gave −1.84 d per 100 m at a
coastal locality and +1.46 d per 100 m inland, opposite in sign, which
suggests elevation partly acts on BioSIM through station selection rather than
through lapse rate alone; the mechanism was not pursued.

**A detection threshold at low catch.** If small catches simply clipped the
tails of the flight curve, they would appear narrower and later-starting. The
opposite held: the share of a season's catch falling after the modelled window
correlated −0.30 with log catch size, and season width correlated −0.23, so
small catches carried *more* late catch and slightly wider curves.

**Immigration at nightly scale.** Across 1,282 within-year locality pairs the
median correlation between nightly catch series was 0.195, against −0.039 for
a null in which localities were paired across shuffled years; 21.0% of real
pairs exceeded r = 0.5, against 6.4% of null pairs. These correlations are
higher than before the zero-fill, since localities now agree on shared empty
nights as well as on catches. Correlation decayed with
inter-trap distance (r = −0.17) and was weakest for north–south pairs
(r = −0.18), and the best lag was zero. Localities therefore share weather at
nightly scale, as expected, but there is no signature of synchronous arrivals
that would explain a systematic southern excess.

**Catch size.** Season total survived as a second covariate on the offset,
nearly independent of latitude (r = 0.22 between the two). Fitted jointly, the
offset increased by +3.35 days per degree north and +1.43 days per tenfold
increase in catch, with 0.212 of the variance unique to latitude and 0.026
unique to catch. The effect is between localities, not within them: mean
offset correlated 0.38 with mean catch across localities, but only 0.06 with
catch within a locality across years. A locality's poor years are not its
low-catch years. The spread of the offset was also wider at small catches (sd
7.3 days below 1,000 moths against 4.5 days above 5,000), so part of the
catch-size effect is estimation precision rather than biology — but the shift
in the mean is not.

## Regional weather outperforms site-specific weather

Because the error is spatial, the site-specific weather driving each
simulation was tested directly. Predicting the observed median catch date with
a leave-one-year-out intercept, BioSIM run at each trap's own coordinates gave
an RMSE of 6.57 days, while a single BioSIM run at the network centroid gave
4.38 days (Table 2). The same held for bayessbw, 6.73 against 5.13 days. The
floor for this comparison — predicting a locality-year from the mean observed
date of the other localities that year — was 3.80 days, so one regional run
came within 0.58 days of knowing the answer from the other traps.

Site-specific weather thus makes both models worse by 1.6–2.2 days of RMSE.
Whatever spatial signal the weather carries into the models is, in this
network, noise with respect to the catches.

**Table 2. Predicting the observed median catch date (203 locality-years).**

| Predictor | varies by site | RMSE, d | MAE, d | r |
|---|---|---|---|---|
| BioSIM, site-specific | yes | 6.57 | 5.01 | 0.49 |
| BioSIM, one regional point | no | 4.38 | 3.45 | 0.65 |
| bayessbw, site-specific | yes | 6.73 | 5.13 | 0.38 |
| bayessbw, one regional run | no | 5.13 | 4.11 | 0.52 |
| Other localities that year | yes | 3.80 | 2.74 | 0.74 |

## A nightly model of the shape of the season

Nine specifications of the nightly count model were compared over 7,481 nights
(Table 3). The regional driver with an empirically reshaped curve (`reg+`)
gave a mean absolute error of 3.21 days on the median catch date, with 80% of
locality-years within 5 days, against 4.54 days for the same model on
site-specific curves and 5.04–5.09 days for the raw simulated curve or a
seasonal climatology fitted from other years. Adding latitude as a covariate
did not help (3.31 days).

The empirical correction is large. The exponent on the regional curve was
0.240, so only about a quarter of the simulated curve's shape survived
fitting; the correction moved the predicted median by +2.77 days on average
(sd 5.7) and widened the predicted central 90% from 13.9 to 15.8 days, against
11.2 observed. Replacing the whole regional curve by its median date alone
(`reg-med`) cost nothing — 3.21 days — which locates the models' contribution
precisely: they supply the date, and the empirical spline supplies the shape.

**Table 3. Leave-one-year-out skill of the nightly count specifications.**

| Specification | mean \|median error\|, d | within 5 d |
|---|---|---|
| Regional curve, reshaped (`reg+`) | 3.21 | 80% |
| Regional curve, median only (`reg-med`) | 3.21 | 80% |
| Regional + latitude | 3.31 | 79% |
| bayessbw pupation, regional | 4.25 | 63% |
| Site-specific curve, reshaped | 4.54 | 64% |
| Site-specific + latitude | 4.58 | 65% |
| bayessbw pupation, site-specific | 4.85 | 61% |
| Raw simulated curve | 5.04 | 62% |
| Seasonal climatology | 5.09 | 53% |

## The conditional distribution of a night's catch

Over the same nights, the Poisson model is untenable by orders of magnitude
(ΔAIC > 2.5 × 10⁶). The Poisson-lognormal fit had the lowest AIC (46,220),
with the type-2 negative binomial 159 AIC behind (46,379), the zero-inflated
negative binomial a further 2 behind and the type-1 negative binomial 1,308
behind (Table 4). Zero inflation is therefore unnecessary even though half the
nights are zeros: the negative binomial already reproduces the 3,715 observed
zero nights (3,685 simulated), and adding a zero component changes nothing but
the parameter count.

The empirical variance–mean slope was 1.83, between the linear form of NB1 and
the quadratic of NB2 and closer to the latter once the fit is weighted by the
bulk of the data. The fitted NB2 dispersion was θ = 0.190, i.e.
Var = μ + 5.26 μ², a multiplicative noise with a coefficient of variation of
2.29: a night expected to yield 100 moths has a standard deviation of 230 and
a 30% chance of catching none. The two leading families separate on the
extremes rather than the centre, and the zero-filled nights widen that gap:
simulating from the Poisson-lognormal fit produced a maximum night of 801,923
moths against 27,500 observed, while NB2 produced 100,658. NB2 was therefore adopted for prediction and simulation, and
the Poisson-lognormal retained as the interpretable description; neither
should be trusted for extreme quantiles.

**Table 4. Distribution of a night's catch, same mean structure (7,481 nights).**

| Family | parameters | AIC | ΔAIC | zeros simulated | max simulated |
|---|---|---|---|---|---|
| Poisson-lognormal | 8 | 46,220 | 0 | 3,543 | 801,923 |
| Negative binomial (NB2) | 8 | 46,379 | 159 | 3,685 | 100,658 |
| Zero-inflated NB2 | 9 | 46,381 | 161 | 3,687 | 89,439 |
| Negative binomial (NB1) | 8 | 47,528 | 1,308 | 3,965 | 6,673 |
| Poisson | 7 | 2,553,113 | 2,506,893 | 1,099 | 17,000 |

Observed: 3,715 zero nights, maximum 27,500 moths.

## Which forecast input is needed

The three developmental distributions that will be available at prediction
time were then scored against each other in the same framework (Table 5).
BioSIM's flight and emergence curves remained close for timing — 3.23 and 3.44
days of median absolute error — while bayessbw pupation now cost about a day
and a half (4.93 days), no better than a temperature-only model with no
developmental input at all (4.82 days).

Nightly flight-hour temperature again improved the nightly numbers more than
the median: adding it moved the log score per moth from −2.962 to −2.901 with
the flight curve and from −3.061 to −2.945 with emergence, while the median
error moved by less than 0.2 days. Nominal 90% prediction intervals covered
93–94% of nights under every specification.

**Table 5. Forecast inputs, leave-one-year-out (203 locality-years).**

| Input | median \|error\|, d | width error, d | log score per moth | 90% coverage |
|---|---|---|---|---|
| BioSIM flight | 3.23 | +2.78 | −2.962 | 0.93 |
| BioSIM flight + temperature | 3.25 | +2.38 | −2.901 | 0.93 |
| BioSIM emergence + temperature | 3.33 | +3.28 | −2.945 | 0.94 |
| BioSIM emergence | 3.44 | +3.58 | −3.061 | 0.94 |
| bayessbw pupation + temperature | 4.68 | +4.53 | −3.075 | 0.93 |
| Temperature only | 4.82 | +5.42 | −3.248 | 0.93 |
| bayessbw pupation | 4.93 | +4.75 | −3.239 | 0.93 |

## The predicted season is too broad, and can be sharpened

Every specification predicted a season broader than observed, and the excess
is real but smaller than a naive comparison suggests. Against an observed
central 90% of 11.6 days, the mean predicted curve spanned 13.6 days under the
best specification, but a simulated realisation of that curve spanned 13.2
days, since any single realisation is narrower than the curve it comes from.
76% of observed widths fell inside the simulated 90% band (71% for the flight
curve without temperature, 68% for pupation, 66% for temperature only), so the
over-width is of the order of 1.6 days rather than the 2 days a
curve-to-realisation comparison would imply.

Sharpening the fitted log-share curve by a power κ, chosen out of sample by
stacking, removed almost all of it — but only under one of two criteria
(Table 6). By log score, κ = 1 for every driver: the proper scoring rule
prefers the broad curve, because sharpening moves mass off the shoulders where
a few moths still arrive. By width matching, κ ranged from 1.25 to 1.51 and
reduced the width error from +2.33 to +0.42 days for the flight curve, at a
cost of 0.01 days on the median and 0.053 on the log score per moth. With the
emergence curve the sharpened fit was marginally better on the median as
well. Both
calibrations were retained, since the choice depends on whether the quantity
of interest is the window or a particular night's count.

**Table 6. Sharpening, chosen out of sample by stacking.**

| Driver | criterion | κ | median error, d | width error, d | log score per moth |
|---|---|---|---|---|---|
| BioSIM flight | none | 1.00 | 3.25 | +2.33 | −2.901 |
| BioSIM flight | log score | 1.00 | 3.25 | +2.33 | −2.901 |
| BioSIM flight | width | 1.25 | 3.26 | +0.42 | −2.953 |
| BioSIM emergence | none | 1.00 | 3.33 | +3.23 | −2.945 |
| BioSIM emergence | log score | 1.00 | 3.33 | +3.23 | −2.945 |
| BioSIM emergence | width | 1.38 | 3.30 | +0.54 | −3.029 |
| bayessbw pupation | none | 1.00 | 4.68 | +4.48 | −3.075 |
| bayessbw pupation | log score | 1.00 | 4.68 | +4.48 | −3.075 |
| bayessbw pupation | width | 1.51 | 4.73 | +0.78 | −3.284 |

## Worked prediction

Applied to two locality-years from their emergence curve, hourly temperature
and season total alone, the calibrated predictor placed the median catch date
1.9 days early at East Allagash Township in 1983 — the largest record in the
dataset, 71,730 moths over 44 nights — and 3.4 days late at Meddybemps
Township in 1981 (10,466 moths over 37 nights), with 95% and 97% of nights
inside the nominal 90% interval. The temperature term is visible in the nightly trace: at
East Allagash, day 187 averaged 6.9 °C through the flight hours and caught a
single moth, against hundreds on the adjacent warm nights.
