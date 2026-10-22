# grmRT

**Joint graded response and distance-difficulty response-time models for Likert-type questionnaires**

`grmRT` fits a joint model for Likert-type item responses and item response
times. Responses follow Samejima's graded response model (GRM). Log response
times follow a *distance-difficulty* (inverted-U) model: respondents answer
more slowly when their trait level lies close to the item location, and more
quickly when the item is clearly (dis)agreeable for them
(Kuncel, 1973; Ferrando & Lorenzo-Seva, 2007; Molenaar et al., 2021).

```
P(X_ij >= k | theta_i) = logistic(a_j (theta_i - b_jk))
log T_ij = lambda_j - tau_i - gamma_j h(theta_i - delta_j) + e_ij
```

with `h(u) = sqrt(u^2 + kappa) - sqrt(kappa)` a smooth version of `|u|`.

## What the package adds

* **Fast frequentist estimation.** Marginal maximum likelihood with an EM
  algorithm (trait grid x adaptive Gauss-Hermite quadrature for speed),
  suitable for tens of thousands of respondents. Standard errors use the
  cross-product of person-level scores.
* **Robust timing model.** Optional Student-*t* errors with item-specific
  degrees of freedom for noisy timings, for example inter-click intervals in
  single-page web questionnaires.
* **Item diagnostics.** Wald and likelihood-ratio tests of `gamma_j = 0` per
  item, with multiplicity adjustment, and tests of asymmetry. Items whose
  times do not slow down near the item location are flagged for review.
* **Collateral information.** Trait scores (EAP) that use response times,
  with the reduction in posterior SD relative to responses-only scores.
* **Asymmetric slopes** below and above the item location, identified up to
  a common shift and reported relative to the item average.

## Installation

```r
# install.packages("remotes")
remotes::install_github("agustinus-bambang/grmRT")
```

## Example

```r
library(grmRT)
set.seed(1)
sim <- simulate_grmrt(n_persons = 1000, n_items = 8,
                      gamma = c(rep(0.2, 7), 0))      # item 8 has no effect
fit <- grmrt_fit(sim$responses, sim$logrt)            # normal errors
fit_t <- grmrt_fit(sim$responses, sim$logrt, error = "t", start = fit)
anova(fit, fit_t)                                      # normal vs t errors
summary(fit_t)
item_diagnostics(fit_t)                                # flags item 8
head(score_persons(fit_t))                             # trait and speed EAPs
plot_time_distance(fit_t)                              # inverted-U check
```

For raw survey data, `prep_rt()` reorders items into their order of
presentation, recodes reverse-keyed items, converts times to log seconds,
and treats implausible timings as missing:

```r
dat <- prep_rt(df, items = item_cols, times = time_cols, order = shown_order,
               time_unit = "ms", reverse = reversed_items,
               min_time = 0.5, max_time = 60, drop_first = TRUE)
fit <- grmrt_fit(dat, error = "t")
```

## Reproducibility

The `inst/replication` folder of the source repository contains the
simulation study and the analysis of the IPIP-FFM data (Open-Source
Psychometrics Project). The data are not redistributed; download
`IPIP-FFM-data-8Nov2018.zip` from <https://openpsychometrics.org/_rawdata/>.

## References

Ferrando, P. J., & Lorenzo-Seva, U. (2007). A measurement model for Likert
responses that incorporates response time. *Multivariate Behavioral Research,
42*(4), 675-706.

Kuncel, R. B. (1973). Response processes and relative location of subject and
item. *Educational and Psychological Measurement, 33*(3), 545-563.

Molenaar, D., Rózsa, S., & Kõ, N. (2021). Modeling asymmetry in the
time-distance relation of ordinal personality items. *Applied Psychological
Measurement, 45*(3), 178-194.
