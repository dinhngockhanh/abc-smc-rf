#   Approximate Bayesian Computation sequential Monte Carlo via random forests

##  Installation

The [ABC-SMC-(D)RF](https://doi.org/10.1007/s11222-025-10748-x) library can be installed with

```r
devtools::install_github("dinhngockhanh/abcsmcrf")
```

##  Default perturbation kernel

By default, `smcrf()` uses an independent [Beaumont](https://doi.org/10.1111/j.1541-0420.2008.01180.x) kernel: each parameter is perturbed with a (truncated) normal whose variance is twice the empirical variance of particles resampled from the previous iteration. Omit `rperturb` and `dperturb` to use this kernel, and pass `parameter_bounds` to truncate parameters to their prior support. Parameters not listed in `parameter_bounds` are treated as unbounded.

```r
smcrf(
    method = "smcrf-multi-param",
    statistics_target = statistics_target,
    model = model,
    rprior = rprior,
    dprior = dprior,
    parameter_bounds = data.frame(
        parameter = c("theta1", "theta2"),
        min = c(-Inf, 0),
        max = c(Inf, Inf)
    ),
    nParticles = rep(1000, 10)
)
```

To use a custom kernel, supply both `rperturb` and `dperturb`. See `?make_beaumont_kernel` and `?smcrf`.

##  Vignettes

The `vignettes` folder contains examples of using ABC-SMC-(D)RF, [ABC-RF](https://doi.org/10.1093/bioinformatics/bty867), [ABC-DRF](https://www.jmlr.org/papers/v23/21-0585.html) and traditional ABC methods to infer parameters in different mathematical models, showcased in the paper [Approximate Bayesian computation sequential Monte Carlo via random forests](https://doi.org/10.1007/s11222-025-10748-x).

##  References
1.  Dinh KN, Liu C, Xiang Z, Liu Z, Tavaré S. [Approximate Bayesian computation sequential Monte Carlo via random forests](https://doi.org/10.1007/s11222-025-10748-x). Statistics and Computing 35, 219 (2025).
2.  Raynal L, Marin JM, Pudlo P, Ribatet M, Robert CP, Estoup A. [ABC random forests for Bayesian parameter inference](https://doi.org/10.1093/bioinformatics/bty867). Bioinformatics 35(10):1720-8 (2019).
3.  Ćevid D, Michel L, Näf J, Bühlmann P, Meinshausen N. [Distributional random forests: heterogeneity adjustment and multivariate distributional regression](https://www.jmlr.org/papers/v23/21-0585.html). Journal of Machine Learning Research 23(333):1-79 (2022).
