# Overview
ESPA is an empirical saddlepoint approximation method for smoothing Kaplan-Meier estimator with the right-censored data.

ESPA provides method to compute the support for the survival time, compute the probability density for the survival time, and compute the distribution for the survival time, compute the median for the survival time, compute the probability of survival at time zero, and plot the survival function.

`ESAP_supp()` computes the support of the survival time.

`ESAP_pdf()` computes the density at a given failure time.

`ESAP_survival()`computes the survival probabiltiy at a given failure time.

`ESAP_median()` computes the median of the survival time

`ESAP_t0()` computes the probability of survival at time zero.

`ESAP_plotSurvival()` plots the survival function over the support computed by `ESAP_Supp()`.

# Installation
```{r }
# Install the the development version from GitHub:
devtools::install_github("PratheepaJega/ESPA")
```

If you find a bug, please report it with a reproducible example on [GitHub](https://github.com/PratheepaJ/ESPA/issues). 

