# Overview
kmsmooth is an empirical method to modify the Kaplan-Meier estimator 
to produce a smooth estimate of the survival functions with the right-censored data:

`ESAP_supp()` computes the support of the survival time.

`ESAP_pdf()` computes the density at a given failure time.

`ESAP_survival()`computes the survival probabiltiy at a given failure time.

`ESAP_median()` computes the median of the survival time

`ESAP_t0()` computes the probability of survival at time zero.

`ESAP_plotSurvival()` plots the survival function over the support computed by `ESAP_Supp()`.

# Installation
```{r }
# Install the the development version from GitHub:
devtools::install_github("PratheepaJega/kmsmooth")
```

If you find a bug, please report it with a reproducible example on [GitHub](https://github.com/PJega31/kmsmooth/issues). 

# Usage
```{r }
library(kmsmooth)
# Simulate a data

# failure time
time=

# status
status=

# Find the support
supp=ESAP_supp(time,status)

#Create a grid to compute pdf and CDF
t.grid=seq(supp[1],supp[2],by=.1)

# Compute the pdf
pdft=ESAP_pdf(time, status,t.grid)

# Compute the survival
sur=ESAP_survival(time, status,t.grid)

# Compute the probability of survival or failure at time zero
prob.tzero=ESAP_t0(time, status,t=0)

# Compute the median survival/failure time
median.sur=ESAP_median(time,status,supp)

# Plot the survival function over the support
survival.gg=ESAP_plotSurvival(time,status,supp=supp)
```
