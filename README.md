# BDFM
Bayesian dynamic factor model estimation in R.

[![Build Status](https://travis-ci.org/christophsax/tsbox.svg?branch=master)](https://travis-ci.org/srlanalytics/BDFM)
[![codecov](https://codecov.io/github/srlanalytics/BDFM/branch/master/graphs/badge.svg)](https://codecov.io/github/srlanalytics/BDFM)

## Description

BDFM estimates dynamic factor models by simulation using the Durbin and Koopman (2012) disturbance smoother and returns estimated factors, predicted values of observables, posterior distributions for predicted values (optional), and forecast updates by series (optional). Maximum likelihood estimation via Watson and Engle (1983) is also supported. Input data may be noisy, have missing values, or "ragged edges" due to different start or end dates.

## Installation

```t
remotes::install_github("srlanalytics/BDFM")
```

## Basic Use

Estimate a model with two factors and two lags in the transition equation
```r
library(BDFM)
Est <- bdfm(Y = Data, factors = 2, lags = 2)
```
Input data is any data type that can be coverted by ```as.matrix() ``` and should index time in rows and series in columns. If the library [tsbox](https://github.com/christophsax/tsbox) is present time seires attributes of input data are preserved for any `ts_boxable()` data format. In this case data need not be entered in tabular format for types such as `data.table` or `data.frame`.

Quick time series plot estimated factors
```r
ts.plot(Est$factors)
```
Quick time series plot of fitted values (including forecasts and nowcasts) for first series in input data
```r
ts.plot(Est$values[,1])
```

## Full Input

```r
Est <- bdfm(Y, factors = 1, lags = 2, forecast = 0, B_prior = NULL, lam_B = 0, H_prior = NULL, lam_H = 0, nu_q = 0, nu_r = NULL, ID = 'PC_full', intercept = TRUE, reps = 1000, burn = 500)
```

```Y```           Input data

```factors```     Number of factors

```lags```        Number of lags in transition equation

```forecast```    Number of periods ahead to forecast

```B_prior```     Prior on transtion equation (default is to shrink values towards zero)

```lam_B```       Prior tightness on transition equation

```H_prior```    Prior on loadings/observation equation (default is zero)

```lam_H```      Prior tightness on loadings, scalar

 ```nu_q```      Prior "degrees of freedom" for transition equation (i.e. shrink covariance of shocks towards zero), scalar

```nu_r```       Prior "degrees of freedom" for observation equation (i.e. shrink covariance of shocks towards zero), vector with elements equaling number of observed series

```ID```        Model identification technique. 'PC_full' identifies on principal components using all series. 'PC-sub' identifies on principal components using a subset of the data that maximizes the number of observations when many observations are missing.

```intercept``` T/F estimate intercept term

```reps```      Number of MCMC draws

## Full Output

```B```       Estimated transition matrix (posterior median)

```q```       Estimated covariance of shocks to transition equation (posterior median)

```H```       Estimated loadings (posterior median)

```R```       Estimated covariance of shocks to observation equation (poserior median, diagonal matrix)

```itc```     Intercept terms for observables

```values```  Fitted (smoothed) values for observables (including nowcasts and forecasts)

```factors``` Estimated factors (posterior median)

```Bstore```  Full posterior distribution for ```B```

```Qstore```  Full posterior distribution for ```Q```

```Hstore```  Full posterior distribution for ```H```

```Rstore```  Full posterior distribution for ```R```

```UDstore``` Forecast updates (prediction error times gain) at each period

## Estimation Details

You can find details on estimation routines and derivations in the short book *Practical Implementation of Factor Models*. [Free Download](http://srlquantitative.com)

## Examples

Examples scripts are included in the [Examples](https://github.com/srlanalytics/BDFM/tree/master/inst/Examples) file of this repository including code to run your own live nowcasts of US industrial production using the [Fred](https://fred.stlouisfed.org/) API:

![](https://github.com/srlanalytics/BDFM/blob/master/img/US_IP.png)
