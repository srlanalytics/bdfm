# BDFM
Bayesian dynamic factor model estimation in R.

## Description

BDFM estimates dynamic factor models by simulation using the Durbin and Koopman (2012) disturbance smoother and returns estimated factors, predicted values of observables, posterior distributions for predicted values (optional), and forecast updates by series (optional). Maximum likelihood estimation via Watson and Engle (1983) is also supported. Input data may be noisy, have missing values, or "ragged edges" due to different start or end dates, but must be uniform frequency for this version. 

## Installation

```t
remotes::install_github("srlanalytics/BDFM")
```

## Basic Use

Estimate a model with two factors and two lags in the transition equation
```r
library(BDFM)
Est <- BDFM(Y = Data, m = 2, p = 2)
```
Input data is any data type that can be coverted by ```as.matrix() ``` and should index time in rows and series in columns.

Quick time series plot estimated factors
```r
ts.plot(Est$X)
```
Quick time series plot of fitted values (including forecasts and nowcasts) for first series in input data
```r
ts.plot(Est$Y[,1])
```

## Full Input

```r
Est <- (Y, m, p, FC = 0, Bp = NULL, lam_B = 1, Hp = NULL, lam_H = 1, nu_q = NULL, nu_r = NULL, ID = "PC_full", ITC = T, reps = 1000)
```

```Y```     Input data

```m```     Number of factors

```p```     Number of lags in transition equation

```FC```    Number of periods ahead to forecast

```Bp```    Prior on transtion equation (default is to shrink values towards zero)

```lam_B``` Prior tightness on transition equation

```Hp```    Prior on loadings/observation equation (default is zero)

```lam_H``` Prior tightness on loadings

```nu_q```  Prior "degrees of freedom" for transition equation (i.e. shrink covariance of shocks towards zero)

```nu_r```  Prior "degrees of freedom" for observation equation (i.e. shrink covariance of shocks towards zero)

```ID```    Model identification technique. 'PC_full' identifies on principal components using all series. 'PC-sub' identifies on principal components using a subset of the data that maximizes the number of observations when many observations are missing.

```ITC``` T/F estimate intercept term

```reps``` Number of MCMC draws

## Full Output

```B``` Estimated transition matrix (posterior median)

```q``` Estimated covariance of shocks to transition equation (posterior median)

```H``` Estimated loadings (posterior median)

```R``` Estimated covariance of shocks to observation equation (poserior median, diagonal matrix)

```itc``` Intercept terms for observables

```Ys``` Fitted (smoothed) values for observables (including nowcasts and forecasts)

```X``` Estimated factors (posterior median)

```Bstore``` Full posterior distribution for ```B```

```Qstore``` Full posterior distribution for ```Q```

```Hstore``` Full posterior distribution for ```H```

```Rstore``` Full posterior distribution for ```R```

```UDstore``` Forecast updates (prediction error times gain) at each period

```Ydist```  Distribution for end of sample fitted values

## Estimation Details

You can find details on estimation routines and derivations in the short book *Practical Implementation of Factor Models*. [Free Download](http://srlquantitative.com)

## Examples

Examples scripts are included in the [Examples](https://github.com/srlanalytics/BDFM/tree/master/inst/Examples) file of this repository including code to run your own live nowcasts of US industrial production using the [Fred](https://fred.stlouisfed.org/) API:

![](https://github.com/srlanalytics/BDFM/blob/master/img/US_IP.png)
