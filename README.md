
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rsk

<!-- badges: start -->
<!-- badges: end -->

The goal of rsk is to facilitate fast reading of pressure and
temperature RBR files. This package is in development and exploration
stage as it is accessing a binary table of the .rsk file and the format
may change in the future. If you need a stable reader for RBR files you
may find the *oce* package useful.

Steps:

1.  read SQLite3 database tables (\*.rsk)
2.  check date times
3.  convert binary data (*downloads* table) to millivolt
4.  apply basic calibration parameters
5.  apply temperature compensation for pressure

## Installation

The development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jkennel/rsk")
```

## Example (reading a 1 GB file)

``` r
library(rsk)
fn  <- '/home/jonathankennel/Storage/data/rbr/rd45a 081871_20191118_1213.rsk'

# this reads in the data
system.time(rbr_rsk <- Rsk$new(file_name = fn, raw = TRUE))
#>    user  system elapsed 
#>   2.078   0.547   2.040
system.time(rbr_rsk <- Rsk$new(file_name = fn, raw = FALSE))
#>    user  system elapsed 
#>   8.008   0.408   7.791
system.time(rbr_oce <- oce::read.rsk(file = fn))
#>    user  system elapsed 
#>   8.782   0.456   9.238
```

## Data

The data is stored in a single *data.table* in wide format. *raw* refers
to the millivolt measurements. For newer duets, temperature is measured
onboard and at the sensor tip. *pressure\_dbar\_comp* refers to the
final temperature compensated pressure.

``` r
rbr_rsk$data
#>                      datetime   temp12   pres26   temp05
#>        1: 2019-07-10 22:58:05 29.50085 9.516360 26.57318
#>        2: 2019-07-10 22:58:06 29.50362 9.517103 26.62170
#>        3: 2019-07-10 22:58:07 29.50651 9.517877 26.59744
#>        4: 2019-07-10 22:58:08 29.48962 9.521399 26.64597
#>        5: 2019-07-10 22:58:09 29.50550 9.519480 26.64597
#>       ---                                               
#> 11292045: 2019-11-18 17:13:16 21.69363 9.492299 20.27555
#> 11292046: 2019-11-18 17:13:17 21.72822 9.491554 20.27555
#> 11292047: 2019-11-18 17:13:18 21.75391 9.495810 20.25231
#> 11292048: 2019-11-18 17:13:19 21.77008 9.496685 20.22907
#> 11292049: 2019-11-18 17:13:20 21.77443 9.490568 20.25231
```

## Plot of channels

``` r
rbr_rsk$glance()
```
