
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rsk

<!-- badges: start -->

<!-- badges: end -->

The goal of rsk is to facilitate fast reading of pressure and
temperature RBR files. This package is in the beginning of development
stages.

Steps:

1.  read SQLite3 database tables (\*.rsk)
2.  check date times
3.  convert binary data (*downloads* table) to millivolt
4.  apply basic calibration parameters
5.  apply temperature compensation for pressure

## Installation

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jkennel/rsk")
```

## Example

``` r

library(rsk)
fn  <- '/home/jonathankennel/Storage/data/rbr/rd45a 081871_20191118_1213.rsk'

# this reads in the data
rbr <- Rsk$new(file_name = fn)
```

## Data

The data is stored in a single dataset in wide format. *raw* refers to
the millivolt measurements. For newer duets, temperature is measured
onboard and at the sensor tip. *pressure\_dbar\_comp* refers to the
final temperature compensated pressure.

``` r
rbr$data
#>                      datetime temperature_raw temperature pressure_raw pressure
#>        1: 2019-07-10 22:58:05       0.3053592    29.50085   0.05286860 9.553530
#>        2: 2019-07-10 22:58:06       0.3053336    29.50362   0.05287504 9.555076
#>        3: 2019-07-10 22:58:07       0.3053068    29.50651   0.05287659 9.555449
#>        4: 2019-07-10 22:58:08       0.3054631    29.48962   0.05289459 9.559774
#>        5: 2019-07-10 22:58:09       0.3053162    29.50550   0.05288661 9.557855
#>       ---                                                                      
#> 11297712: 2019-11-18 17:13:16       0.3838298    21.69363   0.05234885 9.428646
#> 11297713: 2019-11-18 17:13:17       0.3834562    21.72822   0.05234575 9.427901
#> 11297714: 2019-11-18 17:13:18       0.3831789    21.75391   0.05236197 9.431796
#> 11297715: 2019-11-18 17:13:19       0.3830044    21.77008   0.05236411 9.432312
#> 11297716: 2019-11-18 17:13:20       0.3829576    21.77443   0.05234015 9.426555
#>           temperature_onboard_raw temperature_onboard pressure_dbar_comp
#>        1:               0.4533691            26.57318           9.516360
#>        2:               0.4528809            26.62170           9.517103
#>        3:               0.4531250            26.59744           9.517877
#>        4:               0.4526367            26.64597           9.521399
#>        5:               0.4526367            26.64597           9.519480
#>       ---                                                               
#> 11297712:               0.5183105            20.27555           9.492299
#> 11297713:               0.5183105            20.27555           9.491554
#> 11297714:               0.5185547            20.25231           9.495810
#> 11297715:               0.5187988            20.22907           9.496685
#> 11297716:               0.5185547            20.25231           9.490568
```

## Plot of channels

``` r
rbr$glance()
```
