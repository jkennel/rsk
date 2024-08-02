
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rsk

<!-- badges: start -->
<!-- badges: end -->

The goal of rsk is to facilitate fast reading of pressure and
temperature RBR files. This package is in development and exploration
stage as it is accessing a binary table of the .rsk file and the format
may change in the future. I am not associated with RBR Ltd., and the
focus of this package has been on the compact solo and duet devices for
pressure and temperature. If you need a more stable reader for RBR files
you may find the *oce* package more useful.

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
library(oce)
#> Loading required package: gsw
fn = "../../r_scratch/ELR1-R1_Port16_208651_20240626_1035.rsk"
# this reads in the data
system.time(rbr_rsk <- Rsk$new(file_name = fn, raw = TRUE))
#>    user  system elapsed 
#>   0.766   0.195   0.847
system.time(rbr_rsk <- Rsk$new(file_name = fn, raw = FALSE))
#>    user  system elapsed 
#>   2.738   0.148   2.443
system.time(rbr_oce <- oce::read.rsk(file = fn))
#>    user  system elapsed 
#>   5.254   0.290   5.546
```

## Data

The data is stored in a single *data.table* in wide format. *raw* refers
to the millivolt measurements. For newer duets, temperature is measured
onboard and at the sensor tip. *pressure_dbar_comp* refers to the final
temperature compensated pressure.

``` r
rbr_rsk$data
#> Key: <datetime>
#>                     datetime   pres25   temp05
#>                       <POSc>    <num>    <num>
#>       1: 2024-03-31 12:00:00 9.745065 21.08618
#>       2: 2024-03-31 12:00:01 9.744286 21.08618
#>       3: 2024-03-31 12:00:02 9.744805 21.08618
#>       4: 2024-03-31 12:00:03 9.744675 21.08618
#>       5: 2024-03-31 12:00:04 9.744140 21.10955
#>      ---                                      
#> 7526159: 2024-06-26 14:35:58 9.640644 27.28401
#> 7526160: 2024-06-26 14:35:59 9.639737 27.33291
#> 7526161: 2024-06-26 14:36:00 9.640282 27.30846
#> 7526162: 2024-06-26 14:36:01 9.640644 27.28401
#> 7526163: 2024-06-26 14:36:02 9.642498 27.03984
```

## Plot of channels

``` r
rbr_rsk$glance()
```
