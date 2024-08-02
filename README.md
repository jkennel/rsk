
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
fn = "../../r_scratch/c2 081870_20240501_1603.rsk"
# this reads in the data
system.time(rbr_rsk <- Rsk$new(file_name = fn, raw = TRUE))
#>    user  system elapsed 
#>   1.574   0.419   1.653
system.time(rbr_rsk <- Rsk$new(file_name = fn, raw = FALSE))
#>    user  system elapsed 
#>   5.495   0.341   5.162
system.time(rbr_oce <- oce::read.rsk(file = fn))
#>    user  system elapsed 
#>  11.407   0.481  11.890
```

## Data

The data is stored in a single *data.table* in wide format. *raw* refers
to the millivolt measurements. For newer duets, temperature is measured
onboard and at the sensor tip. *pressure_dbar_comp* refers to the final
temperature compensated pressure.

``` r
rbr_rsk$data
#> Key: <datetime>
#>                      datetime   temp12   pres26   temp05
#>                        <POSc>    <num>    <num>    <num>
#>        1: 2023-12-08 19:51:38 21.83675 9.490034 24.72569
#>        2: 2023-12-08 19:51:39 21.79042 9.495672 24.60640
#>        3: 2023-12-08 19:51:40 21.70013 9.492937 24.63025
#>        4: 2023-12-08 19:51:41 21.67313 9.492478 24.65410
#>        5: 2023-12-08 19:51:42 21.67435 9.491035 24.63025
#>       ---                                               
#> 12528671: 2024-05-01 20:02:48 19.21357 9.427519 18.56454
#> 12528672: 2024-05-01 20:02:49 19.19873 9.428241 18.61078
#> 12528673: 2024-05-01 20:02:50 19.18560 9.427547 18.56454
#> 12528674: 2024-05-01 20:02:51 19.18334 9.431602 18.58766
#> 12528675: 2024-05-01 20:02:52 19.18289 9.430419 18.58766
```

## Plot of channels

``` r
rbr_rsk$glance()
```
