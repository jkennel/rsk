// #define ARMA_DONT_PRINT_ERRORS
#define ARMA_NO_DEBUG
// #define ARMA_USE_TBB_ALLOC
#define RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::interfaces(r, cpp)]]



// [[Rcpp::export]]
arma::mat anti_subset(arma::vec& x,
                      arma::uvec& idx,
                      arma::uword nc,
                      arma::uword nr) {


  // fill non value entries with NA
  x.elem(idx).fill(NA_REAL);
  arma::vec y = x.elem(find_finite(x)) / pow(2.0, 30);
  arma::mat out(y);


  if (nc == 1) {
    return(out);
  }


  // generate new matrix
  out.reshape(nc, nr);

  return(out.t());
}



//' rbr_times
//'
//' @param ev_tstamp
//' @param ev_index
//' @param ti
//'
//' @return
//' @export
//'
//' @examples
// [[Rcpp::export]]
Rcpp::DatetimeVector rbr_times(const arma::vec& ev_tstamp,
                               const arma::uvec& ev_index,
                               double ti) {

  size_t n_ev = ev_tstamp.n_elem - 1;
  arma::uword s, e;
  int dif;

  arma::vec a(max(ev_index) - n_ev);

  for(size_t j = 0; j < n_ev; j++) {
    s = ev_index[j] - 1;
    e = ev_index[j+1] - j - 1;
    dif = e - s;
    // Rcpp::Rcout << "The value is " << s << std::endl;
    // Rcpp::Rcout << "The value is " << e << std::endl;
    a.subvec(s, e) = ev_tstamp[j] + arma::linspace(0, dif-1, dif) * ti;
  }

  Rcpp::DatetimeVector dt(Rcpp::NumericVector(a.begin(), a.end()));
  dt.attr("tzone") = "UTC";
  return(dt);
}


//' rbr_raw_to_pressure
//'
//' @param x
//' @param calib
//'
//' @return
//' @export
//'
//' @examples
// [[Rcpp::export]]
arma::vec rbr_raw_to_pressure(const arma::vec& x,
                              const arma::vec& calib) {

  return(arma::polyval(arma::reverse(calib), x));

}


//' rbr_raw_to_temperature
//'
//' @param x
//' @param calib
//'
//' @return
//' @export
//'
//' @examples
// [[Rcpp::export]]
arma::vec rbr_raw_to_temperature(const arma::vec& x,
                                 const arma::vec& calib) {

  double k_to_c = 273.15;
  arma::vec z = arma::log(1.0 / x - 1.0);

  z = arma::polyval(arma::reverse(calib), z);
  z = 1.0 / z - k_to_c;

  return(z);

}


//' rbr_temperature_correction
//'
//' @param pressure
//' @param temperature
//' @param x calibration constants
//'
//' @return
//' @export
//'
//' @examples
// [[Rcpp::export]]
arma::vec rbr_temperature_correction(const arma::vec& pressure,
                                     const arma::vec& temperature,
                                     const arma::vec& x) {

  arma::vec tcal = temperature - x(5);
  arma::vec out = pressure - x(0);
  arma::uword s = 0;
  arma::uword e = 3;

  arma::vec co = arma::reverse(x.subvec(s, e));
  co(3) = 0.0;

  out -= arma::polyval(co, tcal);

  out = x(0) + out / (1.0 + x(4) * tcal);

  return(out);
}





// [[Rcpp::export]]
Rcpp::DataFrame test_all(arma::vec& raw_vector,
                         arma::uvec idx,
                         size_t nc,
                         const arma::vec& ev_tstamp,
                         const arma::uvec& ev_index,
                         double ti,
                         const arma::mat& base_calib,
                         const arma::vec& is_temp,
                         const arma::vec& temp_comp,
                         size_t pres_index,
                         size_t temp_index,
                         Rcpp::CharacterVector names,
                         int subset = 1
) {

  size_t nr = (raw_vector.n_elem - idx.n_elem) / nc;

  size_t ind = 0;
  Rcpp::DataFrame out;
  arma::mat raw_values;

  // remove non-representative values
  if (subset > 1) {
      arma::uvec sub1 = arma::regspace<arma::uvec>(0, 100, nr);
      Rcpp::IntegerVector sub2 = Rcpp::IntegerVector(sub1.begin(), sub1.end());

      raw_values = anti_subset(raw_vector, idx, nc, nr).rows(sub1);
      Rcpp::DatetimeVector dt = rbr_times(ev_tstamp, ev_index, ti);
      // generate date times
      out = Rcpp::DataFrame::create(
          Rcpp::Named("datetime") = dt[sub2]);

  } else {
      raw_values = anti_subset(raw_vector, idx, nc, nr);
      Rcpp::DatetimeVector dt = rbr_times(ev_tstamp, ev_index, ti);

      // generate date times
      out = Rcpp::DataFrame::create(
          Rcpp::Named("datetime") = dt);
  }

  // Rcpp::Rcout << "The value is " << raw_values(0) << std::endl;
  // Rcpp::Rcout << "The value is " << raw_values(nr) << std::endl;
  // Rcpp::Rcout << "The value is " << raw_values(nr*2) << std::endl;


  Rcpp::String nm;

  // copy raw values
  for(size_t j = 0; j < nc; j++) {
    nm = names(ind);

    out.push_back(raw_values.col(j), nm);

    ind = ind + 1;
    nm = names(ind);

    if (is_temp(j)){
      out.push_back(rbr_raw_to_temperature(raw_values.col(j),
                                           base_calib.col(j)),
                                           nm);
    } else {
      out.push_back(rbr_raw_to_pressure(raw_values.col(j),
                                        base_calib.col(j)),
                                        nm);
    }

    ind = ind + 1;
  }

  // do temperature compensation
  if(pres_index != 0 & temp_index != 0) {
    nm = names(ind);

    arma::vec pressure = out[pres_index];
    arma::vec temperature = out[temp_index];
    out.push_back(rbr_temperature_correction(pressure,
                                             temperature,
                                             temp_comp),
                                             nm);

  }

  return (out);

}

/*** R

# library(data.table)
# library(rsk)
#
# b <- rnorm(10)
# a <- rnorm(10)
# d <- rnorm(6)
# rbr_temperature_correction(a, b, d)
# rbr_raw_to_temperature(b, d)
# rbr_raw_to_pressure(b, d)
# rbr_times(c(1, 1000),c(1, 100), 1)
#
# system.time(((aa <- setDT(test_all(
#     z$raw_value,
#     as.numeric(z[['events']][['tstamp']] * 1e-3),
#     z[['events']][['sampleIndex']],
#     z[['continuous']][['samplingPeriod']] / 1000,
#     matrix(c(z$base_coefficients(1), z$base_coefficients(2), z$base_coefficients(3)), ncol = 3),
#     c(TRUE, FALSE, TRUE),
#     4,
#     6,
#     z$comp_coefficients())))))


#
# system.time(a <- test(x,y))
# system.time(aa <- test2(x,y))
# system.time(b <- as.data.table(a))
# system.time(bb <- setDT(aa))
# tmp <- z$raw_value
# a <- tmp
# co <- z$base_coefficients
# is_temp <- c(TRUE, TRUE, FALSE)
# t1 <- system.time({
#
#
#     for(i in seq_along(is_temp)) {
#
#         if(is_temp[i]){
#             a[, i] <- rbr_raw_to_temperature(tmp[,i],
#                                              z$base_coefficients(i))
#         } else {
#             a[,i] <- rbr_raw_to_pressure(tmp[,i],
#                                          z$base_coefficients(i))
#         }
#
#     }
# }
# )
#
# raw_parse = function(z) {
#
#     # read binary data
#     raw_val <- unlist(lapply(z[['blob']][['data']], function(x) {
#         readBin(x,
#                 n = 100000L,
#                 what = 'integer',
#                 size = 4L,
#                 signed = TRUE,
#                 endian = 'little')
#     }))
#
#     h_len  <- max(which(raw_val[1:300] == z[['time_1']]))
#     to_rem <- c(1:h_len, (z[['events']][['sampleIndex']][-1]-1) * z[['n_chan']] + h_len + 1:2)
#
#
#     raw_value <- matrix(raw_val[-to_rem] / 2^30,
#                         ncol =  z[['n_chan']],
#                         byrow = TRUE)
#
# }
#
# t2 <- system.time({
#     raw_parse(z)
# }
# )
#
# t3 <- system.time({
#     rbr_temperature_correction(a[,2],
#                                a[,3],
#                                z$comp_coefficients())
# })
# sum(c(t1['elapsed'], t2['elapsed'], t3['elapsed']))


*/



