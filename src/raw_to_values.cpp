// // #define ARMA_DONT_PRINT_ERRORS
// #define ARMA_NO_DEBUG
// // #define ARMA_USE_TBB_ALLOC
// #define RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR
//
// // [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>
// // [[Rcpp::interfaces(r, cpp)]]
//
//
// //==============================================================================
// //' @title anti_subset
// //'
// //' @param x values to be subset
// //' @param idx indices to remove
// //' @param nc number of columns
// //' @param nr number of rows
// //'
// //' @return
// //' @export
// //'
// //==============================================================================
// // [[Rcpp::export]]
// arma::mat anti_subset(arma::vec& x,
//                       arma::uvec& idx,
//                       arma::uword nc,
//                       arma::uword nr) {
//
//
//   // fill non value entries with NA
//   x.elem(idx).fill(NA_REAL);
//   arma::vec y = x.elem(find_finite(x)) / pow(2.0, 30);
//   arma::mat out(y);
//
//
//   if (nc == 1) {
//     return(out);
//   }
//
//
//   // generate new matrix
//   out.reshape(nc, nr);
//
//   return(out.t());
// }
//
//
// //==============================================================================
// //' @title rbr_times
// //'
// //' @param ev_tstamp
// //' @param ev_index
// //' @param ti
// //'
// //' @return
// //' @export
// //'
// //' @examples
// //==============================================================================
// // [[Rcpp::export]]
// Rcpp::DatetimeVector rbr_times(const arma::vec& ev_tstamp,
//                                const arma::uvec& ev_index,
//                                double ti) {
//
//   size_t n_ev = ev_tstamp.n_elem - 1;
//
//
//   arma::uword s, e;
//   int dif;
//
//   arma::vec a(max(ev_index) - n_ev);
//
//   for(size_t j = 0; j < n_ev; j++) {
//
//     s = ev_index[j] - 1;
//     e = ev_index[j+1] - j - 1;
//     dif = e - s;
//
//     if(dif > 1){
//       a.subvec(s, e) = ev_tstamp[j] + arma::regspace(0, dif-1) * ti;
//     }
//   }
//
//   Rcpp::DatetimeVector dt(Rcpp::NumericVector(a.begin(), a.end()));
//   dt.attr("tzone") = "UTC";
//   return(dt);
// }
//
//
// //==============================================================================
// //' @title rbr_raw_times
// //'
// //' @param raw_tstamp
// //' @param raw_index
// //' @param time_interval
// //'
// //' @return
// //' @export
// //'
// //' @examples
// //==============================================================================
// // [[Rcpp::export]]
// Rcpp::DatetimeVector rbr_raw_times(arma::vec raw_tstamp,
//                                    arma::uvec raw_index,
//                                    size_t n_raw,
//                                    size_t n_to_rem,
//                                    double time_interval) {
//
//   // output length in
//   arma::uword n_times = (n_raw - (n_to_rem * 4)) / 12;
//   arma::uword n_ev = raw_tstamp.n_elem - 1;
//   arma::uvec index(n_ev + 1);
//
//   for(size_t i = 0; i <= n_ev; i++){
//     index[i] = ((raw_index[i] / 4) - i * 2) / 3;
//   }
//
//   arma::uword s = 0;
//   arma::uword e;
//
//   double to_add;
//   double ref = 946684800.0; // 2000-01-01
//
//   arma::vec out_times(n_times);
//   out_times.fill(0);
//
//
//   for(size_t j = 0; j < n_ev; j++) {
//
//     if(index[j+1] != index[j]) {
//       e = index[j+1] - index[j];
//       to_add = raw_tstamp[j] + ref;
//       out_times.subvec(s, s+e-1) = to_add + (arma::regspace(0, e-1) * time_interval);
//       s += e;
//     }
//
//   }
//
//   Rcpp::DatetimeVector dt(Rcpp::NumericVector(out_times.begin(),
//                                               out_times.end()));
//   dt.attr("tzone") = "UTC";
//   return(dt);
// }
//
//
//
//
//
//
// //==============================================================================
// //' @title rbr_raw_to_pressure
// //'
// //' @param x
// //' @param calib
// //'
// //' @return
// //' @export
// //'
// //' @examples
// //==============================================================================
// // [[Rcpp::export]]
// arma::vec rbr_raw_to_pressure(const arma::vec& x,
//                               const arma::vec& calib) {
//
//   return(arma::polyval(arma::reverse(calib), x));
//
// }
//
//
// //==============================================================================
// //' @title rbr_raw_to_temperature
// //'
// //' @param x
// //' @param calib
// //'
// //' @return
// //' @export
// //'
// //' @examples
// //==============================================================================
// // [[Rcpp::export]]
// arma::vec rbr_raw_to_temperature(const arma::vec& x,
//                                  const arma::vec& calib) {
//
//   double k_to_c = 273.15;
//   arma::vec z = arma::log(1.0 / x - 1.0);
//
//   z = arma::polyval(arma::reverse(calib), z);
//   z = 1.0 / z - k_to_c;
//
//   return(z);
//
// }
//
//
// //==============================================================================
// //' @title rbr_temperature_correction
// //'
// //' @param pressure
// //' @param temperature
// //' @param x calibration constants
// //'
// //' @return
// //' @export
// //'
// //' @examples
// //==============================================================================
// // [[Rcpp::export]]
// arma::vec rbr_temperature_correction(const arma::vec& pressure,
//                                      const arma::vec& temperature,
//                                      const arma::vec& x) {
//
//   arma::vec tcal = temperature - x(5);
//   arma::vec out = pressure - x(0);
//   arma::uword s = 0;
//   arma::uword e = 3;
//
//   arma::vec co = arma::reverse(x.subvec(s, e));
//   co(3) = 0.0;
//
//   out -= arma::polyval(co, tcal);
//
//   out = x(0) + out / (1.0 + x(4) * tcal);
//
//   return(out);
// }
//
//
//
//
// //==============================================================================
// //' @title apply_calibration_compensation
// //'
// //' @param raw_vector
// //' @param idx
// //' @param nc
// //' @param ev_tstamp
// //' @param ev_index
// //' @param ti
// //' @param base_calib
// //' @param is_temp
// //' @param temp_comp
// //' @param pres_index
// //' @param temp_index
// //' @param names
// //' @param subset
// //'
// //' @return
// //' @export
// //'
// //' @examples
// //==============================================================================
// // [[Rcpp::export]]
// Rcpp::DataFrame apply_calibration_compensation(
//     arma::vec& raw_vector,
//     arma::uvec idx,
//     size_t nc,
//     const arma::vec& ev_tstamp,
//     const arma::uvec& ev_index,
//     double ti,
//     const arma::mat& base_calib,
//     const arma::vec& is_temp,
//     const arma::vec& temp_comp,
//     size_t pres_index,
//     size_t temp_index,
//     Rcpp::CharacterVector names,
//     int subset = 1
// ) {
//
//   size_t nr = (raw_vector.n_elem - idx.n_elem) / nc;
//
//   size_t ind = 0;
//   Rcpp::DataFrame out;
//   arma::mat raw_values;
//
//   // remove non-representative values
//   if (subset > 1) {
//       arma::uvec sub1 = arma::regspace<arma::uvec>(0, 100, nr);
//       Rcpp::IntegerVector sub2 = Rcpp::IntegerVector(sub1.begin(), sub1.end());
//
//       raw_values = anti_subset(raw_vector, idx, nc, nr).rows(sub1);
//       Rcpp::DatetimeVector dt = rbr_times(ev_tstamp, ev_index, ti);
//
//       // generate date times
//       out = Rcpp::DataFrame::create(
//           Rcpp::Named("datetime") = dt[sub2]);
//
//   } else {
//       raw_values = anti_subset(raw_vector, idx, nc, nr);
//       Rcpp::DatetimeVector dt = rbr_times(ev_tstamp, ev_index, ti);
//
//       // generate date times
//       out = Rcpp::DataFrame::create(
//           Rcpp::Named("datetime") = dt);
//   }
//
//   // Rcpp::Rcout << "The value is " << raw_values(0) << std::endl;
//   // Rcpp::Rcout << "The value is " << raw_values(nr) << std::endl;
//   // Rcpp::Rcout << "The value is " << raw_values(nr*2) << std::endl;
//
//
//   Rcpp::String nm;
//
//   // copy raw values
//   for(size_t j = 0; j < nc; j++) {
//     nm = names(ind);
//
//     out.push_back(raw_values.col(j), nm);
//
//     ind = ind + 1;
//     nm = names(ind);
//
//     if (is_temp(j)){
//       out.push_back(rbr_raw_to_temperature(raw_values.col(j),
//                                            base_calib.col(j)),
//                                            nm);
//     } else {
//       out.push_back(rbr_raw_to_pressure(raw_values.col(j),
//                                         base_calib.col(j)),
//                                         nm);
//     }
//
//     ind = ind + 1;
//   }
//
//   // do temperature compensation
//   if(pres_index != 0 & temp_index != 0) {
//     nm = names(ind);
//
//     arma::vec pressure = out[pres_index];
//     arma::vec temperature = out[temp_index];
//     out.push_back(rbr_temperature_correction(pressure,
//                                              temperature,
//                                              temp_comp),
//                                              nm);
//
//   }
//
//   return (out);
//
// }
//
//
// // [[Rcpp::export]]
// unsigned long long to_unsigned_long(std::string str){
//   std::string::size_type sz = 0;   // alias of size_t
//
//   unsigned long long ull = std::stoull (str,&sz,0);
//   return(ull);
//
// }
//
// // [[Rcpp::export]]
// Rcpp::DataFrame apply_calibration_compensation_raw(
//     arma::mat& raw_matrix,
//     const arma::mat& base_calib,
//     const arma::vec& is_temp,
//     const arma::vec& temp_comp,
//     size_t pres_index,
//     size_t temp_index,
//     arma::vec times,
//     arma::uvec time_index,
//     double time_interval,
//     size_t n_raw,
//     size_t n_to_remove
// ) {
//
//   size_t nc = raw_matrix.n_cols;
//   size_t nr = raw_matrix.n_rows;
//
//   // time column
//   Rcpp::DatetimeVector dt = rbr_raw_times(times,
//                                           time_index,
//                                           n_raw,
//                                           n_to_remove,
//                                           time_interval);
//
//   // generate date times
//   Rcpp::DataFrame out = Rcpp::DataFrame::create(
//     Rcpp::Named("datetime") = dt);
//
//
//
//   for(size_t j = 0; j < nc; j++) {
//
//     out.push_back(raw_matrix.col(j));
//
//     if (is_temp(j)){
//       out.push_back(rbr_raw_to_temperature(raw_matrix.col(j),
//                                            base_calib.col(j)));
//     } else {
//       out.push_back(rbr_raw_to_pressure(raw_matrix.col(j),
//                                         base_calib.col(j)));
//     }
//   }
//
//   // do temperature compensation
//   if(pres_index != 0 & temp_index != 0) {
//
//     arma::vec pressure = out[pres_index];
//     arma::vec temperature = out[temp_index];
//     out.push_back(rbr_temperature_correction(pressure,
//                                              temperature,
//                                              temp_comp));
//
//   }
//
//   return(out);
//
// }
//
//
// // [[Rcpp::export]]
// Rcpp::IntegerVector rsk_find_events(Rcpp::RawVector x){
//
//   Rcpp::IntegerVector z;
//
//   for(int i = 3; i < x.size(); i = i + 4) {
//     if(x[i] == 0xF7){
//       z.push_back(i - 2);
//     }
//   }
//
//   return(z);
//
// }
//
//
//
//
// // [[Rcpp::export]]
// uint32_t raw_to_4byte_unsigned(std::vector<unsigned char> x) {
//
//   uint32_t t = (x[3] << 24) | (x[2] << 16) | (x[1] << 8) | x[0];
//
//   return(t);
//
// }
// // [[Rcpp::export]]
// int32_t raw_to_4byte_signed(std::vector<unsigned char> x) {
//
//   int32_t t = (x[3] << 24) | (x[2] << 16) | (x[1] << 8) | x[0];
//
//   return(t);
//
// }
//
//
// // [[Rcpp::export]]
// arma::uvec raw_to_unsigned_vec(std::vector<unsigned char> x) {
//
//   arma::uvec out(x.size()/4);
//
//   for(size_t i=0; i < x.size(); i+=4) {
//     out(i/4) = (x[3+i] << 24) | (x[2+i] << 16) | (x[1+i] << 8) | x[0+i];
//   }
//
//   return(out);
//
// }
//
//
//
//
// /*** R
// # rv <- as.raw(c('0x10','0x01','0x10','0x01'))#,'0x00','0x10','0x00','0x01'))
// # system.time(
// #   aa <- raw_to_unsigned_vec(hex_dat)
// # )
//
// # cv <- as.character(paste0(c('0x',rv), collapse = ''))
// # readBin(rv, what = 'integer', size = 4L, endian = 'little')
// # microbenchmark::microbenchmark(
// #   a <- readBin(rv, what = 'integer', size = 4L, endian = 'little'),
// #   b <- to_unsigned_long(cv),
// #   c <- rawToChar(rv),
// #   d <- rawToBits(rv),
// #   e <- packBits(rawToBits(rv), type = 'raw'),
// #   f <- raw_slice(rv),
// #   tmp <- as.character(rv)
// # )
// # system.time(
// #   aaa <- raw_slice(rv)
// #   )
// # print(aaa)
//
// # raw_matrix <- matrix(raw_vector[-(idx+1)] / 2^30,
// #             ncol = 3,
// #             nrow = nr,
// #             byrow = TRUE)
// #
// # system.time({
// # tmp <- setDT(apply_matrix_calibration(
// #   raw_matrix,
// #   matrix(coefficients[grepl('c', key)]$value, ncol = 3),
// #   c(TRUE, FALSE, TRUE),
// #   coefficients[grepl('x', key)]$value,
// #   2*2-1,
// #   3*2-1
// # ))
// #
// # })
// # system.time(tmp <- rsk_find_events(dat))
// #
// #
// #
// # get_bin_int_4_byte(dat[tmp[2] + -24:31],n=100)
// # get_bin_int_4_byte(dat[1:28],n=100)
// #
// # times <- dat[rep(tmp, each = 4) + c(4:7)]
// # times <- get_bin_int_4_byte(times, n = length(tmp))
// #
// # system.time(b <- readBin(dat, n = 1e8, what = 'integer', size = 4L))
// # idx <- idx_to_mask(tmp, 3)
// # nr <- (length(b) - length(idx)) / 3
// # system.time(a <- matrix(b[-(idx+1)]/2^30, ncol = 3, nrow = nr, byrow = TRUE))
// # system.time(d <- anti_subset(b, idx, nc = 3, nr = nr))
// #
// #
// # to_exclude = dat[rep(tmp, each = 8) + c(0:7)]
//
// # library(data.table)
// # library(rsk)
// #
// # b <- rnorm(10)
// # a <- rnorm(10)
// # d <- rnorm(6)
// # rbr_temperature_correction(a, b, d)
// # rbr_raw_to_temperature(b, d)
// # rbr_raw_to_pressure(b, d)
// # rbr_times(c(1, 1000),c(1, 100), 1)
// #
// # system.time(((aa <- setDT(test_all(
// #     z$raw_value,
// #     as.numeric(z[['events']][['tstamp']] * 1e-3),
// #     z[['events']][['sampleIndex']],
// #     z[['continuous']][['samplingPeriod']] / 1000,
// #     matrix(c(z$base_coefficients(1), z$base_coefficients(2), z$base_coefficients(3)), ncol = 3),
// #     c(TRUE, FALSE, TRUE),
// #     4,
// #     6,
// #     z$comp_coefficients())))))
//
//
// #
// # system.time(a <- test(x,y))
// # system.time(aa <- test2(x,y))
// # system.time(b <- as.data.table(a))
// # system.time(bb <- setDT(aa))
// # tmp <- z$raw_value
// # a <- tmp
// # co <- z$base_coefficients
// # is_temp <- c(TRUE, TRUE, FALSE)
// # t1 <- system.time({
// #
// #
// #     for(i in seq_along(is_temp)) {
// #
// #         if(is_temp[i]){
// #             a[, i] <- rbr_raw_to_temperature(tmp[,i],
// #                                              z$base_coefficients(i))
// #         } else {
// #             a[,i] <- rbr_raw_to_pressure(tmp[,i],
// #                                          z$base_coefficients(i))
// #         }
// #
// #     }
// # }
// # )
// #
// # raw_parse = function(z) {
// #
// #     # read binary data
// #     raw_val <- unlist(lapply(z[['blob']][['data']], function(x) {
// #         readBin(x,
// #                 n = 100000L,
// #                 what = 'integer',
// #                 size = 4L,
// #                 signed = TRUE,
// #                 endian = 'little')
// #     }))
// #
// #     h_len  <- max(which(raw_val[1:300] == z[['time_1']]))
// #     to_rem <- c(1:h_len, (z[['events']][['sampleIndex']][-1]-1) * z[['n_chan']] + h_len + 1:2)
// #
// #
// #     raw_value <- matrix(raw_val[-to_rem] / 2^30,
// #                         ncol =  z[['n_chan']],
// #                         byrow = TRUE)
// #
// # }
// #
// # t2 <- system.time({
// #     raw_parse(z)
// # }
// # )
// #
// # t3 <- system.time({
// #     rbr_temperature_correction(a[,2],
// #                                a[,3],
// #                                z$comp_coefficients())
// # })
// # sum(c(t1['elapsed'], t2['elapsed'], t3['elapsed']))
//
//
// */
//
//
//
