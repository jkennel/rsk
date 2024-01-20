// #define ARMA_DONT_PRINT_ERRORS
#define ARMA_NO_DEBUG
// #define ARMA_USE_TBB_ALLOC
#define RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::interfaces(r, cpp)]]


// [[Rcpp::export]]
uint32_t raw_to_4byte_unsigned(const Rcpp::RawVector x) {

  uint32_t t = (x[3] << 24) | (x[2] << 16) | (x[1] << 8) | x[0];

  return(t);

}
// [[Rcpp::export]]
int32_t raw_to_4byte_signed(const Rcpp::RawVector x, size_t index) {

  int32_t t = (x[3+index] << 24) | (x[2+index] << 16) | (x[1+index] << 8) | x[index];

  return(t);

}


// [[Rcpp::export]]
Rcpp::NumericMatrix raw_to_unsigned_vec(const Rcpp::RawVector x,
                                        const size_t header_length,
                                        const size_t n_columns,
                                        const Rcpp::IntegerVector to_remove) {

  size_t n_data = x.size() - header_length;
  size_t j=0;
  size_t k=0;

  Rcpp::NumericMatrix out(n_columns, (n_data - to_remove.length() * 4) / (4 * n_columns));

  for(size_t i = header_length; i < x.size(); i += 4) {
    if(i == to_remove[k]) {
      k += 1;
    } else {
      out[j] = ((x[3+i] << 24) | (x[2+i] << 16) | (x[1+i] << 8) | x[i] ) / 1073741824.0;
      // out(j) = raw_to_4byte_signed(x, i);
      j += 1;
    }
  }

  return(out);

}


// [[Rcpp::export]]
Rcpp::List rsk_find_events(const Rcpp::RawVector x,
                                    const size_t header_length){

  Rcpp::IntegerVector time_index;
  Rcpp::IntegerVector times;
  int ind;
  int tm;

  for(int i = 3 + header_length; i < x.size(); i = i + 4) {
    switch(x[i]) {
    case 0xF7:
    case 0xF5:
      tm = ((x[4+i] << 24) | (x[3+i] << 16) | (x[2+i] << 8) | x[1+i] );
      if(tm > 0){
        time_index.push_back(i - 3);
        times.push_back(tm);
      } else {
        Rcpp::warning("Times before 2000-01-01 were detected - be careful with data");
      }
      break;
      default:
        break;
    }
  }
  Rcpp::List times_list = Rcpp::List::create(Rcpp::Named("time_index") = time_index,
                                             Rcpp::_["times"] = times);
  return(times_list);

}

// // [[Rcpp::export]]
// Rcpp::IntegerVector rsk_find_times(const Rcpp::RawVector x,
//                                    const Rcpp::IntegerVector time_index){
//
//   size_t n_times = time_index.length();
//   Rcpp::IntegerVector times(n_times);
//
//   size_t ind;
//
//   for(size_t i = 0; i < n_times; i++) {
//     ind = time_index[i] + 4;
//     times[i] = ((x[3+ind] << 24) | (x[2+ind] << 16) | (x[1+ind] << 8) | x[0+ind] );
//   }
//
//   return(times);
//
// }

// [[Rcpp::export]]
Rcpp::IntegerVector rsk_incomplete_events(const Rcpp::IntegerVector x,
                                          const Rcpp::IntegerVector times,
                                          const size_t n_cols,
                                          const bool f5)
{

  Rcpp::IntegerVector z;

  size_t difference;
  for( size_t i = 0; i < x.size(); i++) {

    if(i == x.size() - 1) {

      z.push_back(x[i]);
      z.push_back(x[i] + 4);
    } else {

      z.push_back(x[i]);
      z.push_back(x[i] + 4);


      if(f5 & (i == 0)) {
        z.push_back(x[i] + 8);
        difference = ((x[i+1] - x[i]) - 12) % (4 * n_cols);

      } else {
        difference = ((x[i+1] - x[i]) - 8) % (4 * n_cols);
      }
      // remove incomplete events
      if(difference != 0) {
        for(int j = difference / 4; j > 0; j--) {
          z.push_back(x[i+1] - 4 * j);
        }

      }
    }
  }



  return(z);

}




//==============================================================================
//' @title rsk_raw_to_pressure
//'
//' @param x
//' @param calib
//'
//' @return
//' @export
//'
//' @examples
//==============================================================================
// [[Rcpp::export]]
arma::vec rsk_raw_to_pressure(const arma::vec x,
                              const arma::vec calib) {

  return(arma::polyval(arma::reverse(calib), x));

}


//==============================================================================
//' @title rsk_raw_to_temperature
//'
//' @param x
//' @param calib
//'
//' @return
//' @export
//'
//' @examples
//==============================================================================
// [[Rcpp::export]]
arma::vec rsk_raw_to_temperature(const arma::vec x,
                                 const arma::vec calib) {

  double k_to_c = 273.15;
  arma::vec z = arma::log(1.0 / x - 1.0);

  z = arma::polyval(arma::reverse(calib), z);
  z = 1.0 / z - k_to_c;

  return(z);

}

//==============================================================================
//' @title rsk_temperature_correction
//'
//' @param pressure
//' @param temperature
//' @param x calibration constants
//'
//' @return
//' @export
//'
//' @examples
//==============================================================================
// [[Rcpp::export]]
arma::vec rsk_temperature_correction(arma::vec out,
                                     arma::vec tcal,
                                     const arma::vec x) {

  // check length - old version had fewer coefficients
  // x2 is x5
  // x3 is x0
  if(x.n_elem == 4) {

    tcal = tcal - x(2);

    out = x(3) + (out - x(3) - x(0) * (tcal)) / (1.0 + x(1) * tcal);

  } else {
    // x0 is the calibration pressure 'Pcal' in dbar,
    // x1, x2, x3, x4  correspond directly to the constants "Kp1" through "Kp4",
    // x5 is the calibration temperature "Tcal" in °C,

    tcal = tcal - x(5);
    out = out - x(0);

    arma::uword s = 0;
    arma::uword e = 3;

    arma::vec co = arma::reverse(x.subvec(s, e));
    co(3) = 0.0;

    out -= arma::polyval(co, tcal);

    out = x(0) + out / (1.0 + x(4) * tcal);

  }

  return(out);
}



//==============================================================================
//' @title rsk_raw_times
//'
//' @param raw_tstamp
//' @param raw_index
//' @param n_times
//' @param n_columns
//' @param measurement_interval
//'
//' @return
//' @export
//'
//' @examples
//==============================================================================
// [[Rcpp::export]]
Rcpp::DatetimeVector rsk_raw_times(const Rcpp::IntegerVector raw_tstamp,
                                   const Rcpp::IntegerVector raw_index,
                                   const size_t n_times,
                                   const size_t n_columns,
                                   double measurement_interval) {

  // output length in
  size_t n_ev = raw_tstamp.length();



  // time index
  Rcpp::IntegerVector index(n_ev);

  for(size_t i = 0; i < n_ev; i++){
    index[i] = ((raw_index[i] / 4) - i * 2) / n_columns;
  }


  // start and end
  size_t s = 0;
  size_t e;


  // reference timess
  double to_add;
  double ref = 946684800.0; // 2000-01-01


  Rcpp::NumericVector out_times(n_times);
  out_times.fill(0);

  Rcpp::NumericVector temp_vec;

  if(n_ev == 1) {
    temp_vec = Rcpp::seq(0, n_times-1); // create integer first then convert
    out_times = raw_tstamp[0] + ref + temp_vec * measurement_interval;
  }


  for(size_t j = 0; j < n_ev - 1; j++) {

    if(index[j+1] != index[j]) {

      e = index[j+1] - index[j];

      to_add = raw_tstamp[j] + ref;
      temp_vec = Rcpp::seq(0, e-1); // create integer first then convert
      temp_vec = to_add + Rcpp::as<Rcpp::NumericVector>(temp_vec) * measurement_interval;

      out_times[Rcpp::Range(s, s + e - 1)] = temp_vec;
      s += e;

    }

  }

  Rcpp::DatetimeVector dt(out_times);

  dt.attr("tzone") = "UTC";
  return(dt);
}


// // This needs review
// // [[Rcpp::export]]
// size_t get_header_length(Rcpp::RawVector x, size_t max_length = 65535) {
//
//   size_t header_length;
//
//   for(size_t i = 3; i < x.size(); i = i + 4) {
//     if(x[i] == 0xF7){
//       header_length = i - 3 ;
//       break;
//     }
//   }
//
//   return(header_length);
//
// }


// This needs review if other formats are possible
// [[Rcpp::export]]
size_t get_header_length(Rcpp::RawVector x) {

  size_t read_0 = ( x[0] );
  size_t header_length;

  if(read_0 == 0) {
    header_length = ( (x[1] << 8) | x[0] );
  } else {
    header_length = ( (x[1+7] << 8) | x[0+7] );
  }



  // if(x[header_length + 3] == 0xF5) {
  //   Rcpp::Rcout << "ncol : " << 1 << "\n";
  //
  //   // header_length -= ;
  // }
  //
  // if(header_length != 0) {
  //   for(size_t i = 3; i < max_length; i = i + 4) {
  //     if(x[i] == 0xF7){
  //       header_length = i - 3 ;
  //       break;
  //     }
  //   }
  // }


  return(header_length);

}



// This function attempts to only interpret the data in the binary format
// excluding all the header info
// [[Rcpp::export]]
Rcpp::DataFrame rsk_read_bin(Rcpp::RawVector x,
                             Rcpp::LogicalVector is_temp,
                             Rcpp::NumericMatrix base_calib,
                             Rcpp::NumericVector temp_comp,
                             double measurement_interval,
                             size_t pressure_index,
                             size_t temperature_index,
                             bool keep_raw) {


  // check temperature columns
  size_t n_channels = is_temp.length();
  bool f5 = false;

  // where does the data start
  size_t header_length = get_header_length(x);

  if(x[header_length + 3] == 0xF5) {
    f5 = true;
  }

  Rcpp::List tms = rsk_find_events(x, header_length);
  Rcpp::IntegerVector time_index = tms[0];
  Rcpp::IntegerVector times = tms[1];
  // Rcpp::IntegerVector times      = rsk_find_times(x, time_index);
  Rcpp::IntegerVector to_remove  = rsk_incomplete_events(time_index, times, n_channels, f5);

  Rcpp::NumericMatrix raw_matrix = raw_to_unsigned_vec(x, header_length, n_channels, to_remove);

  Rcpp::DataFrame out_df = Rcpp::DataFrame::create();

  size_t n_raw = x.length();

  // time column
  out_df.push_back(rsk_raw_times(times,
                                 time_index,
                                 raw_matrix.ncol(),
                                 raw_matrix.nrow(),
                                 measurement_interval));

  // Rcpp::Rcout << "The size of datetime : " << out_df.nrows() << "\n";

  // keep the millivolt readings
  if (keep_raw) {
    pressure_index = pressure_index * 2;
    temperature_index = temperature_index * 2;
  }

  // pressure and temperature columns are handled differently
  for(size_t j = 0; j < n_channels; j++) {

    // raw millivolt readings
    if (keep_raw) {
      out_df.push_back(raw_matrix.row(j));
    }

    // calibrated temperature
    if (is_temp(j)){
      out_df.push_back(rsk_raw_to_temperature(raw_matrix.row(j),
                                              base_calib.row(j)));
      // calibrated pressure
    } else {
      out_df.push_back(rsk_raw_to_pressure(raw_matrix.row(j),
                                           base_calib.row(j)));
    }
  }
  // Rcpp::Rcout << "The size of time : " << temperature_index << "\n";
  // Rcpp::Rcout << "The size of time : " << pressure_index << "\n";

  // do temperature compensation
  if(pressure_index != 0 & temperature_index != 0) {

    arma::vec pressure = out_df[pressure_index];
    arma::vec temperature = out_df[temperature_index];
    out_df.push_back(rsk_temperature_correction(pressure,
                                                temperature,
                                                temp_comp));


  }


  return(out_df);
}


/*** R
# library(data.table)
# x <- '/home/jonathankennel/Desktop/rd45a 081871_20201026_1223.bin'
#
#
# system.time({
#   hex_dat <- readBin(x, n = 2e8, what = 'raw')
#   n   <- setDT(rsk_read_bin(hex_dat,
#                c(TRUE, FALSE, TRUE),
#                t(base_calib),
#                base_comp))
# })


# readBin(hex_dat[661:664], what = 'integer', size = 4)
# readBin(hex_dat[665:668], what = 'integer', size = 4)
# readBin(hex_dat[668:671], what = 'integer', size = 4)
*/
