// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/rsk.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// raw_to_4byte_unsigned
uint32_t raw_to_4byte_unsigned(const Rcpp::RawVector x);
static SEXP _rsk_raw_to_4byte_unsigned_try(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Rcpp::RawVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(raw_to_4byte_unsigned(x));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _rsk_raw_to_4byte_unsigned(SEXP xSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_rsk_raw_to_4byte_unsigned_try(xSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// raw_to_4byte_signed
int32_t raw_to_4byte_signed(const Rcpp::RawVector x, size_t index);
static SEXP _rsk_raw_to_4byte_signed_try(SEXP xSEXP, SEXP indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Rcpp::RawVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< size_t >::type index(indexSEXP);
    rcpp_result_gen = Rcpp::wrap(raw_to_4byte_signed(x, index));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _rsk_raw_to_4byte_signed(SEXP xSEXP, SEXP indexSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_rsk_raw_to_4byte_signed_try(xSEXP, indexSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// raw_to_unsigned_vec
Rcpp::NumericMatrix raw_to_unsigned_vec(const Rcpp::RawVector x, const size_t header_length, const size_t n_columns, const Rcpp::IntegerVector to_remove);
static SEXP _rsk_raw_to_unsigned_vec_try(SEXP xSEXP, SEXP header_lengthSEXP, SEXP n_columnsSEXP, SEXP to_removeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Rcpp::RawVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const size_t >::type header_length(header_lengthSEXP);
    Rcpp::traits::input_parameter< const size_t >::type n_columns(n_columnsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type to_remove(to_removeSEXP);
    rcpp_result_gen = Rcpp::wrap(raw_to_unsigned_vec(x, header_length, n_columns, to_remove));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _rsk_raw_to_unsigned_vec(SEXP xSEXP, SEXP header_lengthSEXP, SEXP n_columnsSEXP, SEXP to_removeSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_rsk_raw_to_unsigned_vec_try(xSEXP, header_lengthSEXP, n_columnsSEXP, to_removeSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// rsk_find_events
Rcpp::IntegerVector rsk_find_events(const Rcpp::RawVector x, const size_t header_length);
static SEXP _rsk_rsk_find_events_try(SEXP xSEXP, SEXP header_lengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Rcpp::RawVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const size_t >::type header_length(header_lengthSEXP);
    rcpp_result_gen = Rcpp::wrap(rsk_find_events(x, header_length));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _rsk_rsk_find_events(SEXP xSEXP, SEXP header_lengthSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_rsk_rsk_find_events_try(xSEXP, header_lengthSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// rsk_find_times
Rcpp::IntegerVector rsk_find_times(const Rcpp::RawVector x, const Rcpp::IntegerVector time_index);
static SEXP _rsk_rsk_find_times_try(SEXP xSEXP, SEXP time_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Rcpp::RawVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type time_index(time_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(rsk_find_times(x, time_index));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _rsk_rsk_find_times(SEXP xSEXP, SEXP time_indexSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_rsk_rsk_find_times_try(xSEXP, time_indexSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// rsk_incomplete_events
Rcpp::IntegerVector rsk_incomplete_events(const Rcpp::IntegerVector x, const size_t n_cols, const bool f5);
static SEXP _rsk_rsk_incomplete_events_try(SEXP xSEXP, SEXP n_colsSEXP, SEXP f5SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const size_t >::type n_cols(n_colsSEXP);
    Rcpp::traits::input_parameter< const bool >::type f5(f5SEXP);
    rcpp_result_gen = Rcpp::wrap(rsk_incomplete_events(x, n_cols, f5));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _rsk_rsk_incomplete_events(SEXP xSEXP, SEXP n_colsSEXP, SEXP f5SEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_rsk_rsk_incomplete_events_try(xSEXP, n_colsSEXP, f5SEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// rsk_raw_to_pressure
arma::vec rsk_raw_to_pressure(const arma::vec x, const arma::vec calib);
static SEXP _rsk_rsk_raw_to_pressure_try(SEXP xSEXP, SEXP calibSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type calib(calibSEXP);
    rcpp_result_gen = Rcpp::wrap(rsk_raw_to_pressure(x, calib));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _rsk_rsk_raw_to_pressure(SEXP xSEXP, SEXP calibSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_rsk_rsk_raw_to_pressure_try(xSEXP, calibSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// rsk_raw_to_temperature
arma::vec rsk_raw_to_temperature(const arma::vec x, const arma::vec calib);
static SEXP _rsk_rsk_raw_to_temperature_try(SEXP xSEXP, SEXP calibSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type calib(calibSEXP);
    rcpp_result_gen = Rcpp::wrap(rsk_raw_to_temperature(x, calib));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _rsk_rsk_raw_to_temperature(SEXP xSEXP, SEXP calibSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_rsk_rsk_raw_to_temperature_try(xSEXP, calibSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// rsk_temperature_correction
arma::vec rsk_temperature_correction(arma::vec out, arma::vec tcal, const arma::vec x);
static SEXP _rsk_rsk_temperature_correction_try(SEXP outSEXP, SEXP tcalSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::vec >::type out(outSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tcal(tcalSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rsk_temperature_correction(out, tcal, x));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _rsk_rsk_temperature_correction(SEXP outSEXP, SEXP tcalSEXP, SEXP xSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_rsk_rsk_temperature_correction_try(outSEXP, tcalSEXP, xSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// rsk_raw_times
Rcpp::DatetimeVector rsk_raw_times(const Rcpp::IntegerVector raw_tstamp, const Rcpp::IntegerVector raw_index, const size_t n_times, const size_t n_columns, double measurement_interval);
static SEXP _rsk_rsk_raw_times_try(SEXP raw_tstampSEXP, SEXP raw_indexSEXP, SEXP n_timesSEXP, SEXP n_columnsSEXP, SEXP measurement_intervalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type raw_tstamp(raw_tstampSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type raw_index(raw_indexSEXP);
    Rcpp::traits::input_parameter< const size_t >::type n_times(n_timesSEXP);
    Rcpp::traits::input_parameter< const size_t >::type n_columns(n_columnsSEXP);
    Rcpp::traits::input_parameter< double >::type measurement_interval(measurement_intervalSEXP);
    rcpp_result_gen = Rcpp::wrap(rsk_raw_times(raw_tstamp, raw_index, n_times, n_columns, measurement_interval));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _rsk_rsk_raw_times(SEXP raw_tstampSEXP, SEXP raw_indexSEXP, SEXP n_timesSEXP, SEXP n_columnsSEXP, SEXP measurement_intervalSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_rsk_rsk_raw_times_try(raw_tstampSEXP, raw_indexSEXP, n_timesSEXP, n_columnsSEXP, measurement_intervalSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// get_header_length
size_t get_header_length(Rcpp::RawVector x);
static SEXP _rsk_get_header_length_try(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(get_header_length(x));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _rsk_get_header_length(SEXP xSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_rsk_get_header_length_try(xSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// rsk_read_bin
Rcpp::DataFrame rsk_read_bin(Rcpp::RawVector x, Rcpp::LogicalVector is_temp, Rcpp::NumericMatrix base_calib, Rcpp::NumericVector temp_comp, double measurement_interval, size_t pressure_index, size_t temperature_index, bool keep_raw);
static SEXP _rsk_rsk_read_bin_try(SEXP xSEXP, SEXP is_tempSEXP, SEXP base_calibSEXP, SEXP temp_compSEXP, SEXP measurement_intervalSEXP, SEXP pressure_indexSEXP, SEXP temperature_indexSEXP, SEXP keep_rawSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::RawVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type is_temp(is_tempSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type base_calib(base_calibSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type temp_comp(temp_compSEXP);
    Rcpp::traits::input_parameter< double >::type measurement_interval(measurement_intervalSEXP);
    Rcpp::traits::input_parameter< size_t >::type pressure_index(pressure_indexSEXP);
    Rcpp::traits::input_parameter< size_t >::type temperature_index(temperature_indexSEXP);
    Rcpp::traits::input_parameter< bool >::type keep_raw(keep_rawSEXP);
    rcpp_result_gen = Rcpp::wrap(rsk_read_bin(x, is_temp, base_calib, temp_comp, measurement_interval, pressure_index, temperature_index, keep_raw));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _rsk_rsk_read_bin(SEXP xSEXP, SEXP is_tempSEXP, SEXP base_calibSEXP, SEXP temp_compSEXP, SEXP measurement_intervalSEXP, SEXP pressure_indexSEXP, SEXP temperature_indexSEXP, SEXP keep_rawSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_rsk_rsk_read_bin_try(xSEXP, is_tempSEXP, base_calibSEXP, temp_compSEXP, measurement_intervalSEXP, pressure_indexSEXP, temperature_indexSEXP, keep_rawSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}

// validate (ensure exported C++ functions exist before calling them)
static int _rsk_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("uint32_t(*raw_to_4byte_unsigned)(const Rcpp::RawVector)");
        signatures.insert("int32_t(*raw_to_4byte_signed)(const Rcpp::RawVector,size_t)");
        signatures.insert("Rcpp::NumericMatrix(*raw_to_unsigned_vec)(const Rcpp::RawVector,const size_t,const size_t,const Rcpp::IntegerVector)");
        signatures.insert("Rcpp::IntegerVector(*rsk_find_events)(const Rcpp::RawVector,const size_t)");
        signatures.insert("Rcpp::IntegerVector(*rsk_find_times)(const Rcpp::RawVector,const Rcpp::IntegerVector)");
        signatures.insert("Rcpp::IntegerVector(*rsk_incomplete_events)(const Rcpp::IntegerVector,const size_t,const bool)");
        signatures.insert("arma::vec(*rsk_raw_to_pressure)(const arma::vec,const arma::vec)");
        signatures.insert("arma::vec(*rsk_raw_to_temperature)(const arma::vec,const arma::vec)");
        signatures.insert("arma::vec(*rsk_temperature_correction)(arma::vec,arma::vec,const arma::vec)");
        signatures.insert("Rcpp::DatetimeVector(*rsk_raw_times)(const Rcpp::IntegerVector,const Rcpp::IntegerVector,const size_t,const size_t,double)");
        signatures.insert("size_t(*get_header_length)(Rcpp::RawVector)");
        signatures.insert("Rcpp::DataFrame(*rsk_read_bin)(Rcpp::RawVector,Rcpp::LogicalVector,Rcpp::NumericMatrix,Rcpp::NumericVector,double,size_t,size_t,bool)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _rsk_RcppExport_registerCCallable() { 
    R_RegisterCCallable("rsk", "_rsk_raw_to_4byte_unsigned", (DL_FUNC)_rsk_raw_to_4byte_unsigned_try);
    R_RegisterCCallable("rsk", "_rsk_raw_to_4byte_signed", (DL_FUNC)_rsk_raw_to_4byte_signed_try);
    R_RegisterCCallable("rsk", "_rsk_raw_to_unsigned_vec", (DL_FUNC)_rsk_raw_to_unsigned_vec_try);
    R_RegisterCCallable("rsk", "_rsk_rsk_find_events", (DL_FUNC)_rsk_rsk_find_events_try);
    R_RegisterCCallable("rsk", "_rsk_rsk_find_times", (DL_FUNC)_rsk_rsk_find_times_try);
    R_RegisterCCallable("rsk", "_rsk_rsk_incomplete_events", (DL_FUNC)_rsk_rsk_incomplete_events_try);
    R_RegisterCCallable("rsk", "_rsk_rsk_raw_to_pressure", (DL_FUNC)_rsk_rsk_raw_to_pressure_try);
    R_RegisterCCallable("rsk", "_rsk_rsk_raw_to_temperature", (DL_FUNC)_rsk_rsk_raw_to_temperature_try);
    R_RegisterCCallable("rsk", "_rsk_rsk_temperature_correction", (DL_FUNC)_rsk_rsk_temperature_correction_try);
    R_RegisterCCallable("rsk", "_rsk_rsk_raw_times", (DL_FUNC)_rsk_rsk_raw_times_try);
    R_RegisterCCallable("rsk", "_rsk_get_header_length", (DL_FUNC)_rsk_get_header_length_try);
    R_RegisterCCallable("rsk", "_rsk_rsk_read_bin", (DL_FUNC)_rsk_rsk_read_bin_try);
    R_RegisterCCallable("rsk", "_rsk_RcppExport_validate", (DL_FUNC)_rsk_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_rsk_raw_to_4byte_unsigned", (DL_FUNC) &_rsk_raw_to_4byte_unsigned, 1},
    {"_rsk_raw_to_4byte_signed", (DL_FUNC) &_rsk_raw_to_4byte_signed, 2},
    {"_rsk_raw_to_unsigned_vec", (DL_FUNC) &_rsk_raw_to_unsigned_vec, 4},
    {"_rsk_rsk_find_events", (DL_FUNC) &_rsk_rsk_find_events, 2},
    {"_rsk_rsk_find_times", (DL_FUNC) &_rsk_rsk_find_times, 2},
    {"_rsk_rsk_incomplete_events", (DL_FUNC) &_rsk_rsk_incomplete_events, 3},
    {"_rsk_rsk_raw_to_pressure", (DL_FUNC) &_rsk_rsk_raw_to_pressure, 2},
    {"_rsk_rsk_raw_to_temperature", (DL_FUNC) &_rsk_rsk_raw_to_temperature, 2},
    {"_rsk_rsk_temperature_correction", (DL_FUNC) &_rsk_rsk_temperature_correction, 3},
    {"_rsk_rsk_raw_times", (DL_FUNC) &_rsk_rsk_raw_times, 5},
    {"_rsk_get_header_length", (DL_FUNC) &_rsk_get_header_length, 1},
    {"_rsk_rsk_read_bin", (DL_FUNC) &_rsk_rsk_read_bin, 8},
    {"_rsk_RcppExport_registerCCallable", (DL_FUNC) &_rsk_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_rsk(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
