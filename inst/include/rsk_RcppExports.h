// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_rsk_RCPPEXPORTS_H_GEN_
#define RCPP_rsk_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace rsk {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("rsk", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("rsk", "_rsk_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in rsk");
            }
        }
    }

    inline uint32_t raw_to_4byte_unsigned(const Rcpp::RawVector x) {
        typedef SEXP(*Ptr_raw_to_4byte_unsigned)(SEXP);
        static Ptr_raw_to_4byte_unsigned p_raw_to_4byte_unsigned = NULL;
        if (p_raw_to_4byte_unsigned == NULL) {
            validateSignature("uint32_t(*raw_to_4byte_unsigned)(const Rcpp::RawVector)");
            p_raw_to_4byte_unsigned = (Ptr_raw_to_4byte_unsigned)R_GetCCallable("rsk", "_rsk_raw_to_4byte_unsigned");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_raw_to_4byte_unsigned(Shield<SEXP>(Rcpp::wrap(x)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<uint32_t >(rcpp_result_gen);
    }

    inline int32_t raw_to_4byte_signed(const Rcpp::RawVector x, size_t index) {
        typedef SEXP(*Ptr_raw_to_4byte_signed)(SEXP,SEXP);
        static Ptr_raw_to_4byte_signed p_raw_to_4byte_signed = NULL;
        if (p_raw_to_4byte_signed == NULL) {
            validateSignature("int32_t(*raw_to_4byte_signed)(const Rcpp::RawVector,size_t)");
            p_raw_to_4byte_signed = (Ptr_raw_to_4byte_signed)R_GetCCallable("rsk", "_rsk_raw_to_4byte_signed");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_raw_to_4byte_signed(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(index)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<int32_t >(rcpp_result_gen);
    }

    inline Rcpp::NumericMatrix raw_to_unsigned_vec(const Rcpp::RawVector x, const size_t header_length, const size_t n_columns, const Rcpp::IntegerVector to_remove) {
        typedef SEXP(*Ptr_raw_to_unsigned_vec)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_raw_to_unsigned_vec p_raw_to_unsigned_vec = NULL;
        if (p_raw_to_unsigned_vec == NULL) {
            validateSignature("Rcpp::NumericMatrix(*raw_to_unsigned_vec)(const Rcpp::RawVector,const size_t,const size_t,const Rcpp::IntegerVector)");
            p_raw_to_unsigned_vec = (Ptr_raw_to_unsigned_vec)R_GetCCallable("rsk", "_rsk_raw_to_unsigned_vec");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_raw_to_unsigned_vec(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(header_length)), Shield<SEXP>(Rcpp::wrap(n_columns)), Shield<SEXP>(Rcpp::wrap(to_remove)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::NumericMatrix >(rcpp_result_gen);
    }

    inline Rcpp::List rsk_find_events(const Rcpp::RawVector x, const size_t header_length) {
        typedef SEXP(*Ptr_rsk_find_events)(SEXP,SEXP);
        static Ptr_rsk_find_events p_rsk_find_events = NULL;
        if (p_rsk_find_events == NULL) {
            validateSignature("Rcpp::List(*rsk_find_events)(const Rcpp::RawVector,const size_t)");
            p_rsk_find_events = (Ptr_rsk_find_events)R_GetCCallable("rsk", "_rsk_rsk_find_events");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_rsk_find_events(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(header_length)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

    inline Rcpp::IntegerVector rsk_incomplete_events(const Rcpp::IntegerVector x, const Rcpp::IntegerVector times, const size_t n_cols, const bool f5) {
        typedef SEXP(*Ptr_rsk_incomplete_events)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_rsk_incomplete_events p_rsk_incomplete_events = NULL;
        if (p_rsk_incomplete_events == NULL) {
            validateSignature("Rcpp::IntegerVector(*rsk_incomplete_events)(const Rcpp::IntegerVector,const Rcpp::IntegerVector,const size_t,const bool)");
            p_rsk_incomplete_events = (Ptr_rsk_incomplete_events)R_GetCCallable("rsk", "_rsk_rsk_incomplete_events");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_rsk_incomplete_events(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(times)), Shield<SEXP>(Rcpp::wrap(n_cols)), Shield<SEXP>(Rcpp::wrap(f5)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::IntegerVector >(rcpp_result_gen);
    }

    inline arma::vec rsk_raw_to_pressure(const arma::vec x, const arma::vec calib) {
        typedef SEXP(*Ptr_rsk_raw_to_pressure)(SEXP,SEXP);
        static Ptr_rsk_raw_to_pressure p_rsk_raw_to_pressure = NULL;
        if (p_rsk_raw_to_pressure == NULL) {
            validateSignature("arma::vec(*rsk_raw_to_pressure)(const arma::vec,const arma::vec)");
            p_rsk_raw_to_pressure = (Ptr_rsk_raw_to_pressure)R_GetCCallable("rsk", "_rsk_rsk_raw_to_pressure");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_rsk_raw_to_pressure(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(calib)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::vec >(rcpp_result_gen);
    }

    inline arma::vec rsk_raw_to_temperature(const arma::vec x, const arma::vec calib) {
        typedef SEXP(*Ptr_rsk_raw_to_temperature)(SEXP,SEXP);
        static Ptr_rsk_raw_to_temperature p_rsk_raw_to_temperature = NULL;
        if (p_rsk_raw_to_temperature == NULL) {
            validateSignature("arma::vec(*rsk_raw_to_temperature)(const arma::vec,const arma::vec)");
            p_rsk_raw_to_temperature = (Ptr_rsk_raw_to_temperature)R_GetCCallable("rsk", "_rsk_rsk_raw_to_temperature");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_rsk_raw_to_temperature(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(calib)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::vec >(rcpp_result_gen);
    }

    inline arma::vec rsk_temperature_correction(arma::vec out, arma::vec tcal, const arma::vec x) {
        typedef SEXP(*Ptr_rsk_temperature_correction)(SEXP,SEXP,SEXP);
        static Ptr_rsk_temperature_correction p_rsk_temperature_correction = NULL;
        if (p_rsk_temperature_correction == NULL) {
            validateSignature("arma::vec(*rsk_temperature_correction)(arma::vec,arma::vec,const arma::vec)");
            p_rsk_temperature_correction = (Ptr_rsk_temperature_correction)R_GetCCallable("rsk", "_rsk_rsk_temperature_correction");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_rsk_temperature_correction(Shield<SEXP>(Rcpp::wrap(out)), Shield<SEXP>(Rcpp::wrap(tcal)), Shield<SEXP>(Rcpp::wrap(x)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::vec >(rcpp_result_gen);
    }

    inline Rcpp::DatetimeVector rsk_raw_times(const Rcpp::IntegerVector raw_tstamp, const Rcpp::IntegerVector raw_index, const size_t n_times, const size_t n_columns, double measurement_interval) {
        typedef SEXP(*Ptr_rsk_raw_times)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_rsk_raw_times p_rsk_raw_times = NULL;
        if (p_rsk_raw_times == NULL) {
            validateSignature("Rcpp::DatetimeVector(*rsk_raw_times)(const Rcpp::IntegerVector,const Rcpp::IntegerVector,const size_t,const size_t,double)");
            p_rsk_raw_times = (Ptr_rsk_raw_times)R_GetCCallable("rsk", "_rsk_rsk_raw_times");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_rsk_raw_times(Shield<SEXP>(Rcpp::wrap(raw_tstamp)), Shield<SEXP>(Rcpp::wrap(raw_index)), Shield<SEXP>(Rcpp::wrap(n_times)), Shield<SEXP>(Rcpp::wrap(n_columns)), Shield<SEXP>(Rcpp::wrap(measurement_interval)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::DatetimeVector >(rcpp_result_gen);
    }

    inline size_t get_header_length(Rcpp::RawVector x) {
        typedef SEXP(*Ptr_get_header_length)(SEXP);
        static Ptr_get_header_length p_get_header_length = NULL;
        if (p_get_header_length == NULL) {
            validateSignature("size_t(*get_header_length)(Rcpp::RawVector)");
            p_get_header_length = (Ptr_get_header_length)R_GetCCallable("rsk", "_rsk_get_header_length");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_get_header_length(Shield<SEXP>(Rcpp::wrap(x)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<size_t >(rcpp_result_gen);
    }

    inline Rcpp::DataFrame rsk_read_bin(Rcpp::RawVector x, Rcpp::LogicalVector is_temp, Rcpp::NumericMatrix base_calib, Rcpp::NumericVector temp_comp, double measurement_interval, size_t pressure_index, size_t temperature_index, bool keep_raw) {
        typedef SEXP(*Ptr_rsk_read_bin)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_rsk_read_bin p_rsk_read_bin = NULL;
        if (p_rsk_read_bin == NULL) {
            validateSignature("Rcpp::DataFrame(*rsk_read_bin)(Rcpp::RawVector,Rcpp::LogicalVector,Rcpp::NumericMatrix,Rcpp::NumericVector,double,size_t,size_t,bool)");
            p_rsk_read_bin = (Ptr_rsk_read_bin)R_GetCCallable("rsk", "_rsk_rsk_read_bin");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_rsk_read_bin(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(is_temp)), Shield<SEXP>(Rcpp::wrap(base_calib)), Shield<SEXP>(Rcpp::wrap(temp_comp)), Shield<SEXP>(Rcpp::wrap(measurement_interval)), Shield<SEXP>(Rcpp::wrap(pressure_index)), Shield<SEXP>(Rcpp::wrap(temperature_index)), Shield<SEXP>(Rcpp::wrap(keep_raw)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::DataFrame >(rcpp_result_gen);
    }

}

#endif // RCPP_rsk_RCPPEXPORTS_H_GEN_
