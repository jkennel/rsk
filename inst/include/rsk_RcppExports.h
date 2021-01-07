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

    inline arma::mat anti_subset(arma::vec& x, arma::uvec& idx, arma::uword nc, arma::uword nr) {
        typedef SEXP(*Ptr_anti_subset)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_anti_subset p_anti_subset = NULL;
        if (p_anti_subset == NULL) {
            validateSignature("arma::mat(*anti_subset)(arma::vec&,arma::uvec&,arma::uword,arma::uword)");
            p_anti_subset = (Ptr_anti_subset)R_GetCCallable("rsk", "_rsk_anti_subset");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_anti_subset(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(idx)), Shield<SEXP>(Rcpp::wrap(nc)), Shield<SEXP>(Rcpp::wrap(nr)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline Rcpp::DatetimeVector rbr_times(const arma::vec& ev_tstamp, const arma::uvec& ev_index, double ti) {
        typedef SEXP(*Ptr_rbr_times)(SEXP,SEXP,SEXP);
        static Ptr_rbr_times p_rbr_times = NULL;
        if (p_rbr_times == NULL) {
            validateSignature("Rcpp::DatetimeVector(*rbr_times)(const arma::vec&,const arma::uvec&,double)");
            p_rbr_times = (Ptr_rbr_times)R_GetCCallable("rsk", "_rsk_rbr_times");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_rbr_times(Shield<SEXP>(Rcpp::wrap(ev_tstamp)), Shield<SEXP>(Rcpp::wrap(ev_index)), Shield<SEXP>(Rcpp::wrap(ti)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::DatetimeVector >(rcpp_result_gen);
    }

    inline arma::vec rbr_raw_to_pressure(const arma::vec& x, const arma::vec& calib) {
        typedef SEXP(*Ptr_rbr_raw_to_pressure)(SEXP,SEXP);
        static Ptr_rbr_raw_to_pressure p_rbr_raw_to_pressure = NULL;
        if (p_rbr_raw_to_pressure == NULL) {
            validateSignature("arma::vec(*rbr_raw_to_pressure)(const arma::vec&,const arma::vec&)");
            p_rbr_raw_to_pressure = (Ptr_rbr_raw_to_pressure)R_GetCCallable("rsk", "_rsk_rbr_raw_to_pressure");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_rbr_raw_to_pressure(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(calib)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::vec >(rcpp_result_gen);
    }

    inline arma::vec rbr_raw_to_temperature(const arma::vec& x, const arma::vec& calib) {
        typedef SEXP(*Ptr_rbr_raw_to_temperature)(SEXP,SEXP);
        static Ptr_rbr_raw_to_temperature p_rbr_raw_to_temperature = NULL;
        if (p_rbr_raw_to_temperature == NULL) {
            validateSignature("arma::vec(*rbr_raw_to_temperature)(const arma::vec&,const arma::vec&)");
            p_rbr_raw_to_temperature = (Ptr_rbr_raw_to_temperature)R_GetCCallable("rsk", "_rsk_rbr_raw_to_temperature");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_rbr_raw_to_temperature(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(calib)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::vec >(rcpp_result_gen);
    }

    inline arma::vec rbr_temperature_correction(const arma::vec& pressure, const arma::vec& temperature, const arma::vec& x) {
        typedef SEXP(*Ptr_rbr_temperature_correction)(SEXP,SEXP,SEXP);
        static Ptr_rbr_temperature_correction p_rbr_temperature_correction = NULL;
        if (p_rbr_temperature_correction == NULL) {
            validateSignature("arma::vec(*rbr_temperature_correction)(const arma::vec&,const arma::vec&,const arma::vec&)");
            p_rbr_temperature_correction = (Ptr_rbr_temperature_correction)R_GetCCallable("rsk", "_rsk_rbr_temperature_correction");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_rbr_temperature_correction(Shield<SEXP>(Rcpp::wrap(pressure)), Shield<SEXP>(Rcpp::wrap(temperature)), Shield<SEXP>(Rcpp::wrap(x)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::vec >(rcpp_result_gen);
    }

    inline Rcpp::DataFrame test_all(arma::vec& raw_vector, arma::uvec idx, size_t nc, const arma::vec& ev_tstamp, const arma::uvec& ev_index, double ti, const arma::mat& base_calib, const arma::vec& is_temp, const arma::vec& temp_comp, size_t pres_index, size_t temp_index, Rcpp::CharacterVector names, int subset = 1) {
        typedef SEXP(*Ptr_test_all)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_test_all p_test_all = NULL;
        if (p_test_all == NULL) {
            validateSignature("Rcpp::DataFrame(*test_all)(arma::vec&,arma::uvec,size_t,const arma::vec&,const arma::uvec&,double,const arma::mat&,const arma::vec&,const arma::vec&,size_t,size_t,Rcpp::CharacterVector,int)");
            p_test_all = (Ptr_test_all)R_GetCCallable("rsk", "_rsk_test_all");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_test_all(Shield<SEXP>(Rcpp::wrap(raw_vector)), Shield<SEXP>(Rcpp::wrap(idx)), Shield<SEXP>(Rcpp::wrap(nc)), Shield<SEXP>(Rcpp::wrap(ev_tstamp)), Shield<SEXP>(Rcpp::wrap(ev_index)), Shield<SEXP>(Rcpp::wrap(ti)), Shield<SEXP>(Rcpp::wrap(base_calib)), Shield<SEXP>(Rcpp::wrap(is_temp)), Shield<SEXP>(Rcpp::wrap(temp_comp)), Shield<SEXP>(Rcpp::wrap(pres_index)), Shield<SEXP>(Rcpp::wrap(temp_index)), Shield<SEXP>(Rcpp::wrap(names)), Shield<SEXP>(Rcpp::wrap(subset)));
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
