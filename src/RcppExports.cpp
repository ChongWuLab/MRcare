// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cML_estimateC
Rcpp::List cML_estimateC(arma::vec& b_exp, arma::vec& b_out, arma::vec& se_exp, arma::vec& se_out, const int K, double initial_theta, arma::vec initial_mu, const int maxit);
RcppExport SEXP _MRcare_cML_estimateC(SEXP b_expSEXP, SEXP b_outSEXP, SEXP se_expSEXP, SEXP se_outSEXP, SEXP KSEXP, SEXP initial_thetaSEXP, SEXP initial_muSEXP, SEXP maxitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type b_exp(b_expSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type b_out(b_outSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se_exp(se_expSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se_out(se_outSEXP);
    Rcpp::traits::input_parameter< const int >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type initial_theta(initial_thetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type initial_mu(initial_muSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit(maxitSEXP);
    rcpp_result_gen = Rcpp::wrap(cML_estimateC(b_exp, b_out, se_exp, se_out, K, initial_theta, initial_mu, maxit));
    return rcpp_result_gen;
END_RCPP
}
// cML_SdThetaC
double cML_SdThetaC(arma::vec& b_exp, arma::vec& b_out, arma::vec& se_exp, arma::vec& se_out, double theta, arma::vec b_vec, arma::vec r_vec);
RcppExport SEXP _MRcare_cML_SdThetaC(SEXP b_expSEXP, SEXP b_outSEXP, SEXP se_expSEXP, SEXP se_outSEXP, SEXP thetaSEXP, SEXP b_vecSEXP, SEXP r_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type b_exp(b_expSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type b_out(b_outSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se_exp(se_expSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se_out(se_outSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b_vec(b_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type r_vec(r_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(cML_SdThetaC(b_exp, b_out, se_exp, se_out, theta, b_vec, r_vec));
    return rcpp_result_gen;
END_RCPP
}
// cML_estimate_randomC
Rcpp::List cML_estimate_randomC(arma::vec& b_exp, arma::vec& b_out, arma::vec& se_exp, arma::vec& se_out, const int K, const int random_start, const int maxit);
RcppExport SEXP _MRcare_cML_estimate_randomC(SEXP b_expSEXP, SEXP b_outSEXP, SEXP se_expSEXP, SEXP se_outSEXP, SEXP KSEXP, SEXP random_startSEXP, SEXP maxitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type b_exp(b_expSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type b_out(b_outSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se_exp(se_expSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se_out(se_outSEXP);
    Rcpp::traits::input_parameter< const int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const int >::type random_start(random_startSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit(maxitSEXP);
    rcpp_result_gen = Rcpp::wrap(cML_estimate_randomC(b_exp, b_out, se_exp, se_out, K, random_start, maxit));
    return rcpp_result_gen;
END_RCPP
}
// mr_cMLC
Rcpp::List mr_cMLC(arma::vec& b_exp, arma::vec& b_out, arma::vec& se_exp, arma::vec& se_out, arma::vec K_vec, const int random_start, const int maxit, const int n);
RcppExport SEXP _MRcare_mr_cMLC(SEXP b_expSEXP, SEXP b_outSEXP, SEXP se_expSEXP, SEXP se_outSEXP, SEXP K_vecSEXP, SEXP random_startSEXP, SEXP maxitSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type b_exp(b_expSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type b_out(b_outSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se_exp(se_expSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se_out(se_outSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type K_vec(K_vecSEXP);
    Rcpp::traits::input_parameter< const int >::type random_start(random_startSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(mr_cMLC(b_exp, b_out, se_exp, se_out, K_vec, random_start, maxit, n));
    return rcpp_result_gen;
END_RCPP
}
// mr_cMLC2
Rcpp::List mr_cMLC2(arma::vec& b_exp, arma::vec& b_out, arma::vec& se_exp, arma::vec& se_out, arma::vec K_vec, const int random_start, const int maxit, const int n, arma::vec Kall);
RcppExport SEXP _MRcare_mr_cMLC2(SEXP b_expSEXP, SEXP b_outSEXP, SEXP se_expSEXP, SEXP se_outSEXP, SEXP K_vecSEXP, SEXP random_startSEXP, SEXP maxitSEXP, SEXP nSEXP, SEXP KallSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type b_exp(b_expSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type b_out(b_outSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se_exp(se_expSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se_out(se_outSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type K_vec(K_vecSEXP);
    Rcpp::traits::input_parameter< const int >::type random_start(random_startSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Kall(KallSEXP);
    rcpp_result_gen = Rcpp::wrap(mr_cMLC2(b_exp, b_out, se_exp, se_out, K_vec, random_start, maxit, n, Kall));
    return rcpp_result_gen;
END_RCPP
}
// care_rep
Rcpp::List care_rep(Rcpp::List& MRdat, double thetaInit, double seInit, int nrep, int random_start, int maxit, int n);
RcppExport SEXP _MRcare_care_rep(SEXP MRdatSEXP, SEXP thetaInitSEXP, SEXP seInitSEXP, SEXP nrepSEXP, SEXP random_startSEXP, SEXP maxitSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type MRdat(MRdatSEXP);
    Rcpp::traits::input_parameter< double >::type thetaInit(thetaInitSEXP);
    Rcpp::traits::input_parameter< double >::type seInit(seInitSEXP);
    Rcpp::traits::input_parameter< int >::type nrep(nrepSEXP);
    Rcpp::traits::input_parameter< int >::type random_start(random_startSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(care_rep(MRdat, thetaInit, seInit, nrep, random_start, maxit, n));
    return rcpp_result_gen;
END_RCPP
}
// cML_estimate_randomCtest
Rcpp::List cML_estimate_randomCtest(arma::vec& b_exp, arma::vec& b_out, arma::vec& se_exp, arma::vec& se_out, const int K, const int random_start, const int maxit);
RcppExport SEXP _MRcare_cML_estimate_randomCtest(SEXP b_expSEXP, SEXP b_outSEXP, SEXP se_expSEXP, SEXP se_outSEXP, SEXP KSEXP, SEXP random_startSEXP, SEXP maxitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type b_exp(b_expSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type b_out(b_outSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se_exp(se_expSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se_out(se_outSEXP);
    Rcpp::traits::input_parameter< const int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const int >::type random_start(random_startSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit(maxitSEXP);
    rcpp_result_gen = Rcpp::wrap(cML_estimate_randomCtest(b_exp, b_out, se_exp, se_out, K, random_start, maxit));
    return rcpp_result_gen;
END_RCPP
}
// cML_SdThetaCtest
Rcpp::List cML_SdThetaCtest(arma::vec& b_exp, arma::vec& b_out, arma::vec& se_exp, arma::vec& se_out, double theta, arma::vec b_vec, arma::vec r_vec);
RcppExport SEXP _MRcare_cML_SdThetaCtest(SEXP b_expSEXP, SEXP b_outSEXP, SEXP se_expSEXP, SEXP se_outSEXP, SEXP thetaSEXP, SEXP b_vecSEXP, SEXP r_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type b_exp(b_expSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type b_out(b_outSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se_exp(se_expSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se_out(se_outSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b_vec(b_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type r_vec(r_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(cML_SdThetaCtest(b_exp, b_out, se_exp, se_out, theta, b_vec, r_vec));
    return rcpp_result_gen;
END_RCPP
}
// mr_CD_measureC
double mr_CD_measureC(arma::vec b_exp, arma::vec b_out, arma::vec se2_exp, arma::vec se2_out, arma::vec w);
RcppExport SEXP _MRcare_mr_CD_measureC(SEXP b_expSEXP, SEXP b_outSEXP, SEXP se2_expSEXP, SEXP se2_outSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type b_exp(b_expSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b_out(b_outSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type se2_exp(se2_expSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type se2_out(se2_outSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(mr_CD_measureC(b_exp, b_out, se2_exp, se2_out, w));
    return rcpp_result_gen;
END_RCPP
}
// cML_CD_estimateC2
Rcpp::List cML_CD_estimateC2(arma::vec& b_exp, arma::vec& b_out, arma::vec& se2_exp, arma::vec& se2_out, arma::vec w, const int K, double initial_theta, const int maxit);
RcppExport SEXP _MRcare_cML_CD_estimateC2(SEXP b_expSEXP, SEXP b_outSEXP, SEXP se2_expSEXP, SEXP se2_outSEXP, SEXP wSEXP, SEXP KSEXP, SEXP initial_thetaSEXP, SEXP maxitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type b_exp(b_expSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type b_out(b_outSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se2_exp(se2_expSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se2_out(se2_outSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< const int >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type initial_theta(initial_thetaSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit(maxitSEXP);
    rcpp_result_gen = Rcpp::wrap(cML_CD_estimateC2(b_exp, b_out, se2_exp, se2_out, w, K, initial_theta, maxit));
    return rcpp_result_gen;
END_RCPP
}
// cML_CD_randomC
Rcpp::List cML_CD_randomC(arma::vec& b_exp, arma::vec& b_out, arma::vec& se2_exp, arma::vec& se2_out, arma::vec& w, const int K, const int random_start, const int maxit, double initial_theta);
RcppExport SEXP _MRcare_cML_CD_randomC(SEXP b_expSEXP, SEXP b_outSEXP, SEXP se2_expSEXP, SEXP se2_outSEXP, SEXP wSEXP, SEXP KSEXP, SEXP random_startSEXP, SEXP maxitSEXP, SEXP initial_thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type b_exp(b_expSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type b_out(b_outSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se2_exp(se2_expSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se2_out(se2_outSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const int >::type random_start(random_startSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type initial_theta(initial_thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(cML_CD_randomC(b_exp, b_out, se2_exp, se2_out, w, K, random_start, maxit, initial_theta));
    return rcpp_result_gen;
END_RCPP
}
// mr_cMLC_CD
arma::vec mr_cMLC_CD(arma::vec& b_exp, arma::vec& b_out, arma::vec& se2_exp, arma::vec& se2_out, arma::vec K_vec, arma::vec w, const int random_start, const int maxit, const int n);
RcppExport SEXP _MRcare_mr_cMLC_CD(SEXP b_expSEXP, SEXP b_outSEXP, SEXP se2_expSEXP, SEXP se2_outSEXP, SEXP K_vecSEXP, SEXP wSEXP, SEXP random_startSEXP, SEXP maxitSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type b_exp(b_expSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type b_out(b_outSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se2_exp(se2_expSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se2_out(se2_outSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type K_vec(K_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< const int >::type random_start(random_startSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(mr_cMLC_CD(b_exp, b_out, se2_exp, se2_out, K_vec, w, random_start, maxit, n));
    return rcpp_result_gen;
END_RCPP
}
// mr_cMLC_CD_boot_efron
arma::mat mr_cMLC_CD_boot_efron(arma::vec& b_exp, arma::vec& b_out, arma::vec& se2_exp, arma::vec& se2_out, Rcpp::List K_vec, arma::mat wAll, const int random_start, const int maxit, const int n, const int nrep);
RcppExport SEXP _MRcare_mr_cMLC_CD_boot_efron(SEXP b_expSEXP, SEXP b_outSEXP, SEXP se2_expSEXP, SEXP se2_outSEXP, SEXP K_vecSEXP, SEXP wAllSEXP, SEXP random_startSEXP, SEXP maxitSEXP, SEXP nSEXP, SEXP nrepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type b_exp(b_expSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type b_out(b_outSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se2_exp(se2_expSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se2_out(se2_outSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type K_vec(K_vecSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type wAll(wAllSEXP);
    Rcpp::traits::input_parameter< const int >::type random_start(random_startSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type nrep(nrepSEXP);
    rcpp_result_gen = Rcpp::wrap(mr_cMLC_CD_boot_efron(b_exp, b_out, se2_exp, se2_out, K_vec, wAll, random_start, maxit, n, nrep));
    return rcpp_result_gen;
END_RCPP
}
// mr_cMLC_CD_bootdir_efron
arma::mat mr_cMLC_CD_bootdir_efron(arma::mat& b_exp, arma::mat& b_out, arma::vec& se2_exp, arma::vec& se2_out, Rcpp::List K_vec, arma::vec w, const int random_start, const int maxit, const int n, const int nrep);
RcppExport SEXP _MRcare_mr_cMLC_CD_bootdir_efron(SEXP b_expSEXP, SEXP b_outSEXP, SEXP se2_expSEXP, SEXP se2_outSEXP, SEXP K_vecSEXP, SEXP wSEXP, SEXP random_startSEXP, SEXP maxitSEXP, SEXP nSEXP, SEXP nrepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type b_exp(b_expSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type b_out(b_outSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se2_exp(se2_expSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se2_out(se2_outSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type K_vec(K_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< const int >::type random_start(random_startSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type nrep(nrepSEXP);
    rcpp_result_gen = Rcpp::wrap(mr_cMLC_CD_bootdir_efron(b_exp, b_out, se2_exp, se2_out, K_vec, w, random_start, maxit, n, nrep));
    return rcpp_result_gen;
END_RCPP
}
// mr_cMLC_CD3
Rcpp::List mr_cMLC_CD3(arma::vec& b_exp, arma::vec& b_out, arma::vec& se2_exp, arma::vec& se2_out, arma::vec K_vec, arma::vec w, const int random_start, const int maxit, const int n);
RcppExport SEXP _MRcare_mr_cMLC_CD3(SEXP b_expSEXP, SEXP b_outSEXP, SEXP se2_expSEXP, SEXP se2_outSEXP, SEXP K_vecSEXP, SEXP wSEXP, SEXP random_startSEXP, SEXP maxitSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type b_exp(b_expSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type b_out(b_outSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se2_exp(se2_expSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se2_out(se2_outSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type K_vec(K_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< const int >::type random_start(random_startSEXP);
    Rcpp::traits::input_parameter< const int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(mr_cMLC_CD3(b_exp, b_out, se2_exp, se2_out, K_vec, w, random_start, maxit, n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MRcare_cML_estimateC", (DL_FUNC) &_MRcare_cML_estimateC, 8},
    {"_MRcare_cML_SdThetaC", (DL_FUNC) &_MRcare_cML_SdThetaC, 7},
    {"_MRcare_cML_estimate_randomC", (DL_FUNC) &_MRcare_cML_estimate_randomC, 7},
    {"_MRcare_mr_cMLC", (DL_FUNC) &_MRcare_mr_cMLC, 8},
    {"_MRcare_mr_cMLC2", (DL_FUNC) &_MRcare_mr_cMLC2, 9},
    {"_MRcare_care_rep", (DL_FUNC) &_MRcare_care_rep, 7},
    {"_MRcare_cML_estimate_randomCtest", (DL_FUNC) &_MRcare_cML_estimate_randomCtest, 7},
    {"_MRcare_cML_SdThetaCtest", (DL_FUNC) &_MRcare_cML_SdThetaCtest, 7},
    {"_MRcare_mr_CD_measureC", (DL_FUNC) &_MRcare_mr_CD_measureC, 5},
    {"_MRcare_cML_CD_estimateC2", (DL_FUNC) &_MRcare_cML_CD_estimateC2, 8},
    {"_MRcare_cML_CD_randomC", (DL_FUNC) &_MRcare_cML_CD_randomC, 9},
    {"_MRcare_mr_cMLC_CD", (DL_FUNC) &_MRcare_mr_cMLC_CD, 9},
    {"_MRcare_mr_cMLC_CD_boot_efron", (DL_FUNC) &_MRcare_mr_cMLC_CD_boot_efron, 10},
    {"_MRcare_mr_cMLC_CD_bootdir_efron", (DL_FUNC) &_MRcare_mr_cMLC_CD_bootdir_efron, 10},
    {"_MRcare_mr_cMLC_CD3", (DL_FUNC) &_MRcare_mr_cMLC_CD3, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_MRcare(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
