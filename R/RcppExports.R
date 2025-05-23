# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

cML_estimateC <- function(b_exp, b_out, se_exp, se_out, K, initial_theta, initial_mu, maxit) {
    .Call(`_MRcare_cML_estimateC`, b_exp, b_out, se_exp, se_out, K, initial_theta, initial_mu, maxit)
}

cML_SdThetaC <- function(b_exp, b_out, se_exp, se_out, theta, b_vec, r_vec) {
    .Call(`_MRcare_cML_SdThetaC`, b_exp, b_out, se_exp, se_out, theta, b_vec, r_vec)
}

cML_estimate_randomC <- function(b_exp, b_out, se_exp, se_out, K, random_start, maxit) {
    .Call(`_MRcare_cML_estimate_randomC`, b_exp, b_out, se_exp, se_out, K, random_start, maxit)
}

mr_cMLC <- function(b_exp, b_out, se_exp, se_out, K_vec, random_start, maxit, n) {
    .Call(`_MRcare_mr_cMLC`, b_exp, b_out, se_exp, se_out, K_vec, random_start, maxit, n)
}

mr_cMLC2 <- function(b_exp, b_out, se_exp, se_out, K_vec, random_start, maxit, n, Kall) {
    .Call(`_MRcare_mr_cMLC2`, b_exp, b_out, se_exp, se_out, K_vec, random_start, maxit, n, Kall)
}

care_rep <- function(MRdat, thetaInit, seInit, nrep, random_start, maxit, n) {
    .Call(`_MRcare_care_rep`, MRdat, thetaInit, seInit, nrep, random_start, maxit, n)
}

cML_estimate_randomCtest <- function(b_exp, b_out, se_exp, se_out, K, random_start, maxit) {
    .Call(`_MRcare_cML_estimate_randomCtest`, b_exp, b_out, se_exp, se_out, K, random_start, maxit)
}

cML_SdThetaCtest <- function(b_exp, b_out, se_exp, se_out, theta, b_vec, r_vec) {
    .Call(`_MRcare_cML_SdThetaCtest`, b_exp, b_out, se_exp, se_out, theta, b_vec, r_vec)
}

mr_CD_measureC <- function(b_exp, b_out, se2_exp, se2_out, w) {
    .Call(`_MRcare_mr_CD_measureC`, b_exp, b_out, se2_exp, se2_out, w)
}

cML_CD_estimateC2 <- function(b_exp, b_out, se2_exp, se2_out, w, K, initial_theta, maxit) {
    .Call(`_MRcare_cML_CD_estimateC2`, b_exp, b_out, se2_exp, se2_out, w, K, initial_theta, maxit)
}

cML_CD_randomC <- function(b_exp, b_out, se2_exp, se2_out, w, K, random_start, maxit, initial_theta) {
    .Call(`_MRcare_cML_CD_randomC`, b_exp, b_out, se2_exp, se2_out, w, K, random_start, maxit, initial_theta)
}

mr_cMLC_CD <- function(b_exp, b_out, se2_exp, se2_out, K_vec, w, random_start, maxit, n) {
    .Call(`_MRcare_mr_cMLC_CD`, b_exp, b_out, se2_exp, se2_out, K_vec, w, random_start, maxit, n)
}

mr_cMLC_CD_boot_efron <- function(b_exp, b_out, se2_exp, se2_out, K_vec, wAll, random_start, maxit, n, nrep) {
    .Call(`_MRcare_mr_cMLC_CD_boot_efron`, b_exp, b_out, se2_exp, se2_out, K_vec, wAll, random_start, maxit, n, nrep)
}

mr_cMLC_CD_bootdir_efron <- function(b_exp, b_out, se2_exp, se2_out, K_vec, w, random_start, maxit, n, nrep) {
    .Call(`_MRcare_mr_cMLC_CD_bootdir_efron`, b_exp, b_out, se2_exp, se2_out, K_vec, w, random_start, maxit, n, nrep)
}

mr_cMLC_CD3 <- function(b_exp, b_out, se2_exp, se2_out, K_vec, w, random_start, maxit, n) {
    .Call(`_MRcare_mr_cMLC_CD3`, b_exp, b_out, se2_exp, se2_out, K_vec, w, random_start, maxit, n)
}

