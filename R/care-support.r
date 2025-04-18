#' Robust Invalid Variance Weighted Estimator
#'
#' Performs the robust invalid variance weighted estimator for Mendelian randomization.
#'
#' @param gamma.exp_sel Vector of genetic associations with the exposure
#' @param gamma.out_sel Vector of genetic associations with the outcome
#' @param se.exp_sel Vector of standard errors for the exposure associations
#' @param se.out_sel Vector of standard errors for the outcome associations
#' @param etamean Mean of the random variable for re-randomization
#' @param sel.pthr P-value threshold for selecting significant SNPs
#' @param overlap.mat Overlap matrix for correlated samples (NULL for independent samples)
#'
#' @return A list containing the causal estimate, standard error, p-value, and number of SNPs
#'
#'  @export

RIVW <- function(gamma.exp_sel, gamma.out_sel, se.exp_sel, se.out_sel, etamean = 0.5, sel.pthr = 5e-5,overlap.mat = NULL) {
    
    if(!is.null(overlap.mat)) {
        c1 = overlap.mat[1,1]
        c2 = overlap.mat[2,2]
        rho = overlap.mat[1,2]
        
        se2.out_sel = c2 * se.out_sel^2
        se2.exp_sel = c1 * se.exp_sel^2
        
        se.out_sel = sqrt(se2.out_sel)
        se.exp_sel = sqrt(se2.exp_sel)
    }

    # Step 2: Select significant SNPs:
    C_sel = qnorm(sel.pthr/2,lower.tail = FALSE) #* sqrt(1 + eta^2)
    
    # Step 3. Construct the unbiased carved estimator (also the UMVUE)
    alpha1 = (-C_sel - gamma.exp_sel/se.exp_sel) / etamean
    alpha2 = (C_sel - gamma.exp_sel/se.exp_sel) / etamean
    gamma.carve = gamma.exp_sel - (se.exp_sel/etamean) * ( (dnorm(alpha2) - dnorm(alpha1)) / (pnorm(alpha1) + 1 - pnorm(alpha2)) )
    sigma2.carve = (1 - ((alpha2*dnorm(alpha2) - alpha1*dnorm(alpha1)) / (1 - pnorm(alpha2) + pnorm(alpha1) ) - ((dnorm(alpha2) - dnorm(alpha1))/(1 - pnorm(alpha2) + pnorm(alpha1)))^2) / etamean^2 ) * se.exp_sel^2
    
    beta = sum(gamma.carve*gamma.out_sel*(1/se.out_sel^2)) / sum((gamma.carve^2 - sigma2.carve)/se.out_sel^2)
    
    # estimation based on regression residuals
    
    RIVW.var = sum( (gamma.out_sel * gamma.carve - beta * (gamma.carve^2 - sigma2.carve) )^2 / se.out_sel^4) / (sum((gamma.carve^2 - sigma2.carve) / se.out_sel^2) )^2
    sd = sqrt(RIVW.var)
    
    p = pnorm(abs(beta/sd), lower.tail = F) * 2
    return(list(beta = beta, sd = sd, p = p,nIV = length(gamma.exp_sel)))
}

#' Correct Winner's Curse Bias in GWAS Summary Statistics
#'
#' Corrects bias in genetic associations due to winner's curse.
#'
#' @param gamma.exp Vector of genetic associations with the exposure
#' @param gamma.out Vector of genetic associations with the outcome
#' @param se.exp Vector of standard errors for the exposure associations
#' @param se.out Vector of standard errors for the outcome associations
#' @param pthr P-value threshold for selecting significant SNPs
#'
#' @return A list with bias-corrected associations
#'
#' @keywords internal
correct_bias <- function(gamma.exp,gamma.out,se.exp,se.out,pthr) {
    flip = which(gamma.exp < 0)
    gamma.exp[flip] = -gamma.exp[flip]
    gamma.out[flip] = -gamma.out[flip]
    
    C_sel = qnorm(pthr/2,lower.tail = FALSE)
    ind_filter = which(abs(gamma.exp / se.exp) >= C_sel)
    numIV_sel = length(ind_filter)
    gamma.exp_sel = gamma.exp[ind_filter]
    gamma.out_sel = gamma.out[ind_filter]
    se.exp_sel = se.exp[ind_filter]
    se.out_sel = se.out[ind_filter]
    gamma.exp.pval_sel = 2 * pnorm(-abs(gamma.exp_sel) / se.exp_sel)
    
    # gammahat correction process
    gamma.exp_selMLE  = NULL
    for(j in 1:length(gamma.exp_sel)) {
        rootfinding <- function(x) {
            beta = x
            std = se.exp_sel[j]
            alpha2 = beta/std - C_sel
            alpha1 = -beta/std - C_sel
            beta + std * ((dnorm(alpha2)-dnorm(alpha1)) / (pnorm(alpha1)+pnorm(alpha2))) - gamma.exp_sel[j]
        }
        gamma.exp_selMLE = c(gamma.exp_selMLE, nleqslv(gamma.exp_sel[j], rootfinding)$x) #gamma check
    }
    gamma.exp_sel = gamma.exp_selMLE #gamma.exp_sel - biasterm #gamma correct
    
    out = list(gamma.exp = gamma.exp_sel,gamma.out = gamma.out_sel, se.exp = se.exp_sel, se.out = se.out_sel,ind_filter = ind_filter)
    
    return(out)
}




#' Re-randomization Bias Correction
#'
#' Corrects bias using a re-randomization approach.
#'
#' @param gamma.exp Vector of genetic associations with the exposure
#' @param gamma.out Vector of genetic associations with the outcome
#' @param se.exp Vector of standard errors for the exposure associations
#' @param se.out Vector of standard errors for the outcome associations
#' @param etamean Mean of the random variable for re-randomization
#' @param pthr P-value threshold for selecting significant SNPs
#' @param var.type Variance type ("upper" or "org")
#'
#' @return A list with bias-corrected associations
#'
#' @keywords internal
rerandomization_bias <- function(gamma.exp,gamma.out,se.exp,se.out,etamean = 0.5,pthr = 5e-5,var.type = "upper") {
    

    gamma.exp_sel = gamma.exp
    gamma.out_sel = gamma.out
    se.exp_sel = se.exp
    se.out_sel = se.out

    C_sel = qnorm(pthr/2,lower.tail = FALSE)

    # Step 2. Construct the unbiased carved estimator (also the UMVUE)
    alpha1 = (-C_sel - gamma.exp_sel/se.exp_sel) / etamean
    alpha2 = (C_sel - gamma.exp_sel/se.exp_sel) / etamean
    gamma.carve = gamma.exp_sel - (se.exp_sel/etamean) * ( (dnorm(alpha2) - dnorm(alpha1)) / (pnorm(alpha1) + 1 - pnorm(alpha2)) )
    
    
    if(var.type == "upper") {
        sigma2.carve = (1 + etamean^2 ) * se.exp_sel^2
        se.exp_sel = sqrt(sigma2.carve)
    }
    
    gamma.exp_sel = gamma.carve
    
    out = list(gamma.exp = gamma.exp_sel,gamma.out = gamma.out_sel, se.exp = se.exp_sel, se.out = se.out_sel)
    
    return(out)
}


#' CARE Bootstrap Analysis
#'
#' Performs CARE analysis with bootstrapping for robust estimation.
#'
#' @param gamma.exp_sel Vector of genetic associations with the exposure
#' @param gamma.out_sel Vector of genetic associations with the outcome
#' @param se.exp_sel Vector of standard errors for the exposure associations
#' @param se.out_sel Vector of standard errors for the outcome associations
#' @param nx Sample size for exposure
#' @param ny Sample size for outcome
#' @param nrep Number of bootstrap repetitions
#' @param algorithm Algorithm to use: "Lasso" or "CD" (coordinate descent)
#' @param random_start Number of random starting points
#' @param biascorrect Bias correction method: "no", "direct", "rerand", "rerand2", or "rerand3"
#' @param etamean Mean of the random variable for re-randomization
#' @param pthr P-value threshold for selecting significant SNPs
#'
#' @return A list containing CARE analysis results
#'
#' @export
#'
#' @examples
#' \dontrun{
#' results <- CARE2_boot(beta_exp, beta_out, se_exp, se_out, 10000, 10000)
#' }
CARE2_boot <- function(gamma.exp_sel, gamma.out_sel, se.exp_sel, se.out_sel, nx,ny,nrep = 1000, algorithm = "Lasso", random_start = 0, biascorrect = "no",etamean = 0.5, pthr = 5e-5){
    
  # select IV and correct winner curse bias
  if(biascorrect == "no") {
    
  } else if (biascorrect == "direct") {
    dat = correct_bias(gamma.exp_sel,gamma.out_sel,se.exp_sel,se.out_sel,pthr)
    
    gamma.exp_sel = dat$gamma.exp
    gamma.out_sel = dat$gamma.out
    se.exp_sel = dat$se.exp
    se.out_sel = dat$se.out
    
  } else if (biascorrect == "rerand") {
    dat = rerandomization_bias(gamma.exp_sel,gamma.out_sel,se.exp_sel,se.out_sel,etamean,pthr,"upper")
    
    gamma.exp_sel = dat$gamma.exp
    gamma.out_sel = dat$gamma.out
    se.exp_sel = dat$se.exp
    se.out_sel = dat$se.out
  } else if (biascorrect == "rerand2") { #we use the original se for exposure, mimic the situation for direct correction
    # TODO: we may delete these two..
    dat = rerandomization_bias(gamma.exp_sel,gamma.out_sel,se.exp_sel,se.out_sel,etamean,pthr,"org")
    
    gamma.exp_sel = dat$gamma.exp
    gamma.out_sel = dat$gamma.out
    se.exp_sel = dat$se.exp
    se.out_sel = dat$se.out
  } else if (biascorrect == "rerand3") { #we remove IVs that are weak after winner's curse correction
    dat = rerandomization_bias(gamma.exp_sel,gamma.out_sel,se.exp_sel,se.out_sel,etamean,pthr,"org")
    
    gamma.exp_sel = dat$gamma.exp
    
    gamma.out_sel = dat$gamma.out
    se.exp_sel = dat$se.exp
    se.out_sel = dat$se.out
    
    sel.indx = abs(gamma.exp_sel) <= 3 * se.exp_sel
    
    gamma.exp_sel = gamma.exp_sel[sel.indx]
    gamma.out_sel = gamma.out_sel[sel.indx]
    se.exp_sel = se.exp_sel[sel.indx]
    se.out_sel = se.out_sel[sel.indx]
  }
    

    nIV = length(gamma.exp_sel)

    # Calculate F-statistic
    F = sum(gamma.exp_sel^2 / se.exp_sel^2 )/length(gamma.exp_sel) - 1
    
    
    sim.setting = c(nIV, F) # the statistics for this simulation setting
    names(sim.setting) = c("nIV","F")
    
    se2.out_sel = se.out_sel^2
    se2.exp_sel = se.exp_sel^2
    
    se.out_sel = sqrt(se2.out_sel)
    se.exp_sepl = sqrt(se2.exp_sel)

    ##############################################################
    ### Bootstrap to consider post-selection bias
    ##############################################################
    wAll = rmultinom(nrep, size = nIV, prob = rep(1/nIV,nIV))
    
    colSums(wAll!=0)
    if(algorithm == "Lasso") {
      #lambda.max   <- 600#max(abs(gamma.out_sel/se2.out_sel))
      #lambda.min   <- 250#lambda.max * 1E-5
      #Lambda <- exp(1) ^ seq(log(lambda.max), log(lambda.min), length = 100)
      w = rep(1,length(gamma.exp_sel))
      
      # Find optimal lambda values
      Lambda = findLambda2(b_exp = gamma.exp_sel, b_out = gamma.out_sel, se2_exp = se2.exp_sel, se2_out = se2.out_sel, w = w,  initial_theta = 0.0, maxit = 1000,nRefine = 10)
      
      MRcML_result = mr_cMLC_dIVW_lasso_efron(b_exp = gamma.exp_sel, b_out = gamma.out_sel, se2_exp = se2.exp_sel, se2_out = se2.out_sel, Lambda = Lambda,wAll = wAll, initial_theta = 0.0, maxit = 1000,  n = min(nx,ny), nrep = nrep)
    } else {
      Kvec = colSums(wAll!=0)
      KvecList = list()
      for(j in 1:dim(wAll)[2]) {
        tmp = floor(seq(0,Kvec[j] - 2,length.out = 50))
        tmp = unique(tmp)
        KvecList[[j]] = tmp
      }
      
      MRcML_result = mr_cMLC_CD_boot_efron(b_exp = gamma.exp_sel, b_out = gamma.out_sel, se2_exp = se2.exp_sel, se2_out = se2.out_sel, K_vec = KvecList,wAll = wAll, random_start = random_start, maxit = 1000,  n = min(nx,ny), nrep = nrep)
    }
   
    wALLmean = rowSums(wAll) / nrep
    
    wAll2 = wAll - matrix(rep(wALLmean,nrep),dim(wAll)[1],dim(wAll)[2],byrow = FALSE)
    
    t = MRcML_result[,1]
    t = t - mean(t)
    
    wAll3 = wAll2 %*% t /nrep
    MA_se_efron = sqrt(sum(wAll3^2))
    
    
    t = MRcML_result[,2]
    t = t - mean(t)
    
    wAll3 = wAll2 %*% t /nrep
    BIC_se_efron = sqrt(sum(wAll3^2))
    
    
    t = MRcML_result[,3]
    t = t - mean(t)
    
    wAll3 = wAll2 %*% t /nrep
    MA_AIC_se_efron = sqrt(sum(wAll3^2))
    
    t = MRcML_result[,4]
    t = t - mean(t)
    
    wAll3 = wAll2 %*% t /nrep
    AIC_se_efron = sqrt(sum(wAll3^2))
    
    
    MA_theta = mean(MRcML_result[,1])
    MA_se = sqrt(var(MRcML_result[,1]))
    MA_p = pnorm(abs(MA_theta/MA_se), lower.tail = FALSE) * 2
    
    BIC_theta = mean(MRcML_result[,2])
    BIC_se = sqrt(var(MRcML_result[,2]))
    BIC_p = pnorm(abs(BIC_theta/BIC_se), lower.tail = FALSE) * 2
    
    MA_AIC_theta = mean(MRcML_result[,3])
    MA_AIC_se = sqrt(var(MRcML_result[,3]))
    MA_AIC_p = pnorm(abs(MA_AIC_theta/MA_AIC_se), lower.tail = FALSE) * 2
    
    AIC_theta = mean(MRcML_result[,4])
    AIC_se = sqrt(var(MRcML_result[,4]))
    AIC_p = pnorm(abs(AIC_theta/AIC_se), lower.tail = FALSE) * 2
    
    
    w = rep(1,length(gamma.exp_sel))
    se2.exp_sel = se.exp_sel^2
    se2.out_sel = se.out_sel^2
    
    if(algorithm == "Lasso") {
      mrCMLraps = cML_dIVW_LassoTest(b_exp = gamma.exp_sel, b_out = gamma.out_sel, se2_exp = se2.exp_sel, se2_out = se2.out_sel, w = w, Lambda = Lambda, initial_theta = 0.0, maxit = 1000,  n = min(nx,ny))
    } else {
      tmp = floor(seq(0,length(gamma.exp_sel) - 2,length.out = 50))
      K_vec = unique(tmp)
      
      mrCMLraps = mr_cMLC_CD3(b_exp = gamma.exp_sel, b_out = gamma.out_sel, se2_exp = se2.exp_sel, se2_out = se2.out_sel, K_vec = K_vec, w= w, random_start = random_start, maxit = 1000,  n = min(nx,ny))
    }
    
  
    theta2 = mrCMLraps[[1]]
    
    MA_p2 = pnorm(abs(theta2[1]/MA_se), lower.tail = FALSE) * 2
    BIC_p2 = pnorm(abs(theta2[2]/BIC_se), lower.tail = FALSE) * 2
    
    MA_AIC_p2 = pnorm(abs(theta2[3]/MA_AIC_se), lower.tail = FALSE) * 2
    AIC_p2 = pnorm(abs(theta2[4]/AIC_se), lower.tail = FALSE) * 2
    
    
    MA_p3 = pnorm(abs(MA_theta/MA_se_efron), lower.tail = FALSE) * 2
    BIC_p3 = pnorm(abs(BIC_theta/BIC_se_efron), lower.tail = FALSE) * 2
    
    MA_AIC_p3 = pnorm(abs(MA_AIC_theta/MA_AIC_se_efron), lower.tail = FALSE) * 2
    AIC_p3 = pnorm(abs(AIC_theta/AIC_se_efron), lower.tail = FALSE) * 2
    
    
    MA_BIC = c(MA_theta,MA_se,MA_p,theta2[1],MA_p2,MA_se_efron,MA_p3)
    names(MA_BIC) = c("tilde_theta","se","p","hat_theta","hat_p","Efron_se","Efron_p")
    
    BIC = c(BIC_theta,BIC_se,BIC_p,theta2[2],BIC_p2,BIC_se_efron,BIC_p3)
    names(BIC) = c("tilde_theta","se","p","hat_theta","hat_p","Efron_se","Efron_p")
    
    MA_AIC = c(MA_AIC_theta,MA_AIC_se,MA_AIC_p,theta2[3],MA_AIC_p2,MA_AIC_se_efron,MA_AIC_p3)
    names(MA_AIC) = c("tilde_theta","se","p","hat_theta","hat_p","Efron_se","Efron_p")
    
    AIC = c(AIC_theta,AIC_se,AIC_p,theta2[4],AIC_p2,AIC_se_efron,AIC_p3)
    names(AIC) = c("tilde_theta","se","p","hat_theta","hat_p","Efron_se","Efron_p")
    
    ## Number of invalid IVs
    BIC_IV = mrCMLraps[[2]]
    BIC_IV = sum(BIC_IV==0)
    
    AIC_IV = mrCMLraps[[3]]
    AIC_IV = sum(BIC_IV==0)
    
    IV = c(BIC_IV,AIC_IV)
    names(IV) = c("BIC","AIC")
    
    
    res = list(MA_BIC=MA_BIC, BIC=BIC, MA_AIC = MA_AIC, AIC = AIC, setting = sim.setting, IV = IV)
    
    out = list(res = res, BIC_IV = mrCMLraps[[2]], AIC_IV = mrCMLraps[[3]], MRcML_boot = MRcML_result)
    return(out)
}




