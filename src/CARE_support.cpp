// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <Rcpp.h>
#include <vector>

using namespace Rcpp;
using namespace arma;




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List  cML_estimateC(arma::vec &b_exp, arma::vec &b_out, arma::vec &se_exp,  arma::vec & se_out, const int K,  double initial_theta, arma::vec initial_mu, const int maxit) {
    
    const int p = b_exp.size();

    double theta = initial_theta;
    double thetaold = theta - 1;
    
    arma::vec muvec;
    muvec.zeros(p);
    
    arma::vec vimportance;
    vimportance.zeros(p);
    
    arma::vec vbg;
    vbg.zeros(p);
    
    muvec = initial_mu;
    
    int iteindx = 0;
    
    double dlx = 1.0;
    while (dlx > 1e-7 & iteindx < maxit) {
        thetaold = theta;
        iteindx = iteindx + 1;
        
        if(K > 0) {
            
            vimportance = arma::pow(b_out - theta * muvec,2) / arma::pow(se_out,2);
            arma::uvec tmpindx = arma::sort_index(vimportance,"descend");
            arma::uvec tmpindx2 = tmpindx.tail(p-K);
            vbg = b_out - theta * muvec;
            
            vbg(tmpindx2).zeros();
        }
        
        // update mu
        muvec = (b_exp / arma::pow(se_exp,2) + theta * (b_out - vbg) / arma::pow(se_out,2) ) / (1 / arma::pow(se_exp,2) + theta * theta / arma::pow(se_out,2));

        // third, update theta
        theta = arma::accu( ((b_out - vbg) % muvec) / arma::pow(se_out,2))  / arma::accu(arma::pow(muvec,2) / arma::pow(se_out,2));
        
        dlx = thetaold - theta;
        if(dlx <0) {
            dlx = -1 * dlx;
        }
    }
    
    // update vbg and muvec once we find optimal theta
    if(K > 0) {
        arma::uvec nonzeroind = arma::find(vbg != 0);
        muvec.elem(nonzeroind) = b_exp.elem(nonzeroind);
        
        arma::vec tmpvbg = b_out - theta * muvec;
        vbg.elem(nonzeroind) = tmpvbg.elem(nonzeroind);
    }
    
    
    
    Rcpp::List out;
    out["theta"] = theta;
    out["bvec"] = muvec;
    out["rvec"] = vbg;
    out["iteindx"] = iteindx;

    return (out);

}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double  cML_SdThetaC(arma::vec &b_exp, arma::vec &b_out, arma::vec &se_exp,  arma::vec & se_out, double theta, arma::vec b_vec, arma::vec r_vec) {
    
    arma::uvec zeroind = arma::find(r_vec == 0.0);
    
    arma::vec b_vec2 = b_vec.elem(zeroind);
    arma::vec se_out2 = se_out.elem(zeroind);
    arma::vec b_out2 = b_out.elem(zeroind);
    arma::vec se_exp2 = se_exp.elem(zeroind);
    
    double varTheta;
    
    varTheta = 1 / ( arma::accu( arma::pow(b_vec2,2) / arma::pow(se_out2,2))  - arma::accu( arma::pow(2 * theta * b_vec2 - b_out2,2) / ( arma::pow(se_out2,4 ) % (1 / arma::pow(se_exp2,2) + theta * theta / arma::pow(se_out2,2)) ) ) )  ;

    
    if(varTheta <=0) {
        varTheta = -1;
    } else {
        varTheta = sqrt(varTheta);
    }
    
    return (varTheta);

}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List  cML_estimate_randomC(arma::vec &b_exp, arma::vec &b_out, arma::vec &se_exp,  arma::vec & se_out, const int K,  const int random_start, const int maxit) {
    
    const int p = b_exp.size();

    arma:: vec tmp = b_out / b_exp;
    double minTheta = tmp.min();
    double maxTheta = tmp.max();
    
    double rangeTheta = maxTheta - minTheta;
    double initial_theta = 0.0;
    
    arma::vec initial_mu;
    initial_mu.zeros(p);
    
    arma::vec invalidIV;
    invalidIV.zeros(p);
    
    double NegL;
    double theta;
    double sdTheta;
    
    int outj = 0;
    for (int j = 0; j <= random_start; j++) {
        
        if (j > 0) {
            arma::vec tmpInitial_theta = rangeTheta * arma::randu(1) + minTheta;
            initial_theta = tmpInitial_theta(0);
            initial_mu = arma::randn(p) % se_exp + b_exp;
        }
        
        Rcpp::List MLEresult = cML_estimateC(b_exp, b_out, se_exp, se_out, K, initial_theta, initial_mu, maxit);
        
        double thetatmp = MLEresult[0];
        arma::vec b_vec = MLEresult[1];
        arma::vec r_vec = MLEresult[2];
        
        double NegLtmp = arma::accu( arma::pow(b_exp - b_vec,2) / (2 * arma::pow(se_exp,2) ) ) + arma::accu(arma::pow(b_out - thetatmp * b_vec - r_vec,2) / (2 * arma::pow(se_out,2) ) );
        
        double sdThetaTmp = cML_SdThetaC(b_exp, b_out, se_exp,  se_out, thetatmp, b_vec, r_vec);
        
        
        if(sdThetaTmp == -1) {
            continue;
        }
        
        if(outj == 0) {
            NegL = NegLtmp;
            theta = thetatmp;
            invalidIV = r_vec;
            sdTheta = sdThetaTmp;
        }
        
        if(NegLtmp < NegL) {
            NegL = NegLtmp;
            theta = thetatmp;
            invalidIV = r_vec;
            sdTheta = sdThetaTmp;
        }
        
        outj = outj + 1;
        
    }
    
    Rcpp::List out;
    out["theta"] = theta;
    out["se"] = sdTheta;
    out["l"] = NegL;
    out["r_est"] = invalidIV;
    out["outj"] = outj;

    return (out);
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List  mr_cMLC(arma::vec &b_exp, arma::vec &b_out, arma::vec &se_exp,  arma::vec & se_out, arma::vec K_vec,  const int random_start, const int maxit, const int n) {
    
    const int M = K_vec.size();
    const int p = b_exp.size();

    arma::vec theta_v;
    arma::vec sd_v;
    arma::vec l_v;
    arma::mat invalid_mat(p,M);
    
    theta_v.zeros(M);
    sd_v.zeros(M);
    l_v.zeros(M);
    invalid_mat.fill(0);
    
    for (int j = 0; j < M; j++) {
        int K = K_vec(j);
        Rcpp::List randRes = cML_estimate_randomC(b_exp, b_out, se_exp, se_out, K, random_start, maxit);
        
        theta_v(j) = randRes[0];
        sd_v(j) = randRes[1];
        l_v(j) = randRes[2];
        
        arma::vec tmpInvalid = randRes[3];
        invalid_mat.col(j) = tmpInvalid;
    }
    
    /*
    arma::vec BIC_v = arma::log(n) * K_vec + 2 * l_v;
    BIC_v = BIC_v - BIC_v.min();
    arma::vec weight_v = arma::exp(-0.5 * BIC_v);
    weight_v = weight_v / arma::accu(weight_v);
    
    double MA_BIC_Theta = arma::accu(theta_v * weight_v);
    
*/
    
    Rcpp::List out;
    out["theta"] = theta_v;
    out["se"] = sd_v;
    out["NegL"] = l_v;
    out["r_est"] = invalid_mat;

    return (out);
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List  mr_cMLC2(arma::vec &b_exp, arma::vec &b_out, arma::vec &se_exp,  arma::vec & se_out, arma::vec K_vec,  const int random_start, const int maxit, const int n, arma::vec Kall) {
    
    const int M = K_vec.size();
    const int p = b_exp.size();

    arma::vec theta_v;
    arma::vec sd_v;
    arma::vec l_v;
    arma::mat invalid_mat(p,M);
    
    theta_v.zeros(M);
    sd_v.zeros(M);
    l_v.zeros(M);
    invalid_mat.fill(0);
    
    for (int j = 0; j < M; j++) {
        int K = K_vec(j);
        Rcpp::List randRes = cML_estimate_randomC(b_exp, b_out, se_exp, se_out, K, random_start, maxit);

        theta_v(j) = randRes[0];
        sd_v(j) = randRes[1];
        l_v(j) = randRes[2];
        
        arma::vec tmpInvalid = randRes[3];
        invalid_mat.col(j) = tmpInvalid;
    }
    
    arma::uvec validind = arma::find(sd_v != -1);
    
    arma::vec theta_v2 = theta_v.elem(validind);
    arma::vec sd_v2 = sd_v.elem(validind);
    arma::vec l_v2 = l_v.elem(validind);
    
    arma::vec K_vec2 = K_vec.elem(validind);
    arma::vec K_ind = Kall.elem(validind);
    
    
    arma::vec BIC_v = log(n) * K_vec2 + 2 * l_v2;
    BIC_v = BIC_v - BIC_v.min();
    arma::vec weight_v = arma::exp(-0.5 * BIC_v);
    weight_v = weight_v / arma::accu(weight_v);
    
    double MA_BIC_Theta = arma::accu(theta_v2 % weight_v);
    double MA_BIC_se = arma::accu(weight_v % sqrt(arma::pow(sd_v2,2) + arma::pow(theta_v2 - MA_BIC_Theta,2) ) );
    
    arma::uword i = BIC_v.index_min();
    
    double BIC_Theta = theta_v2(i);
    double BIC_se = sd_v2(i);
    
    arma::uword i2 = K_ind(i);
    
    arma::vec invalidIVout = invalid_mat.col(i2);
    
    BIC_v = log(n) * K_vec2 + 2 * l_v2;
    
    Rcpp::List out;
    out["MA_BIC_theta"] = MA_BIC_Theta;
    out["MA_BIC_se"] = MA_BIC_se;
    out["BIC_theta"] = BIC_Theta;
    out["BIC_se"] = BIC_se;
    out["BIC_invalid"] = invalidIVout;
    out["BIC_vec"] = BIC_v;

    return (out);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List  care_rep(Rcpp::List &MRdat, double thetaInit, double seInit, int nrep,int random_start,int maxit, int n) {// This function contains bugs, we use the R version first to see if the idea are correct...
    
    double theta2 = thetaInit - 1;
    double se2 = seInit;
    
    arma::mat FinalOut(nrep,2);
    FinalOut.fill(0);
    
    int outcount = 0;
    //arma::mat InvalidMat(p,nrep);
    //InvalidMat.fill(0);
    
    while( theta2 - thetaInit > 0.01 || theta2 - thetaInit < -0.01 & outcount < 5) {
        
        for (int i = 0; i < nrep; i++) {
            arma::mat MRdattmp = MRdat[i];
            
            arma::vec gamma_exp = MRdattmp.col(0);
            arma::vec se_exp = MRdattmp.col(1);
            arma::vec gamma_out = MRdattmp.col(2);
            arma::vec se_out = MRdattmp.col(3);
            
            Rcout << "read data" << i << std::endl;

            arma::vec Tpleio = (gamma_out - theta2 * gamma_exp) / sqrt(arma::pow(se_out,2) + theta2 * theta2 * arma::pow(se_exp,2) + se2 * se2 * arma::pow(gamma_exp,2) );
            
            arma::uvec negind = arma::find(Tpleio <0);
            Tpleio.elem(negind) = -1 * Tpleio.elem(negind);
                    
            arma::uvec validind = arma::find(Tpleio >1.96);
            
            arma::vec gamma_exp2 = gamma_exp.elem(validind);
            arma::vec se_exp2 = se_exp.elem(validind);
            arma::vec gamma_out2 = gamma_out.elem(validind);
            arma::vec se_out2 = se_out.elem(validind);
            
            Rcout << "Finish test" << i << std::endl;

            int p = se_out2.size();
            arma::vec K_vec = arma::regspace(0,p-2);
            arma::vec Kall = arma::regspace(0,p-2);

            Rcpp::List randRes = mr_cMLC2(gamma_exp2, gamma_out2, se_exp2, se_out2, K_vec, random_start, maxit,n,Kall);
            
            FinalOut(i,0) = randRes[0];
            FinalOut(i,1) = randRes[1];
            //arma::vec tmpInvalid = randRes[4];
            //InvalidMat.col(i) = tmpInvalid;
        }
        
        
        thetaInit = theta2;
        seInit = se2;
        
        theta2 = arma::accu(FinalOut.col(0)) / nrep;
        se2 = arma::accu(arma::pow(FinalOut.col(1),2) ) / nrep;
        outcount = outcount + 1;
    }
    
    
    Rcpp::List out;
    out["theta"] = theta2;
    out["se2"] = se2;
    out["outcount"] = outcount;
    //out["InvalidIV"] = InvalidMat;

    return (out);
}




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List  cML_estimate_randomCtest(arma::vec &b_exp, arma::vec &b_out, arma::vec &se_exp,  arma::vec & se_out, const int K,  const int random_start, const int maxit) {
    
    const int p = b_exp.size();

    arma:: vec tmp = b_out / b_exp;
    double minTheta = tmp.min();
    double maxTheta = tmp.max();
    
    double rangeTheta = maxTheta - minTheta;
    double initial_theta = 0.0;
    
    arma::vec initial_mu;
    initial_mu.zeros(p);
    
    arma::vec invalidIV;
    invalidIV.zeros(p);
    
    double NegL;
    double theta;
    double sdTheta;
    
    int outj = 0;
    for (int j = 0; j <= random_start; j++) {
        
        if (j > 0) {
            arma::vec tmpInitial_theta = rangeTheta * arma::randu(1) + minTheta;
            initial_theta = tmpInitial_theta(0);
            initial_mu = arma::randn(p) % se_exp + b_exp;
        }
        
        Rcpp::List MLEresult = cML_estimateC(b_exp, b_out, se_exp, se_out, K, initial_theta, initial_mu, maxit);
        
        double thetatmp = MLEresult[0];
        arma::vec b_vec = MLEresult[1];
        arma::vec r_vec = MLEresult[2];
        
        double NegLtmp = arma::accu( arma::pow(b_exp - b_vec,2) / (2 * arma::pow(se_exp,2) ) ) + arma::accu(arma::pow(b_out - thetatmp * b_vec - r_vec,2) / (2 * arma::pow(se_out,2) ) );
        
        double sdThetaTmp = cML_SdThetaC(b_exp, b_out, se_exp,  se_out, thetatmp, b_vec, r_vec);;
        
        
        if(outj == 0) {
            NegL = NegLtmp;
            theta = thetatmp;
            invalidIV = r_vec;
            sdTheta = sdThetaTmp;
        }
        
        if(NegLtmp < NegL) {
            NegL = NegLtmp;
            theta = thetatmp;
            invalidIV = r_vec;
            sdTheta = sdThetaTmp;
        }
        
        outj = outj + 1;
        
    }
    
    Rcpp::List out;
    out["theta"] = theta;
    out["se"] = sdTheta;
    out["NegL"] = NegL;
    out["r_est"] = invalidIV;
    out["outj"] = outj;

    return (out);
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List  cML_SdThetaCtest(arma::vec &b_exp, arma::vec &b_out, arma::vec &se_exp,  arma::vec & se_out, double theta, arma::vec b_vec, arma::vec r_vec) {
    
    arma::uvec zeroind = arma::find(r_vec == 0.0);
    
    arma::vec b_vec2 = b_vec.elem(zeroind);
    arma::vec se_out2 = se_out.elem(zeroind);
    arma::vec b_out2 = b_out.elem(zeroind);
    arma::vec se_exp2 = se_exp.elem(zeroind);
    
    double varTheta;
    
    varTheta = 1 / ( arma::accu( arma::pow(b_vec2,2) / arma::pow(se_out2,2))  - arma::accu( arma::pow(2 * theta * b_vec2 - b_out2,2) / ( arma::pow(se_out2,4 ) % (1 / arma::pow(se_exp2,2) + theta * theta / arma::pow(se_out2,2)) ) ) )  ;
    
    
    double varTheta2 = arma::accu( arma::pow(b_vec2,2) / arma::pow(se_out2,2)) ;
    double varTheta3 = arma::accu( arma::pow(2 * theta * b_vec2 - b_out2,2) / ( arma::pow(se_out2,4 ) % (1 / arma::pow(se_exp2,2) + theta * theta / arma::pow(se_out2,2)) ) );
    
    Rcpp::List out;
    
    out["varTheta"] = varTheta;
    out["varTheta2"] = varTheta2;
    out["varTheta3"] = varTheta3;
    
    return (out);

}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double mr_CD_measureC(arma::vec b_exp, arma::vec b_out, arma::vec se2_exp,  arma::vec se2_out, arma::vec w) {
    
    // CHECK this
    double a = arma::accu(w % ((arma::pow(b_exp,2) - se2_exp) / se2_out) );
    double b = arma::accu(w % ((b_exp % b_out) / se2_out) ) ;
    
    double out = b / a;
    return(out);
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List  cML_CD_estimateC2(arma::vec &b_exp, arma::vec &b_out, arma::vec & se2_exp,  arma::vec & se2_out, arma::vec w, const int K,  double initial_theta, const int maxit) {
    
    const int p = b_exp.size();
    
    double theta = initial_theta;
    double thetaold = theta - 1;
    
    
    arma::vec vimportance;
    vimportance.zeros(p);
    
    arma::vec vbg;
    vbg.zeros(p);
    
    int iteindx = 0;
    
    double dlx = 1.0;
    if(K > 0) {
        while ( (dlx > 1e-5)  & (iteindx < maxit) ) {
            thetaold = theta;
            iteindx = iteindx + 1;
            
            vimportance = (arma::pow(b_out - theta * b_exp,2) ) / (se2_out) - theta * theta * se2_exp / se2_out;
            
            vimportance = w % vimportance;
            
            arma::uvec tmpindx = arma::sort_index(vimportance,"descend");
            arma::uvec tmpindx2 = tmpindx.tail(p-K);
            vbg = b_out - theta * b_exp;
            
            vbg(tmpindx2).zeros();
            
            
            arma::vec b_exp2 = b_exp.elem(tmpindx2);
            arma::vec se2_out2 = se2_out.elem(tmpindx2);
            arma::vec b_out2 = b_out.elem(tmpindx2);
            arma::vec se2_exp2 = se2_exp.elem(tmpindx2);
            arma::vec w2 = w.elem(tmpindx2);
            
            //update theta
            theta = mr_CD_measureC(b_exp2, b_out2, se2_exp2, se2_out2, w2);
            
            dlx = thetaold - theta;
            if(dlx <0) {
                dlx = -1 * dlx;
            }
        }
        
        // update vbg and muvec once we find optimal theta
        arma::uvec nonzeroind = arma::find(vbg != 0);
        //muvec.elem(nonzeroind) = b_exp.elem(nonzeroind);
        
        arma::vec tmpvbg = b_out - theta * b_exp;
        vbg.elem(nonzeroind) = tmpvbg.elem(nonzeroind);
        
    } else {
        theta = mr_CD_measureC(b_exp, b_out, se2_exp, se2_out, w);
    }
    
    Rcpp::List out;
    out["theta"] = theta;
    out["rvec"] = vbg;
    out["iteindx"] = iteindx;
    
    return (out);
    
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List  cML_CD_randomC(arma::vec &b_exp, arma::vec &b_out, arma::vec &se2_exp,  arma::vec & se2_out, arma::vec & w, const int K,  const int random_start, const int maxit,double initial_theta) {
    
    const int p = b_exp.size();
    
    arma:: vec tmp = b_out / b_exp;
    double minTheta = tmp.min();
    double maxTheta = tmp.max();
    
    double rangeTheta = maxTheta - minTheta;
    //double initial_theta = 0.0;
    
    arma::vec initial_mu;
    initial_mu.zeros(p);
    
    arma::vec invalidIV;
    invalidIV.zeros(p);
    
    double NegL = -999;
    double theta = -999;
    //double sdTheta;
    
    int outj = 0;
    for (int j = 0; j <= random_start; j++) {
        
        if (j > 0) {
            arma::vec tmpInitial_theta = rangeTheta * arma::randu(1) + minTheta;
            initial_theta = tmpInitial_theta(0);
        }
        
        Rcpp::List MLEresult = cML_CD_estimateC2(b_exp, b_out, se2_exp, se2_out, w, K, initial_theta, maxit);
        
        double thetatmp = MLEresult[0];
        //double sdThetaTmp = MLEresult[1];
        arma::vec r_vec = MLEresult[1];
        
        
        // this needs to be revised
        arma::uvec zeroind = arma::find(r_vec == 0.0);
        
        arma::vec b_exp2 = b_exp.elem(zeroind);
        arma::vec b_out2 = b_out.elem(zeroind);
        arma::vec se2_exp2 = se2_exp.elem(zeroind);
        arma::vec se2_out2 = se2_out.elem(zeroind);
        arma::vec w2 = w.elem(zeroind);
        
        
        // CHECK this
        //double NegLtmp = 0.5 * arma::accu( w2 % (arma::pow(b_out2 - thetatmp * b_exp2 ,2) / se2_out2)) - 0.5 * arma::accu( thetatmp * thetatmp * (w2 % se2_exp2 / se2_out2));
        
        double NegLtmp = 0.5 * arma::accu( w2 % arma::pow(b_out2 - thetatmp * b_exp2,2) / (thetatmp * thetatmp * se2_exp2 + se2_out2));

        
        if(outj == 0) {
            NegL = NegLtmp;
            theta = thetatmp;
            invalidIV = r_vec;
            //sdTheta = sdThetaTmp;
        }
        
        if(NegLtmp < NegL) {
            NegL = NegLtmp;
            theta = thetatmp;
            invalidIV = r_vec;
            //sdTheta = sdThetaTmp;
        }
        
        outj = outj + 1;
        
    }
    
    Rcpp::List out;
    out["theta"] = theta;
    //out["se"] = sdTheta;
    out["NegL"] = NegL;
    out["r_est"] = invalidIV;
    
    return (out);
}




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec  mr_cMLC_CD(arma::vec &b_exp, arma::vec &b_out, arma::vec &se2_exp,  arma::vec & se2_out, arma::vec K_vec, arma::vec w,  const int random_start, const int maxit, const int n) {
    
    const int M = K_vec.size();
    //const int p = b_exp.size();
    
    arma::vec theta_v;
    arma::vec sd_v;
    arma::vec l_v;
    // arma::mat invalid_mat(p,M);
    
    double initial_theta = 0.0;
    
    
    theta_v.zeros(M);
    sd_v.zeros(M);
    l_v.zeros(M);
    // invalid_mat.fill(0);
    
    for (int j = 0; j < M; j++) {
        int K = K_vec(j);
        Rcpp::List randRes = cML_CD_randomC(b_exp, b_out, se2_exp, se2_out, w, K, random_start, maxit,initial_theta);
        
        theta_v(j) = randRes[0];
        //initial_theta = theta_v(j);
        //sd_v(j) = randRes[1];
        l_v(j) = randRes[1];
        
        //arma::vec tmpInvalid = randRes[3];
        //invalid_mat.col(j) = tmpInvalid;
    }
    
    
    arma::vec BIC_v = log(n) * K_vec + 2 * l_v;
    BIC_v = BIC_v - BIC_v.min();
    arma::vec weight_v = arma::exp(-0.5 * BIC_v);
    weight_v = weight_v / arma::accu(weight_v);
    
    double MA_BIC_Theta = arma::accu(theta_v % weight_v);
    //double MA_BIC_se = arma::accu(weight_v % sqrt(arma::pow(sd_v2,2) + arma::pow(theta_v2 - MA_BIC_Theta,2) ) );
    
    arma::uword i = BIC_v.index_min();
    
    double BIC_Theta = theta_v(i);
    
    
    arma::vec AIC_v = 2 * K_vec + 2 * l_v;
    AIC_v = AIC_v - AIC_v.min();
    arma::vec weight_v2 = arma::exp(-0.5 * AIC_v);
    weight_v2 = weight_v2 / arma::accu(weight_v2);
    
    double MA_AIC_Theta = arma::accu(theta_v % weight_v2);
    //double MA_BIC_se = arma::accu(weight_v % sqrt(arma::pow(sd_v2,2) + arma::pow(theta_v2 - MA_BIC_Theta,2) ) );
    
    arma::uword i2 = AIC_v.index_min();
    
    double AIC_Theta = theta_v(i2);
    
    arma::vec out;
    out.zeros(4);
    out(0) = MA_BIC_Theta;
    out(1) = BIC_Theta;
    out(2) = MA_AIC_Theta;
    out(3) = AIC_Theta;
    
    return (out);
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat  mr_cMLC_CD_boot_efron(arma::vec &b_exp, arma::vec &b_out, arma::vec &se2_exp,  arma::vec & se2_out, Rcpp::List K_vec, arma::mat wAll,  const int random_start, const int maxit, const int n, const int nrep) {
    
    // Check input dimensions
    if (nrep > wAll.n_cols) {
        Rcpp::stop("nrep exceeds wAll dimensions");
    }
    
    if (K_vec.size() < nrep) {
        Rcpp::stop("K_vec list has fewer elements than nrep");
    }
    
    arma::mat out(nrep,4);
    
    for (int j = 0; j < nrep; j++) {
        arma::vec w = wAll.col(j);
        
        arma::uvec zeroind = arma::find(w != 0.0);
        
        arma::vec b_exp2 = b_exp.elem(zeroind);
        arma::vec b_out2 = b_out.elem(zeroind);
        arma::vec se2_exp2 = se2_exp.elem(zeroind);
        arma::vec se2_out2 = se2_out.elem(zeroind);
        arma::vec w2 = w.elem(zeroind);
        
        //int p = b_exp2.size();
        
        arma::vec K_vec2 = K_vec[j];
        
        arma::vec tmp = mr_cMLC_CD(b_exp2, b_out2, se2_exp2, se2_out2,  K_vec2, w2, random_start, maxit, n);
        //Rcout << "Finish test" << j << std::endl;
        
        out.row(j) = tmp.t();
    }
    
    return(out);
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat  mr_cMLC_CD_bootdir_efron(arma::mat &b_exp, arma::mat &b_out, arma::vec &se2_exp,  arma::vec & se2_out, Rcpp::List K_vec, arma::vec w,  const int random_start, const int maxit, const int n, const int nrep) {
    
    
    arma::mat out(nrep,4);
    
    for (int j = 0; j < nrep; j++) {
                
        arma::vec b_exp2 = b_exp.col(j);
        arma::vec b_out2 = b_out.col(j);
       
        //int p = b_exp2.size();
        
        arma::vec K_vec2 = K_vec[j];
        
        arma::vec tmp = mr_cMLC_CD(b_exp2, b_out2, se2_exp, se2_out,  K_vec2, w, random_start, maxit, n);
        //Rcout << "Finish test" << j << std::endl;
        
        out.row(j) = tmp.t();
    }
    
    return(out);
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List mr_cMLC_CD3(arma::vec &b_exp, arma::vec &b_out, arma::vec &se2_exp,  arma::vec & se2_out, arma::vec K_vec, arma::vec w,  const int random_start, const int maxit, const int n) {// output the invalid IVs
    
    const int M = K_vec.size();
    const int p = b_exp.size();
    
    arma::vec theta_v;
    arma::vec sd_v;
    arma::vec l_v;
    arma::mat invalid_mat(p,M);
    
    theta_v.zeros(M);
    sd_v.zeros(M);
    l_v.zeros(M);
    invalid_mat.fill(0);
    
    double initial_theta = 0.0;

    //Rcout << "Finish test" << std::endl;
    
    for (int j = 0; j < M; j++) {
        int K = K_vec(j);
        Rcpp::List randRes = cML_CD_randomC(b_exp, b_out, se2_exp, se2_out, w, K, random_start, maxit,initial_theta);
        
        theta_v(j) = randRes[0];
        //initial_theta = theta_v(j);
        //sd_v(j) = randRes[1];
        l_v(j) = randRes[1];
        
        arma::vec tmpInvalid = randRes[2];
        invalid_mat.col(j) = tmpInvalid;
    }
    
    //Rcout << "Finish test" << std::endl;
    
    arma::vec BIC_v = log(n) * K_vec + 2 * l_v;
    BIC_v = BIC_v - BIC_v.min();
    arma::vec weight_v = arma::exp(-0.5 * BIC_v);
    weight_v = weight_v / arma::accu(weight_v);
    
    double MA_BIC_Theta = arma::accu(theta_v % weight_v);
    //double MA_BIC_se = arma::accu(weight_v % sqrt(arma::pow(sd_v2,2) + arma::pow(theta_v2 - MA_BIC_Theta,2) ) );
    
    arma::uword i = BIC_v.index_min();
    
    double BIC_Theta = theta_v(i);
    
    
    arma::vec AIC_v = 2 * K_vec + 2 * l_v;
    AIC_v = AIC_v - AIC_v.min();
    arma::vec weight_v2 = arma::exp(-0.5 * AIC_v);
    weight_v2 = weight_v2 / arma::accu(weight_v2);
    
    double MA_AIC_Theta = arma::accu(theta_v % weight_v2);
    //double MA_BIC_se = arma::accu(weight_v % sqrt(arma::pow(sd_v2,2) + arma::pow(theta_v2 - MA_BIC_Theta,2) ) );
    
    arma::uword i2 = AIC_v.index_min();
    
    double AIC_Theta = theta_v(i2);
    
    arma::vec out;
    out.zeros(4);
    out(0) = MA_BIC_Theta;
    out(1) = BIC_Theta;
    out(2) = MA_AIC_Theta;
    out(3) = AIC_Theta;
    
    Rcpp::List out2;
    
    out2["theta"] = out;
    out2["Invalid_BIC"] = invalid_mat.col(i);
    out2["Invalid_AIC"] = invalid_mat.col(i2);
    
    return (out2);
}
