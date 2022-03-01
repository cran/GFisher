#' Survival function of the generalized Fisher's p-value combination statistic.
#' @param q - observed GFisher statistic.
#' @param df - vector of degrees of freedom for inverse chi-square transformation for each p-value. If all df's are equal, it can be defined by the constant.
#' @param w - vector of weights.
#' @param M - correlation matrix of the input statistics.
#' @param p.type - "two" = two-sided p-values, "one" = one-sided p-values.
#' @param method - "MR" = simulation-assisted moment ratio matching, "HYB" = moment ratio matching by quadratic approximation, "GB" = Brown's method with calculated variance. See details in the reference.
#' @param nsim - number of simulation used in the "MR" method, default = 5e4.
#' @return p-value of the observed GFisher statistic.
#' @references Hong Zhang and Zheyang Wu. "Accurate p-Value Calculation for Generalized Fisher's Combination Tests Under Dependence", <arXiv:2003.01286>.
#' @examples
#' set.seed(123)
#' n = 10
#' M = matrix(0.3, n, n) + diag(0.7, n, n)
#' zscore = matrix(rnorm(n),nrow=1)%*%chol(M)
#' pval = 2*(1-pnorm(abs(zscore)))
#' gf1 = stat.GFisher(pval, df=2, w=1)
#' gf2 = stat.GFisher(pval, df=1:n, w=1:n)
#' p.GFisher(gf1, df=2, w=1, M=M, method="HYB")
#' \donttest{p.GFisher(gf1, df=2, w=1, M=M, method="MR", nsim=5e4)}
#' p.GFisher(gf2, df=1:n, w=1:n, M=M, method="HYB")
#' \donttest{p.GFisher(gf2, df=1:n, w=1:n, M=M, method="MR", nsim=5e4)}
#' @export
#' @import stats
#' @importFrom Matrix nearPD
#' @importFrom methods is
#'
p.GFisher = function(q, df, w, M, p.type="two", method="HYB", nsim=NULL){
  n = dim(M)[1]
  if(length(df)==1){df = rep(df, n)}
  if(length(w)==1){w = rep(w, n)}
  w = n*w/sum(w)
  if(method=="MR"){
    M_chol = try(chol(M))
    if(is(M_chol, "try-error")){
      M = as.matrix(nearPD(M, corr=TRUE)$mat)
      M_chol = chol(M)
    }
    if(!is(nsim, "numeric")){
      nsim = 5e4
    }
    znull = matrix(rnorm(n*nsim), ncol=n, nrow=nsim)%*%M_chol
    if(p.type=="two"){
      pnull = 2*pnorm(-abs(znull),0,1)
    }else{
      pnull = pnorm(znull,0,1, lower.tail = FALSE)
    }
    if(sd(df)==0){
      if(all(df==2)){
        pp_trans = -2*log(pnull)
      }else{
        pp_trans = qchisq(1-pnull, df=df[1])
      }
    }else{
      pp_trans = t(apply(pnull, 1, function(x)qchisq(1-x, df=df)))
    }
    fishernull = apply(pp_trans, 1, function(x)sum(x*w))
    MM = sapply(1:4, function(x)mean(fishernull^x))
    mu = MM[1]; sigma2 = (MM[2]-MM[1]^2)*nsim/(nsim-1)
    cmu3 = MM[3]-3*MM[2]*MM[1]+2*MM[1]^3
    cmu4 = MM[4]-4*MM[3]*MM[1]+6*MM[2]*MM[1]^2-3*MM[1]^4
    gm = cmu3/sigma2^(3/2)
    kp = cmu4/sigma2^2
    a = 9*gm^2/(kp-3)^2
    pval = pgamma((q-mu)/sqrt(sigma2)*sqrt(a)+a, shape=a,  scale = 1, lower.tail = FALSE)
  }else{
    GM = getGFisherGM(df, w, M, p.type)
    mu = sum(w*df); sigma2 = sum(GM)
    if(method=="HYB"){
      out_GM = getGFisherlam(df, w, M, GM, p.type)
      lam = out_GM$lam;
      c2 = sum(lam^2); c3 = sum(lam^3); c4 = sum(lam^4)
      gm = sqrt(8)*c3/(c2)^(3/2)
      kp = 12*c4/c2^2 + 3
      a = 9*gm^2/(kp-3)^2
    }else{ #Brown's method
      a = mu^2/sigma2
    }
    pval = pgamma((q-mu)/sqrt(sigma2)*sqrt(a)+a, shape=a,  scale = 1, lower.tail = FALSE)
  }
  return(pval)
}
