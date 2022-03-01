#' P-value of the omnibus generalized Fisher's p-value combination test.
#' @param p - vector of input p-values.
#' @param DF - matrix of degrees of freedom for inverse chi-square transformation for each p-value. Each row represents a GFisher test.
#' @param W - matrix of weights. Each row represents a GFisher test.
#' @param M - correlation matrix of the input statistics.
#' @param p.type - "two" = two-sided p-values, "one" = one-sided p-values.
#' @param method - "MR" = simulation-assisted moment ratio matching, "HYB" = moment ratio matching by quadratic approximation, "GB" = Brown's method with calculated variance. See details in the reference.
#' @param combine - "cct" = oGFisher using the Cauchy combination method, "mvn" = oGFisher using multivariate normal distribution.
#' @param nsim - number of simulation used in the "MR" method, default = 5e4.
#' @return 1. p-value of the oGFisher test. 2. individual p-value of each GFisher test.
#' @references Hong Zhang and Zheyang Wu. "Accurate p-Value Calculation for Generalized Fisher's Combination Tests Under Dependence", <arXiv:2003.01286>.
#' @examples
#' set.seed(123)
#' n = 10
#' M = matrix(0.3, n, n) + diag(0.7, n, n)
#' zscore = matrix(rnorm(n),nrow=1)%*%chol(M)
#' pval = 2*(1-pnorm(abs(zscore)))
#' DF = rbind(rep(1,n),rep(2,n))
#' W = rbind(rep(1,n), 1:10)
#' p.oGFisher(pval, DF, W, M, p.type="two", method="HYB", combine="cct")
#' @export
#' @importFrom mvtnorm pmvnorm
#' @importFrom methods is

p.oGFisher = function(p, DF, W, M, p.type="two", method="HYB", combine="cct", nsim=NULL){
  out = stat.oGFisher(p, DF, W, M, p.type="two", method="HYB", nsim=NULL)
  if(combine=="cct"){
    thr = out$cct
    pval = pcauchy(thr,lower.tail = F)
  }else{
    thr = out$minp
    nd = dim(DF)[1]
    COR_GFisher = getGFisherCOR(DD=DF, W=W, M=M, p.type=p.type)
    pval = 1 - pmvnorm(lower=rep(-Inf,nd), upper=qnorm(1-thr), mean=rep(0,nd), corr=COR_GFisher)[1]
  }
  return(list(pval=pval, pval_indi=out$PVAL))
}
