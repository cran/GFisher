#' Generalized Fisher's p-value combination statistic.
#' @param p - vector of input p-values.
#' @param df - vector of degrees of freedom for inverse chi-square transformation for each p-value. If all df's are equal, it can be defined by the constant.
#' @param w - vector of weights.
#' @return GFisher statistic sum_i w_i*qchisq(1 - p_i, df_i).
#' @references Hong Zhang and Zheyang Wu. "Accurate p-Value Calculation for Generalized Fisher's Combination Tests Under Dependence", <arXiv:2003.01286>.
#' @examples
#' n = 10
#' pval = runif(n)
#' stat.GFisher(pval, df=2, w=1)
#' stat.GFisher(pval, df=rep(2,n), w=rep(1,n))
#' stat.GFisher(pval, df=1:n, w=1:n)
#' @export
#' @import stats
#' @importFrom methods is

stat.GFisher = function(p, df=2, w=1){
  if(length(w)>1){
    w = length(p)*w/sum(w)
  }else{
    w = 1
  }
  if(all(df==2)){
    pp_trans = -2*log(p)
  }else{
    pp_trans = qchisq(p, df=df, lower.tail=F)
    if(anyNA(pp_trans)){
      na_id = which(is.na(pp_trans))
      pp_trans[na_id] = qchisq(1-p[na_id], df=df[na_id])}
  }
  fisherstat = sum(w*pp_trans)
  return(fisherstat)
}

