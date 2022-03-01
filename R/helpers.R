getGFishercoef = function(D, M, p.type="two"){
  n = dim(M)[1]
  if(p.type=="two"){
    coeff2 = sapply(D, function(y)integrate(function(x)qchisq(pchisq(x^2,df=1),df=y)*dnorm(x)*(x^2-1), -8,8)$value)
    coeff4 = sapply(D, function(y)integrate(function(x)qchisq(pchisq(x^2,df=1),df=y)*dnorm(x)*(x^4-6*x^2+3), -8,8)$value)
    coeff6 = sapply(D, function(y)integrate(function(x)qchisq(pchisq(x^2,df=1),df=y)*dnorm(x)*(x^6-15*x^4+45*x^2-15), -8,8)$value)
    coeff8 = sapply(D, function(y)integrate(function(x)qchisq(pchisq(x^2,df=1),df=y)*dnorm(x)*(x^8-28*x^6+210*x^4-420*x^2+105), -8,8)$value)
    if(length(D)==1){
      coeff2 = rep(coeff2, n)
      coeff4 = rep(coeff4, n)
      coeff6 = rep(coeff6, n)
      coeff8 = rep(coeff8, n)
    }
    return(list(coeff2=coeff2,coeff4=coeff4,coeff6=coeff6,coeff8=coeff8))
  }else{
    coeff1 = sapply(D, function(y)integrate(function(x)qchisq(pnorm(x),df=y)*dnorm(x)*(x), -8,8)$value)
    coeff2 = sapply(D, function(y)integrate(function(x)qchisq(pnorm(x),df=y)*dnorm(x)*(x^2-1), -8,8)$value)
    coeff3 = sapply(D, function(y)integrate(function(x)qchisq(pnorm(x),df=y)*dnorm(x)*(x^3-3*x), -8,8)$value)
    coeff4 = sapply(D, function(y)integrate(function(x)qchisq(pnorm(x),df=y)*dnorm(x)*(x^4-6*x^2+3), -8,8)$value)
    if(length(D)==1){
      coeff1 = rep(coeff1, n)
      coeff2 = rep(coeff2, n)
      coeff3 = rep(coeff3, n)
      coeff4 = rep(coeff4, n)
    }
    return(list(coeff1=coeff1,coeff2=coeff2,coeff3=coeff3,coeff4=coeff4))
  }
}
getGFishercov = function(D1, D2,W1,W2, M, p.type="two", var.correct = TRUE){
  n = dim(M)[1]
  res1.coeff = getGFishercoef(D1, M, p.type=p.type)
  res2.coeff = getGFishercoef(D2, M, p.type=p.type)
  if(p.type=="two"){
    GM_cross = (M^2/2*res1.coeff$coeff2%*%t(res2.coeff$coeff2)+
                  M^4/24*res1.coeff$coeff4%*%t(res2.coeff$coeff4)+
                  M^6/720*res1.coeff$coeff6%*%t(res2.coeff$coeff6)+
                  M^8/40320*res1.coeff$coeff8%*%t(res2.coeff$coeff8))
    if(var.correct==TRUE){
      dd = diag(M)
      v1 = (res1.coeff$coeff2)^2/2+
        (res1.coeff$coeff4)^2/24+
        (res1.coeff$coeff6)^2/720+
        (res1.coeff$coeff8)^2/40320
      v2 = (res2.coeff$coeff2)^2/2+
        (res2.coeff$coeff4)^2/24+
        (res2.coeff$coeff6)^2/720+
        (res2.coeff$coeff8)^2/40320
      GM_cross = diag(sqrt(2*D1/v1),n,n)%*%GM_cross%*%diag(sqrt(2*D2/v2),n,n)
    }
  }else{
    GM_cross = (M*res1.coeff$coeff1%*%t(res2.coeff$coeff1)+
                  M^2/2*res1.coeff$coeff2%*%t(res2.coeff$coeff2)+
                  M^3/6*res1.coeff$coeff3%*%t(res2.coeff$coeff3)+
                  M^4/24*res1.coeff$coeff4%*%t(res2.coeff$coeff4))
    if(var.correct==TRUE){
      dd = diag(M)
      v1 = (res1.coeff$coeff1)^2+
        (res1.coeff$coeff2)^2/2+
        (res1.coeff$coeff3)^2/6+
        (res1.coeff$coeff4)^2/24
      v2 = (res2.coeff$coeff1)^2+
        (res2.coeff$coeff2)^2/2+
        (res2.coeff$coeff3)^2/6+
        (res2.coeff$coeff4)^2/24
      GM_cross = diag(sqrt(2*D1/v1),n,n)%*%GM_cross%*%diag(sqrt(2*D2/v2),n,n)
    }
  }
  return(sum(GM_cross*(W1%*%t(W2))))
}
getGFisherCOR = function(DD, W, M, var.correct=TRUE, p.type="two"){ # DD is a matrix, row represents a GFisher df
  m = dim(DD)[1] # m is the numer of GFisher
  COV = matrix(NA, ncol=m, nrow=m)
  for(i in 1:m){
    for(j in i:m){
      COV[i,j] = getGFishercov(DD[i,], DD[j,],W[i,],W[j,], M, var.correct=var.correct,
                               p.type=p.type)
    }
  }
  COV[lower.tri(COV)] = t(COV)[lower.tri(COV)]
  return(cov2cor(COV))
}
getGFisherGM = function(D, w, M, p.type="two"){
  n = dim(M)[1]
  if(p.type=="two"){
    coeff2 = sapply(D, function(y)integrate(function(x)qchisq(pchisq(x^2,df=1),df=y)*dnorm(x)*(x^2-1), -8,8)$value)
    coeff4 = sapply(D, function(y)integrate(function(x)qchisq(pchisq(x^2,df=1),df=y)*dnorm(x)*(x^4-6*x^2+3), -8,8)$value)
    coeff6 = sapply(D, function(y)integrate(function(x)qchisq(pchisq(x^2,df=1),df=y)*dnorm(x)*(x^6-15*x^4+45*x^2-15), -8,8)$value)
    coeff8 = sapply(D, function(y)integrate(function(x)qchisq(pchisq(x^2,df=1),df=y)*dnorm(x)*(x^8-28*x^6+210*x^4-420*x^2+105), -8,8)$value)
    GM = M^2/2*coeff2%*%t(coeff2)+M^4/24*coeff4%*%t(coeff4) +
      M^6/720*coeff6%*%t(coeff6)+M^8/40320*coeff8%*%t(coeff8)
    GM = diag(sqrt(2*D)*w,n,n)%*%cov2cor(GM)%*%diag(sqrt(2*D)*w,n,n)
  }else{
    coeff3 = sapply(D, function(y)integrate(function(x)qchisq(pnorm(x),df=y)*dnorm(x)*(x^3-3*x), -8,8)$value)
    coeff4 = sapply(D, function(y)integrate(function(x)qchisq(pnorm(x),df=y)*dnorm(x)*(x^4-6*x^2+3), -8,8)$value)
    coeff2 = sapply(D, function(y)integrate(function(x)qchisq(pnorm(x),df=y)*dnorm(x)*(x^2-1), -8,8)$value)
    coeff1 = sapply(D, function(y)integrate(function(x)qchisq(pnorm(x),df=y)*dnorm(x)*(x), -8,8)$value)
    GM = M*coeff1%*%t(coeff1)+M^2/2*coeff2%*%t(coeff2)+
      M^3/6*coeff3%*%t(coeff3)+M^4/24*coeff4%*%t(coeff4)
    GM = diag(sqrt(2*D)*w,n,n)%*%cov2cor(GM)%*%diag(sqrt(2*D)*w,n,n)
    #GM[GM<0]=0
  }
  return(GM)
}
getGFisherlam = function(D, w, M, GM, p.type="two"){
  n = dim(M)[1]
  DM = pmin(matrix(rep(D,n),nrow=n,byrow = FALSE), matrix(rep(D,n),nrow=n,byrow = TRUE))
  M_tilde = sqrt(abs(GM)/DM/2)*sign(M);
  if(any(M_tilde>1)){
    M_tilde[M_tilde>1] = 0.999
  }
  if(any(eigen(M_tilde)$values<1e-10)){
    M_tilde = as.matrix(nearPD(M_tilde, corr = TRUE)$mat)
  }
  M_tilde_chol = chol(M_tilde)

  WM_chol = M_tilde_chol%*%diag(sqrt(w))#M_chol#
  tb = table(D)
  tbval = as.numeric(names(tb))
  lam =eigen(t(WM_chol)%*%WM_chol, symmetric = TRUE, only.values = TRUE)$values
  if(max(D)>1){
    for(i in 2:max(D)){
      Ai = diag(n)
      Did = which(D<i)
      if(length(Did)>0){
        diag(Ai)[Did] = 0
      }
      lam = c(lam,eigen((WM_chol)%*%Ai%*%t(WM_chol), symmetric = TRUE, only.values = TRUE)$values)
      #print(diag(diag(sqrt(w))%*%(Ai)%*%diag(sqrt(w))))
    }
  }
  lam = lam[lam>1e-10]
  return(list(lam=lam))
}

stat.oGFisher = function(p, DF, W, M, p.type="two", method="HYB", nsim=NULL){
  STAT = sapply(1:dim(DF)[1],function(x) stat.GFisher(p = p, df=DF[x,], w=W[x,]))
  PVAL = sapply(1:length(STAT),function(x)p.GFisher(STAT[x], df=DF[x,], w=W[x,], M, p.type, method, nsim))
  minp = min(PVAL)
  # check if there are very small non-zero p values and calcuate the cauchy statistics
  is.small<-(PVAL<1e-15)
  CCTSTAT = PVAL
  CCTSTAT[!is.small]<-tan((0.5-CCTSTAT[!is.small])*pi)
  CCTSTAT[is.small]<-1/CCTSTAT[is.small]/pi
  cct = mean(CCTSTAT)
  #cct = mean(tan(pi*(0.5-PVAL)))
  return(list(STAT=STAT, PVAL=PVAL, minp=minp, cct=cct))
}
