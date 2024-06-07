

##reference: https://projecteuclid.org/euclid.aoas/1542078052

library(zoo)

##calculate \sum y_t*gamma^(t-a)
Ygam_sum <- function(Y, gam){
  rst = 0
  for(i in 1:length(Y)){
    rst = rst + Y[i]*(gam^(i-1))
  }
  return(rst)
}

##Estimate \hat{Ca}
Ca<- function(Y, gam){
  n <- length(Y)
  if (n==1) {rst <- (Ygam_sum(Y, gam)-Y*gam)}
  else
  {rst <- (Ygam_sum(Y, gam)-mean(Y)*(1-gam^n)/(1-gam))/((1-gam^(2*n))/(1-gam^2)-(1-gam^n)^2/(1-gam)^2/n)}
  return(rst)
}

##Estimate beta_0
beta_0 <- function(Y,gam){
  n <- length(Y)
  if (n==1) {return(Y-Ca(Y,gam)*gam)}
  return(mean(Y)-Ca(Y,gam)*(1-gam^n)/(1-gam)/n)
}

##calculate Dy.beta
Dy.beta<- function(Y,gam){
  n <- length(Y)
  return(sum(Y^2)/2 - sum(Y)*beta_0(Y,gam)-Ygam_sum(Y,gam)*Ca(Y,gam)+Ca(Y,gam)^2*(1-gam^(2*n))/(1-gam^2)/2+
           n*beta_0(Y,gam)^2/2+Ca(Y,gam)*beta_0(Y,gam)*(1-gam^n)/(1-gam))
}


##use auto-covariance function to estimate gamma
#reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4881387/pdf/nihms778164.pdf (equation 3)
est.gam <- function(Y){
  Y = dat[trial,]
  n = length(Y)
  auto1 = sum((Y[1:(n-1)]-mean(Y))*(Y[2:n]-mean(Y)))/n
  auto2 = sum((Y[1:(n-2)]-mean(Y))*(Y[3:n]-mean(Y)))/n
  return(auto2/auto1)
}

##function for cv(choosing lambda)
#algorithm 3 in the paper
choose.lam <- function(dat, gam, lamset, trial = trial){
  Y = dat[trial,]
  M = length(lamset)
  n = length(Y)
  cvMSE = matrix(0,M,2)
  for(fold in 1:2){
    if(fold == 1){
      Y_train = Y[seq(1,length(Y),2)]
      Y_test = Y[seq(2,length(Y),2)]
    }
    else{
      Y_train = Y[seq(2,length(Y),2)]
      Y_test = Y[seq(1,length(Y),2)]
    }
    for(m in 1:M){
      #c_train = estspike2(t(Y_train), gam, lamset[m])$ct
      c_train = estspike.gaussian(t(Y_train), gam, lamset[m], trial = 1, power = power, st_gauss = st_gauss)$ct
      c_test = (c_train[1:(n/2-1)] + c_train[2:(n/2)])/2
      #if(fold == 1)
      {c_test = c(c_test, c_train[n/2])}
      #else{c_test = c(c_train[1], c_test)}
      cvMSE[m, fold] = 2/n*sum((Y_test-c_test)^2)
    }
  }
  cv_bar = rowMeans(cvMSE)
  mhat = which.min(cv_bar)
  #se_cv = sqrt((cvMSE[,1] - cv_bar)^2/2 + (cvMSE[,2] - cv_bar)^2/2)
  
  return(list(lam = lamset[mhat], cv_MSE = cv_bar))
}

##simulate homogeneous data with n time points and p trials
simulate <- function(n,p, gam, poisMean, sd, seed){
  set.seed(seed = seed)
  c = matrix(0,p,n)
  f = matrix(0,p,n)
  s = rpois(n, poisMean)
  for(i in 1:n){
    if (i > 1) c[,i] = gam %*% c[,(i-1)] + matrix(s[i],p,1)
    else c[,i] = c[,i] + s[i]
    
    f[,i] = c[,i] + matrix(rnorm(p,0,sd),p,1)
    #f[,i] = c[,i]
  }
  return(list(f = f, c = c, true_cp = which(s!=0)))
}

##calculate the victor-purpura metric
#reference: http://www-users.med.cornell.edu/~jdvicto/spkdm.html
vp.dis <- function(s1, s2, q){
  ##reference code in matlab: http://www-users.med.cornell.edu/~jdvicto/spkdm.html
  n1 = length(s1)
  n2 = length(s2)
  if(q == 0) return(abs(n1 - n2))
  else if(q == Inf) return(n1 + n2 - 2*length(intersect(s1, s2)))
  
  d_mat = matrix(0, n1+1, n2+1)
  d_mat[,1] = as.matrix(0:n1)
  d_mat[1,] = t(as.matrix(0:n2))
  if(sum(s1)&sum(s2)){
    for(i in 2:(n1+1))
      for(j in 2:(n2+1))
        d_mat[i,j]=min(d_mat[i-1,j]+1, d_mat[i,j-1]+1, d_mat[i-1,j-1]+q*abs(s1[i-1]-s2[j-1]));
  }
  return(d_mat[n1+1, n2+1])
}

##calculate the van-rossum distance
#reference: https://www.mitpressjournals.org/doi/pdf/10.1162/089976601300014321
#python reference: https://pythonhosted.org/fit_neuron/_modules/fit_neuron/evaluate/spkd_lib.html#van_rossum_dist
#source_python('van_rossum_dist.py')


estspike.inter <- function(dat, gam, lam, trial = trial, power = power, st_gauss = 0){
  ##initialize the change point sets
  Y = dat[trial, ]
  Fset = c(-lam, rep(0,length(Y)))
  cp = list()
  cp[[1]] = 0
  n = length(Y)
  eps_s = 1
  
  #if(power !=0){st_gauss = st_gauss}
  
  pen = lam
  pen1 = lam
  
  w_t = exp(-st_gauss^power)
  #w_t = 1/(st_gauss+1)
  lam_t = w_t/sum(w_t)*lam*n
  #lam_t = lam/(1+st_gauss^power)
  ##use time-varying penalty term 
  #if(power != 0){pen = lam*exp(-((st_gauss[1])^power))}
  if(power != 0){pen1 = lam_t[1]}
  #pen = lam_t[1]
  
  for (i in 2:(n+1)){
    Fmin = Fset[1] + Dy.beta(Y[1:(i-1)],gam) + pen1
    sprime = i-1
    eps_s = c(eps_s,(i-1))
    for (j in 1:(length(eps_s)-1)){
      #if(power != 0){pen = lam*exp(-((st_gauss[j])^power))}
      if(power != 0){pen = lam_t[eps_s[j]]}
      #pen = lam_t[eps_s[j]]
      Fset.temp = Fset[eps_s[j]] + Dy.beta(Y[eps_s[j]:(i-1)],gam) + pen
      if(Fset.temp <= Fmin) {Fmin = Fset.temp; sprime = eps_s[j];}
    }
    
    ex_idx = NULL
    for (j in 1:(length(eps_s)-1)){
      F_tau = Fset[eps_s[j]] + Dy.beta(Y[eps_s[j]:(i-1)],gam)
      if(F_tau >= Fmin) ex_idx = c(ex_idx, j)
    }
    if(length(ex_idx)>0) eps_s = eps_s[-ex_idx]
    
    Fset[i] = Fmin
    cp[[i]] = unique(c(cp[[sprime]], sprime-1))
    #print(Dy.beta(Y[eps_s[1]:(i-1)],gam))
  }
  
  cpset = cp[[n+1]]
  cpset = cpset+1
  #cpset = cpset[2:length(cpset)]
  
  if(length(cpset) <= 1){
    beta_vec = rep(beta_0(Y,gam),length(Y))
    ct = rep(0, length(Y))
    ct[1] = Ca(Y, gam)
    ct[2:length(Y)] = ct[1]*(gam^seq(1,length(Y)-1))
    cpset = NULL
  }else{
    cpset = cpset[2:length(cpset)]
    beta_vec = rep(0,length(Y))
    beta_vec[1:cpset[1]] = beta_0(Y[1:cpset[1]], gam)
    #beta_vec[1:cpset[1]] = 0
    ct = rep(0, length(Y))
    ct[1] = Ca(Y[1:cpset[1]], gam)
    ct[2:cpset[1]] = ct[1]*(gam^seq(1,cpset[1]-1))
    for(i in 1:length(cpset)){
      if(i == length(cpset)){
        beta_vec[(cpset[i]+1):length(Y)] = beta_0(Y[(cpset[i]+1):length(Y)], gam)
        ct[cpset[i]+1] = Ca(Y[(cpset[i]+1):length(Y)], gam)
        ct[(cpset[i]+1):length(Y)] = ct[cpset[i]+1]*(gam^seq(0,length(Y)-cpset[i]-1))
      }
      else{
        beta_vec[(cpset[i]+1):cpset[i+1]] = beta_0(Y[(cpset[i]+1):cpset[i+1]], gam)
        ct[cpset[i]+1] = Ca(Y[(cpset[i]+1):cpset[i+1]], gam)
        ct[(cpset[i]+1):cpset[i+1]] = ct[cpset[i]+1]*(gam^seq(0,cpset[i+1]-cpset[i]-1))
      }
    }
  }
  if(length(cpset) > 0){
    rm_idx = NULL
    for(i in 1:length(cpset)){
      if((ct[cpset[i]+1]+beta_vec[cpset[i]+1]-ct[cpset[i]]-beta_vec[cpset[i]])<0){rm_idx = c(rm_idx,i)}
    }
    if(length(rm_idx)>0)
    {cpset = cpset[-rm_idx]}
  }
  
  st = ct[2:n] - gam*ct[1:(n-1)]
  st = c(0,st)
  return(list(cp = cpset, ct = ct, beta_t = beta_vec))
}


