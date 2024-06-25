
##reference: https://projecteuclid.org/euclid.aoas/1542078052 Jewell and Witten (2018)
## The code implements algorithm 2 and 3 in the reference

library(zoo)

# ##calculate numerator of Cy (proposition 2 Jewell and Witten (2018))
# numer_C <- function(Y, gam){
#   rst = 0
#   for(i in 1:length(Y)){
#     rst = rst + Y[i]*(gam^(i-1))
#   }
#   return(rst)
# }
# 
# ##calculate Cy
# Cy<- function(Y, gam){
#   return(numer_C(Y, gam)/(1-gam^(2*length(Y)))*(1-gam^2))
# }
# 
# ##calculate Dy
# Dy<- function(Y,gam){
#   return(sum(Y^2)/2 - Cy(Y, gam)*numer_C(Y, gam) + (Cy(Y, gam))^2/2*(1-gam^(2*length(Y)))/(1-gam^2))
# }

#function of detecting spikes using constant penalty
#algorithm2 in Jewell and Witten (2018)
# estspike2 <- function(dat, gam, lam, trial){
#   ##input: decay parameter gamma(gam), penalty term lambda(lam), the trial of data (dat,trial)
#   ##initialize the change point sets
#   Y = dat[trial, ]
#   Fset = c((-1)*lam, rep(0,length(Y)))
#   cp = list()
#   cp[[1]] = 0
#   n = length(Y)
#   eps_s = 1
#   
#   pen = lam
#   
#   for (i in 2:(n+1)){
#     #print(i)
#     Fmin = Fset[1] + Dy(Y[1:(i-1)],gam) + pen
#     sprime = i-1
#     eps_s = c(eps_s,(i-1))
#     for (j in 1:(length(eps_s)-1)){
#       Fset.temp = Fset[eps_s[j]] + Dy(Y[eps_s[j]:(i-1)],gam) + pen
#       if(Fset.temp <= Fmin) {Fmin = Fset.temp; sprime = eps_s[j]}
#     }
#     
#     ex_idx = NULL
#     for (j in 1:(length(eps_s)-1)){
#       F_tau = Fset[eps_s[j]] + Dy(Y[eps_s[j]:(i-1)],gam)
#       if(F_tau >= Fmin) ex_idx = c(ex_idx, j)
#     }
#     if(length(ex_idx)>0) eps_s = eps_s[-ex_idx]
#     
#     Fset[i] = Fmin
#     cp[[i]] = unique(c(cp[[sprime]], sprime-1))
#   }
#   
#   cpset = cp[[n+1]]
#   
#   ##estimate calcium concentration ct
#   if(length(cpset) <= 1){
#     ct = Y
#     cpset = NULL
#   }
#   else{
#     cpset = cpset[2:length(cpset)]
#     ct = rep(0, length(Y))
#     for(i in 1:length(cpset)){
#       if(i == length(cpset)){
#         ct[cpset[i]+1] = Cy(Y[(cpset[i]+1):length(Y)], gam)
#         ct[(cpset[i]+2):length(Y)] = ct[cpset[i]+1]*(gam^seq(1,length(Y)-cpset[i]-1))
#       }
#       else{
#       ct[cpset[i]+1] = Cy(Y[(cpset[i]+1):cpset[i+1]], gam)
#       ct[(cpset[i]+2):cpset[i+1]] = ct[cpset[i]+1]*(gam^seq(1,cpset[i+1]-cpset[i]-1))
#       }
#     }
#     cpset = cpset+1
#   }
#   #use thresholded estimate of spike
#   st = ct[2:n] - gam*ct[1:(n-1)]
#   st = c(0,st)
#   #cpset = cpset[st[cpset]>L]
#   return(list(cp = cpset, ct = ct, st = st))
# }

##use auto-covariance function to estimate gamma
#reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4881387/pdf/nihms778164.pdf (equation 3)
# est.gam <- function(Y){
#   #Y = dat[trial,]
#   n = length(Y)
#   auto1 = sum((Y[1:(n-1)]-mean(Y))*(Y[2:n]-mean(Y)))/n
#   auto2 = sum((Y[1:(n-2)]-mean(Y))*(Y[3:n]-mean(Y)))/n
#   return(auto2/auto1)
# }

##function for cv(choosing lambda)
#algorithm 3 in Jewell and Witten (2018)
# choose.lam <- function(dat, gam, lamset, trial = trial){
#   Y = dat[trial,]
#   M = length(lamset)
#   n = length(Y)
#   cvMSE = matrix(0,M,2)
#   for(fold in 1:2){
#     if(fold == 1){
#       Y_train = Y[seq(1,length(Y),2)]
#       Y_test = Y[seq(2,length(Y),2)]
#     }
#     else{
#       Y_train = Y[seq(2,length(Y),2)]
#       Y_test = Y[seq(1,length(Y),2)]
#     }
#     for(m in 1:M){
#       #c_train = estspike2(t(Y_train), gam, lamset[m])$ct
#       c_train = estspike.gaussian(t(Y_train), gam, lamset[m], trial = 1, power = power, st_gauss = st_gauss)$ct
#       c_test = (c_train[1:(n/2-1)] + c_train[2:(n/2)])/2
#       #if(fold == 1)
#         {c_test = c(c_test, c_train[n/2])}
#       #else{c_test = c(c_train[1], c_test)}
#       cvMSE[m, fold] = 2/n*sum((Y_test-c_test)^2)
#     }
#   }
#   cv_bar = rowMeans(cvMSE)
#   mhat = which.min(cv_bar)
#   #se_cv = sqrt((cvMSE[,1] - cv_bar)^2/2 + (cvMSE[,2] - cv_bar)^2/2)
#   
#   return(list(lam = lamset[mhat], cv_MSE = cv_bar))
# }

##simulate homogeneous data with n time points and p trials
# simulate <- function(n,p, gam, poisMean, sd, seed){
#   set.seed(seed = seed)
#   c = matrix(0,p,n)
#   f = matrix(0,p,n)
#   s = rpois(n, poisMean)
#   for(i in 1:n){
#     if (i > 1) c[,i] = gam %*% c[,(i-1)] + matrix(s[i],p,1)
#     else c[,i] = c[,i] + s[i]
#     
#     f[,i] = c[,i] + matrix(rnorm(p,0,sd),p,1)
#     #f[,i] = c[,i]
#   }
#   return(list(f = f, c = c, true_cp = which(s!=0)))
# }

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

##estimate spikes using time-varying penalty
estspike.gaussian <- function(dat, gam, lam, trial = trial, power = power, st_gauss = st_gauss){
  ##power: the value "a" in the penalty function
  ##st_gauss: the input estimated firing rate
  
  ##initialize the change point sets
  Y = dat[trial, ]
  Fset = c(-lam, rep(0,length(Y))) # F: cost function
  cp = list()
  cp[[1]] = 0 # selected changepoints
  n = length(Y)
  eps_s = 1 # candidate changepoint set
  
  pen = lam
  pen1 = lam
  w_t = exp(-st_gauss^power) # calculate the weight function 
                             # based on the estimated firing rate
  lam_t = w_t/sum(w_t)*lam*n # update the penalty term
                             # based on the weight function
  ##use time-varying penalty term 
  if(power != 0){pen1 = lam_t[1]}
  
  for (i in 2:(n+1)){ # i==s+1: i = 2,3, ..., (n+1) => s = 1,2, ... , n
    Fmin = Fset[1] + Dy(Y[1:(i-1)],gam) + pen1 # compute F(0) + D(y(1:s)) + lambda(0)
    sprime = i-1 # set the initial s'
    eps_s = c(eps_s,(i-1)) # update epsilon_s set
    
    for (j in 1:(length(eps_s)-1)){ # go over only the candidate changepoints (pruned)
      if(power != 0){pen = lam_t[eps_s[j]]}
      Fset.temp = Fset[eps_s[j]] + Dy(Y[eps_s[j]:(i-1)],gam) + pen
      if(Fset.temp <= Fmin) {Fmin = Fset.temp; sprime = eps_s[j];} # obtain s'
    }
    
    ex_idx = NULL
    for (j in 1:(length(eps_s)-1)){
      F_tau = Fset[eps_s[j]] + Dy(Y[eps_s[j]:(i-1)],gam) # pruned
      if(F_tau >= Fmin) ex_idx = c(ex_idx, j)
    }
    if(length(ex_idx)>0) eps_s = eps_s[-ex_idx]
    
    Fset[i] = Fmin # calculate F(s); s = 1,2, ..., n
    cp[[i]] = unique(c(cp[[sprime]], sprime-1)) # update cp(s); s = 1,2, ..., n
  }
  
  cpset = cp[[n+1]] # cpset: changepoints (defined at the last time step)
  cpset = cpset+1 # spike times
  ##estimate calcium concentration ct
  if(length(cpset) <= 1){
    ct = rep(0, length(Y))
    ct[1] = Cy(Y, gam)
    ct[2:length(Y)] = ct[1]*(gam^seq(1,length(Y)-1))
    cpset = NULL
  }else{
    cpset = cpset[2:length(cpset)]
    ct = rep(0, length(Y))
    ct[1] = Cy(Y[1:cpset[1]], gam)
    ct[2:cpset[1]] = ct[1]*(gam^seq(1,cpset[1]-1))
    for(i in 1:length(cpset)){
      if(i == length(cpset)){
        ct[cpset[i]+1] = Cy(Y[(cpset[i]+1):length(Y)], gam)
        ct[(cpset[i]+1):length(Y)] = ct[cpset[i]+1]*(gam^seq(0,length(Y)-cpset[i]-1))
      }else{
        ct[cpset[i]+1] = Cy(Y[(cpset[i]+1):cpset[i+1]], gam)
        ct[(cpset[i]+1):cpset[i+1]] = ct[cpset[i]+1]*(gam^seq(0,cpset[i+1]-cpset[i]-1))
      }
    }
  }
  
  ## remove the spike location i where ct_{i+1} - ct_{i} < 0
  if(length(cpset) > 0){
    rm_idx = NULL
    for(i in 1:length(cpset)){
      if((ct[cpset[i]+1]-ct[cpset[i]])<0){rm_idx = c(rm_idx,i)}
    }
    if(length(rm_idx)>0)
    {cpset = cpset[-rm_idx]}
  }
  
  st = ct[2:n] - gam*ct[1:(n-1)]
  st = c(0,st)
  # cp: spike times
  return(list(cp = cpset, ct = ct, st = st,lam_t = lam_t,Fset=Fset,cpall=cp,eps_s=eps_s))
}


# ##estimate spikes using time-varying penalty (vanilla dynamic programming)
# estspike.vanilla <- function(dat, gam, lam, trial = trial, power = power, st_gauss = st_gauss){
#   ##power: the value "a" in the penalty function
#   ##st_gauss: the input estimated firing rate
#   
#   # Initialize the change point sets
#   Y = dat[trial,]
#   n = length(Y)
#   Fset = c(-lam, rep(0, n)) # length = n+1
#   cp = vector('list', n + 1) # preallocate the list
#   cp[[1]] = 0
#   
#   pen = lam
#   w_t = exp(-st_gauss^power)
#   lam_t = w_t/sum(w_t)*lam*n
#   # Use time-varying penalty term if power != 0
#   if(power != 0){
#     pen1 = lam_t[1]
#   }else{
#     pen1 = lam
#   }
#   
#   for (i in 2:(n+1)){ # i==s+1: i = 2,3, ..., (n+1) => s = 1,2, ... , n
#     Fmin = Fset[1] + Dy(Y[1:(i-1)], gam) + pen1 # D(y(1:s)) => s == i-1
#     sprime = 1  ##### fixed 06/18/2024: set the initial s' = 1 instead of i-1
#     
#     if (i > 2){ # i>2 => s>1 => s = 2,3, ...
#       for (j in 2:(i-1)){ # j == tau+1 
#         if(power != 0){
#           pen = lam_t[j]
#         }
#         Fset.temp = Fset[j] + Dy(Y[j:(i-1)],gam) + pen
#         if(Fset.temp <= Fmin){Fmin = Fset.temp; sprime = j;}
#       }
#     }
#     Fset[i] = Fmin
#     cp[[i]] = unique(c(cp[[sprime]], sprime-1))
#   }
#   
#   cpset = cp[[n+1]]
#   cpset = cpset+1
#   
#   # Estimate calcium concentration ct
#   if(length(cpset) <= 1){
#     ct = rep(0, length(Y))
#     ct[1] = Cy(Y, gam)
#     ct[2:length(Y)] = ct[1]*(gam^seq(1,length(Y)-1))
#     cpset = NULL
#   }else{
#     cpset = cpset[2:length(cpset)]
#     ct = rep(0, length(Y))
#     ct[1] = Cy(Y[1:cpset[1]], gam)
#     ct[2:cpset[1]] = ct[1]*(gam^seq(1,cpset[1]-1))
#     for(i in 1:length(cpset)){
#       if(i == length(cpset)){
#         ct[cpset[i]+1] = Cy(Y[(cpset[i]+1):length(Y)], gam)
#         ct[(cpset[i]+1):length(Y)] = ct[cpset[i]+1]*(gam^seq(0,length(Y)-cpset[i]-1))
#       }else{
#         ct[cpset[i]+1] = Cy(Y[(cpset[i]+1):cpset[i+1]], gam)
#         ct[(cpset[i]+1):cpset[i+1]] = ct[cpset[i]+1]*(gam^seq(0,cpset[i+1]-cpset[i]-1))
#       }
#     }
#   }
#   
#   ## Remove the spike location i where ct_{i+1} - ct_{i} < 0
#   if(length(cpset) > 0){
#     rm_idx <- which((ct[cpset + 1] - ct[cpset]) < 0)
#     if(length(rm_idx)>0){
#       cpset = cpset[-rm_idx]
#     }
#   }
#   
#   st = ct[2:n] - gam*ct[1:(n-1)]
#   st = c(0,st)
#   # cp: spike times
#   return(list(cp = cpset, ct = ct, st = st,lam_t = lam_t,Fset=Fset,cpall=cp))
# }
