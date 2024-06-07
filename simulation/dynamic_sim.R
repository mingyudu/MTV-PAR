







require(doSNOW)
require(doParallel)
#require(fields)
#require(Kendall)

cores <- parallel::detectCores()
cl <- makeSOCKcluster(cores)
registerDoSNOW(cl)
pb <- txtProgressBar(min=1, max=100, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)
result <- 
  foreach(B=1:100,.options.snow=opts) %dopar% {
    #tryCatch({
    library(smoother)
    library(MASS)
    library(zoo)
    
    source("function_l0spike.R")
    source('simulate_sptrain.R')
    
    ##set simulation parameters
    #p: number of trials
    p = 50
    #n_set: length of series; gam_set: gamma value
    #sd_set: standard deviation; pois_set: poisson intensity
    n_set=c(1000);gam_set=0.96;sd_set=0.15;pois_set = c(0.2); 
    #lam_set" lambda values
    #lam_set = c(0.04,0.05,0.1,0.3,1,2,3)
    #lam_set = c(0.05,0.07,0.08,0.1,0.2,0.3,0.5)
    #lam_set = c(0.1,0.15,0.2,0.3,0.5,1)
    #lam_set = c(0.04,0.05,0.07,0.08,0.1,0.2,0.3)
    #lam_set = c(0.055,0.06,0.07,0.08,0.1,0.15,0.2)
    lam_set = c(0.04,0.05,0.07,0.1,0.2,0.3,1)
    
    ##if using constant penalty, set power=0
    power=1
    
    ##create arrays to store results
    dims = c(length(n_set), length(gam_set), length(sd_set), length(pois_set), length(lam_set))
    # victor-purpura distance values
    vp.rst = array(data = NA, dim = dims)
    #firing rate 
    fr.rst = array(data = NA, dim = dims)
    cp.sum = array(data = NA, dim = dims)
    
    
    ###run simulation
    for(q in 1:length(n_set)){
      for(r in 1:length(gam_set)){
        for(s in 1:length(sd_set)){
          for(t in 1:length(pois_set)){
            ##generate data using gaussian kernel rates
            set.seed(B)
            bump = c(300,700)
            sim = simulate.BP2(n_set[q],p,gam_set[r],sd_set[s],pois_set[t], bump=bump)
            dat = sim$y
            true_cp = sim$true_cp
            true_ct = sim$c
            true_st = sim$true_st
            ##calculate true firing rate using gaussian smoothing estimate
            n = n_set[q]
            count0 = matrix(0,p,n)
            true.sum = 0
            for(trial in 1:p){
              temp = true_cp[[trial]]
              true.sum = true.sum + length(temp)
              count0[trial,temp] = count0[trial,temp]+1
            }
            
            true_fr = matrix(0,p,n)
            for(i in 1:p){
              true_fr[i,] = smth.gaussian(count0[i,], window =100,tails = TRUE)*50
            }
            
            ## between-trial smoothing for each time point
            for(i in 1:p){
              temp = true_fr[i,]
              for(j in 1:p){
                if(abs(i-j)<5&&i!=j){temp=rbind(temp,true_fr[j,])}
              }
              true_fr[i,] = colMeans(temp)
            }
            
            #true_fr = lam_s*10
            
            for(u in 1:length(lam_set)){
              #print(c(q,r,s,t,u))
              
              #mse.set = NULL
              vp.set = NULL
              #vp.set1 = NULL
              cp.set = NULL
              
              ##do 2d smoothing using initial estimates
              sp0 = matrix(0,p,n)
              for(trial in 1:p){
                print(trial)
                cp_init = estspike2(dat, gam = gam_set[r], lam = lam_set[u], trial = trial, power = 0)$cp
                sp0[trial,cp_init] = sp0[trial,cp_init]+1}
              
              ##temporal smoothing using gaussian kernel
              st_gauss = matrix(0,p,n)
              for(i in 1:p){
                st_gauss[i,] = smth.gaussian(sp0[i,], window = 100,tails = TRUE)
              }
              
              ## between-trial smoothing using boxcar
              for(i in 1:p){
                temp = st_gauss[i,]
                for(j in 1:p){
                  if(abs(i-j)<5&&i!=j){temp=rbind(temp,st_gauss[j,])}
                }
                st_gauss[i,] = colMeans(temp)
              }
              
              ##scale values 
              st_gauss = st_gauss/max(st_gauss)
              
              
              
              for(trial in 1:p){
                #print(trial)
                gam = gam_set[r]
                lam = lam_set[u]
                
                #estimate trace and change points
                rst = estspike.gaussian(dat, gam = gam, lam = lam, trial = trial,  power = power, st_gauss=st_gauss[trial,])
                est.ct = rst$ct
                est.cp = rst$cp
                est.st = rst$st
                cp.set[[trial]] = est.cp
                
                
                vp.set = c(vp.set, vp.dis(true_cp[[trial]], est.cp, 0.01))
              }
              count = matrix(0,p,n)
              sum1 = 0
              for(trial in 1:p){
                temp = cp.set[[trial]]
                sum1 = sum1 + length(temp)
                count[trial,temp] = count[trial,temp]+1
              }
              
              est_fr = matrix(0,p,n)
              for(i in 1:p){
                est_fr[i,] = smth.gaussian(count[i,], window = 100,tails = TRUE)*50
              }
              
              ## between-trial smoothing using boxcar
              for(i in 1:p){
                temp = est_fr[i,]
                for(j in 1:p){
                  if(abs(i-j)<5&&i!=j){temp=rbind(temp,est_fr[j,])}
                }
                est_fr[i,] = colMeans(temp)
              }
              
              vp.rst[q,r,s,t,u] = mean(na.omit(vp.set))
              fr.rst[q,r,s,t,u] = mean(na.omit((est_fr - true_fr)^2))
              cp.sum[q,r,s,t,u] = sum1
            }
            
          }
        }
      }
    }
    ##output: vp.rst (vp distance); fr.rst(l2 norm), est_fr(estimated firing rate), true_fr(true firing rate)
    comb.result = list(vp.rst = vp.rst[1,1,1,1,],lam_set = lam_set, fr.rst = fr.rst[1,1,1,1,], est_fr=est_fr,true_fr = true_fr,count=count0)}
#    , error=function(e){})}
close(pb)
stopCluster(cl)
save.image("dynamic_sim.RData")





