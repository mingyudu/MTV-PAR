





require(doSNOW)
require(doParallel)
#require(fields)
#require(Kendall)

cores <- parallel::detectCores()
cl <- makeSOCKcluster(cores)
registerDoSNOW(cl)
#mydata <- matrix(rnorm(8000*500), ncol=500)
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
    #lam_set = c(0.15,0.2,0.3,0.5,1,2)
    #lam_set = c(0.07,0.08,0.1,0.2,0.3,0.5,1)
    lam_set = c(0.04,0.05,0.07,0.08,0.1,0.2,0.3)
    #lam_set = c(0.04,0.05,0.07,0.08,0.1,0.2)
    
    ##if using constant penalty, set power=0
    power=1
    
    ##create arrays to store results
    dims = c(length(n_set), length(gam_set), length(sd_set), length(pois_set), length(lam_set))
    # victor-purpura distance values
    vp.rst = array(data = NA, dim = dims)
    #firing rate 
    fr.rst = array(data = NA, dim = dims)
    cp.sum = array(data = NA, dim = dims)
    
    rate_func <- function(x){
      draw_samples(x, N=1,  kernel_fn = se_kernel, length = max(x)/10)
    }
    
    
    ###run simulation
    for(q in 1:length(n_set)){
      for(r in 1:length(gam_set)){
        for(s in 1:length(sd_set)){
          for(t in 1:length(pois_set)){
            ##generate data using gaussian kernel rates
            set.seed(B)
            bump = c(300,700)
            sim = simulate.BP(n_set[q],p,gam_set[r],sd_set[s],pois_set[t], bump=bump)
            dat = sim$y
            true_cp = sim$true_cp
            true_ct = sim$c
            true_st = sim$true_st
            ##calculate true firing rate using gaussian smoothing estimate
            n = n_set[q]
            count0 = rep(0,n)
            true.sum = 0
            for(trial in 1:p){
              temp = true_cp[[trial]]
              true.sum = true.sum + length(temp)
              count0[temp] = count0[temp]+1
            }
            true_fr = smth.gaussian(count0, window = 51,tails = TRUE)/p*50
            
            mean_lam = NULL
            
            for(u in 1:length(lam_set)){
              #print(c(q,r,s,t,u))
              
              vp.set = NULL
              cp.set = NULL
              
              ##do gaussian smoothing using initial estimates
              st_gauss0 = rep(0,n)
              for(trial in 1:p){
                cp_init = estspike.gaussian(dat, gam = gam_set[r], lam = lam_set[u], trial = trial, power = 0,st_gauss=0)$cp
                st_gauss0[cp_init] = st_gauss0[cp_init]+1}
              st_gauss = smth.gaussian(st_gauss0, window = 51,tails = TRUE)
              
              ##scale values 
              st_gauss = st_gauss/max(st_gauss)
       
              
              for(trial in 1:p){
                #print(trial)
                gam = gam_set[r]
                lam = lam_set[u]
                
                #estimate trace and change points
                rst = estspike.gaussian(dat, gam = gam, lam = lam, trial = trial,  power = power, st_gauss=st_gauss)
                est.ct = rst$ct
                est.cp = rst$cp
                est.st = rst$st
                cp.set[[trial]] = est.cp
                
                
                vp.set = c(vp.set, vp.dis(true_cp[[trial]], est.cp, 0.01))
              }
              count = rep(0,n)
              sum1 = 0
              for(trial in 1:p){
                temp = cp.set[[trial]]
                sum1 = sum1 + length(temp)
                count[temp] = count[temp]+1
              }
              est_fr = smth.gaussian(count, window = 51,tails = TRUE)/p*50
              
              vp.rst[q,r,s,t,u] = mean(na.omit(vp.set))
              fr.rst[q,r,s,t,u] = mean(na.omit((est_fr - true_fr)^2))
              cp.sum[q,r,s,t,u] = sum1
            }
            
          }
        }
      }
    }
    ##output: vp.rst (vp distance); fr.rst(l2 norm), est_fr(estimated firing rate), true_fr(true firing rate)
    comb.result = list(vp.rst = vp.rst[1,1,1,1,], fr.rst = fr.rst[1,1,1,1,],lam_set = lam_set,est_fr=est_fr, true_fr = true_fr)}
#    , error=function(e){})}
close(pb)
stopCluster(cl)
save.image("static_sim.RData")








