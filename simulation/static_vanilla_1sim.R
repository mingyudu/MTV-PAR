rm(list = ls())
i = 6 # lambda_set index
B = 17 # random seed index
vanilla = TRUE # vanilla or pruned DP
setwd('/home/exx/Desktop/MTV-PAR/simulation/')

library(smoother)
library(MASS)
library(zoo)
library(tictoc)
source("function_l0spike.R")
source('simulate_sptrain.R')
print('static simulation vanilla!')

##set simulation parameters
P=50 #the number of trials
T=1000 #the number of time points

#AR1 parameters
gam=0.96
sd=0.3
max_fr=0.2 #the maximum firing rate per time bin
bump = c(300,700) #the peaks of a firing rate function
lam_set = c(0.2, 0.3, 0.4, 0.5, 0.75, 1)

##create arrays to store results
# victor-purpura distance values for uniform-penalty
unif_vp = rep(0, length(lam_set))
tv_singletrial_vp = rep(0, length(lam_set))
tv_vp = rep(0, length(lam_set))

#firing rate for uniform-penalty
unif_fr2 = rep(NA, length(lam_set))
tv_fr2 = rep(NA, length(lam_set))
tv_singletrial_fr2 = rep(NA, length(lam_set))

unif_fr = array(NA, dim = c(length(lam_set),T))
tv_fr = array(NA, dim = c(length(lam_set),T))
tv_singletrial_fr = array(NA, dim = c(length(lam_set),T))

# bandwidth selected
win_len = rep(NA, length(lam_set))

# run simulation
set.seed(B)
cat('seed: ', B, '\n')

sim = simulate.BP(T, P, gam, sd, max_fr, bump=bump)
dat = sim$y #p-by-matrix of calcium traces
true_cp = sim$true_cp #a list of P vectors, with each vector records the spike times
true_ct = sim$c #P-by-T matrix, storing the true traces
true_st = sim$true_st #P-by-T matrix, storing the true spike status at each time point
obj=simulate.bump(T, lam_star = max_fr,  bump = bump)
true_fr=obj$lam_s*10
cat('seed ', B, ' simulation finished...\n')

#### perform method 1: use the uniform-penalty L0 algorithm to obtain an initial estimate
cat('seed ', B, ' lam_set ', i, '\n')
tmp_st=matrix(0, P, T)
tmp_st_tv_singletrial=matrix(0, P, T)
unif_est = NULL #initialize
tv_singletrial_est=NULL #initialize

for(trial in 1:P){
  tic(paste0('trial', trial))
  cat('seed ', B, 'lam_set ', i, 'trial ', trial, 'constant penalty\n')
  if(vanilla){
    unif_est[[trial]] = estspike.vanilla(dat, gam = gam, lam = lam_set[i], trial = trial, power = 0, st_gauss=0)
  }else{
    unif_est[[trial]] = estspike.gaussian(dat, gam = gam, lam = lam_set[i], trial = trial, power = 0, st_gauss=0)
  }
  toc()
  tmp_st[trial, unif_est[[trial]]$cp] =1 #matrix of indicators for change points
  #### evaluate the performance of the uniform-penalty L0 algorithm
  unif_vp[i] = unif_vp[i] + vp.dis(unif_est[[trial]]$cp, true_cp[[trial]], 0.01)
  #### conduct time-varying estimation based on current trial
  #first estimate firing rate based on the spike train of current trial, i.e., tmp_st[trial, ]
  win = bw.SJ(unif_est[[trial]]$cp) # bandwidth selection
  st_gauss=smth.gaussian(tmp_st[trial,], window = 2*win, alpha = 1, tails = TRUE)*50
  st_gauss = st_gauss/max(st_gauss)
  tic(paste0('trial', trial))
  cat('seed ', B, 'lam_set ', i, 'trial ', trial, 'TV-1\n')
  if(vanilla){
    tv_singletrial_est[[trial]] = estspike.vanilla(dat, gam = gam, lam = lam_set[i], trial = trial, power = 1, st_gauss=st_gauss)
  }else{
    tv_singletrial_est[[trial]] = estspike.gaussian(dat, gam = gam, lam = lam_set[i], trial = trial, power = 1, st_gauss=st_gauss)
  }
  toc()
  tmp_st_tv_singletrial[trial, tv_singletrial_est[[trial]]$cp]=1
  tv_singletrial_vp[i] = tv_singletrial_vp[i] + vp.dis(tv_singletrial_est[[trial]]$cp, true_cp[[trial]], 0.01)
}

#empirical firing rate based on uniform-penalty
unif_spk = unlist(lapply(1:P, function(i){unif_est[[i]]$cp}))
win_len[i] = bw.SJ(unif_spk)
unif_fr[i,] = smth.gaussian(colMeans(tmp_st), window = 2*win_len[i], alpha = 1, tails = TRUE)*50 #here means was used because tmp_st is binary

#empirical firing rate based on spikes estimates using individual trials
tv_singletrial_fr[i,] = smth.gaussian(colMeans(tmp_st_tv_singletrial), window = 2*win_len[i], alpha = 1, tails = TRUE)*50 #here means was used because tmp_st is binary
#compute L2 error for uniform-penalty
unif_fr2[i] = mean(na.omit((unif_fr[i,] - true_fr)^2))#L2 error
#compute L2 error for time-varying penalty based on individual trials
tv_singletrial_fr2[i] = mean(na.omit((tv_singletrial_fr[i,] - true_fr)^2))

#### perform the time-varying L0 estimation
#### compute the time-varying penalty
##scale values. Don't forget to rescale it
st_gauss=unif_fr[i,]/max(unif_fr[i,])

tmp_st_tv=matrix(0, P, T)
tv_est=NULL
for(trial in 1:P){
  tic(paste0('trial', trial))
  cat('seed ', B, 'lam_set ', i, 'trial ', trial, 'TV-all\n')
  if(vanilla){
    tv_est[[trial]] = estspike.vanilla(dat, gam = gam, lam = lam_set[i], trial = trial, power = 1, st_gauss=st_gauss)
  }else{
    tv_est[[trial]] = estspike.gaussian(dat, gam = gam, lam = lam_set[i], trial = trial, power = 1, st_gauss=st_gauss)
  }
  toc()
  tmp_st_tv[trial, tv_est[[trial]]$cp] =1 #matrix of indicators for change points
  #### evaluate the performance of the uniform varying L0 algorithm
  tv_vp[i] = tv_vp[i] + vp.dis(tv_est[[trial]]$cp, true_cp[[trial]], 0.01)
}

#empirical firing rate based on time-varying penalty from all trials
tv_fr[i,] = smth.gaussian(colMeans(tmp_st_tv), window = 2*win_len[i], alpha = 1, tails = TRUE)*50 #here means was used because tmp_st is binary
#compute L2 error for time-varying penalty based on individual trials
tv_fr2[i] = mean(na.omit((tv_fr[i,] - true_fr)^2))

if(vanilla){
  save.image(paste0('./result/static_vanilla_1sim_20240615_seed=', B, '.RData'))
}else{
  save.image(paste0('./result/static_pruned_1sim_20240615_seed=', B, '.RData'))
}

