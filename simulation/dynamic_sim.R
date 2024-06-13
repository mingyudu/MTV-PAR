require(doSNOW)
require(doParallel)

cores <- parallel::detectCores()
cl <- makeSOCKcluster(cores, outfile="")
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
    
    #simulation parameters
    P=50 #the number of trials
    T=1000 #the number of time points
    window_length = 51;
    
    #AR1 parameters
    gam=0.96
    sd=0.3
    max_fr=0.2 #the maximum firing rate per time bin
    bump = c(300,700) #the peaks of a firing rate function
    lam_set = c(0.2, 0.3, 0.4, 0.5, 0.75, 1.0) # lambda values
    
    #true firing rate functions
    true_fr = matrix(0, P, T)
    pow_vec = exp(-(1:P-(P/2))^2/T)
    for(trial in 1:P){
      pow = pow_vec[trial]
      # obj=simulate.bump2(T, lam_star = max_fr,  bump = bump, pow=pow_vec[trial])
      obj=simulate.bump2(T, bump = bump, pow=pow_vec[trial]) # modify the simulate.bump2 function, removing the argument: lam_star
      true_fr[trial,]=obj$lam_s*50  
    }#P-by-T matrix, storing the true firing rate at each time point
    
    ##create arrays to store results
    # victor-purpura distance values for uniform-penalty
    unif_vp = rep(0, length(lam_set))
    tv_singletrial_vp = rep(0, length(lam_set))
    tv_vp = rep(0, length(lam_set))
    
    #firing rate for uniform-penalty 
    unif_fr2 = rep(NA, length(lam_set))
    tv_fr2 = rep(NA, length(lam_set))
    tv_singletrial_fr2 = rep(NA, length(lam_set))
    
    # estimated firing rate: dynamic
    unif_fr = array(NA, dim = c(length(lam_set),P,T))
    tv_fr = array(NA, dim = c(length(lam_set),P,T))
    tv_singletrial_fr = array(NA, dim = c(length(lam_set),P,T))
    
    # bandwidth selected
    win_len = array(NA, dim = c(length(lam_set),P))
    
    # run simulation
    set.seed(B)
    cat('seed: ', B, '\n')
    #simulate.BP2 is for dynamic firing rate functions
    sim = simulate.BP2(T, P, gam, sd, max_fr, bump=bump)
    
    dat = sim$y #p-by-matrix of calcium traces
    true_cp = sim$true_cp #a list of P vectors, with each vector records the spike times
    true_ct = sim$c #P-by-T matrix, storing the true traces
    true_st = sim$true_st #P-by-T matrix, storing the true spike status at each time point
    cat('seed ', B, ' simulation finished...\n')
    
    #three methods
    #method 1: uniform-penalty L0 using individual trials
    #method 2: time-varying based on firing rate individual trials
    #method 3: time-varying based on firing rate from multiple trials
    
    
    for(i in 1:length(lam_set)){ #the ith lam value
      #### perform method 1: use the uniform-penalty L0 algorithm to obtain an initial estimate
      cat('seed ', B, ' lam_set ', i, '\n')
      tmp_st=matrix(0, P, T)
      tmp_st_tv_singletrial=matrix(0, P, T)
      unif_est = NULL #initialize 
      tv_singletrial_est=NULL #initialize
      for(trial in 1:P){
        cat('seed ', B, 'lam_set ', i, 'trial ', trial, 'uniform initialzed\n')
        unif_est[[trial]] = estspike.gaussian(dat, gam = gam, lam = lam_set[i], trial = trial, power = 0, st_gauss=0)
        tmp_st[trial, unif_est[[trial]]$cp] =1 #matrix of indicators for change points
        #### evaluate the performance of the uniform-penalty L0 algorithm
        unif_vp[i] = unif_vp[i] + vp.dis(unif_est[[trial]]$cp, true_cp[[trial]], 0.01)
        
        #### conduct time-varying estimation based on current trial
        #first estimate firing rate based on the spike train of current trial, i.e., tmp_st[trial, ]
        win = bw.SJ(unif_est[[trial]]$cp) # bandwidth selection
        st_gauss=smth.gaussian(tmp_st[trial,], window = 2*win, alpha = 1, tails = TRUE)*50
        # st_gauss=smth.gaussian(tmp_st[trial,], window = window_length, tails = TRUE)*50
        st_gauss = st_gauss/max(st_gauss)
        tv_singletrial_est[[trial]] = estspike.gaussian(dat, gam = gam, lam = lam_set[i], trial = trial, power = 1, st_gauss=st_gauss)
        tmp_st_tv_singletrial[trial, tv_singletrial_est[[trial]]$cp]=1
        tv_singletrial_vp[i] = tv_singletrial_vp[i] + vp.dis(tv_singletrial_est[[trial]]$cp, true_cp[[trial]], 0.01)
      }
      
      for(trial in 1:P)
      {
        cat('seed ', B, 'lam_set ', i, 'trial ', trial, 'uniform smoothing\n')
        trial_index =c(1:P)[abs(1:P-trial)<5]
        # bandwidth selection
        unif_spk = unlist(lapply(trial_index, function(i){unif_est[[i]]$cp}))
        win_len[i, trial] = bw.SJ(unif_spk)
        #empirical firing rate based on uniform-penalty
        # unif_fr[i, trial, ] = smth.gaussian(colMeans(tmp_st[trial_index, ]), window = window_length, tails = TRUE)*50 #here means was used because tmp_st is binary
        unif_fr[i, trial, ] = smth.gaussian(colMeans(tmp_st[trial_index, ]), window = 2*win_len[i, trial], alpha = 1, tails = TRUE)*50
        #empirical firing rate based on spikes estimates using individual trials
        # tv_singletrial_fr[i, trial, ] = smth.gaussian(colMeans(tmp_st_tv_singletrial[trial_index,]), window = window_length, tails = TRUE)*50 #here means was used because tmp_st is binary
        tv_singletrial_fr[i, trial, ] = smth.gaussian(colMeans(tmp_st_tv_singletrial[trial_index,]), window = 2*win_len[i, trial], alpha = 1, tails = TRUE)*50
      }
      
      #compute L2 error for uniform-penalty
      unif_fr2[i] = mean(na.omit((unif_fr[i,,] - true_fr)^2))#L2 error
      #compute L2 error for time-varying penalty based on individual trials
      tv_singletrial_fr2[i] = mean(na.omit((tv_singletrial_fr[i,,] - true_fr)^2))
      
      #### perform the time-varying L0 estimation 
      
      #### compute the time-varying penalty
      ##scale values. Don't forget to rescale it
      tmp_st_tv=matrix(0, P, T)
      tv_est=NULL
      
      for(trial in 1:P)
      {
        #### use dynamic estimate 
        cat('seed ', B, 'lam_set ', i, 'trial ', trial, 'tv initialzed\n')
        trial_index =c(1:P)[abs(1:P-trial)<5]
        # st_gauss = smth.gaussian(colMeans(tmp_st[trial_index, ]), window = window_length, tails = TRUE)*50 #here means was used because tmp_st is binary
        tv_spk = unlist(lapply(trial_index, function(i){unif_est[[i]]$cp}))
        win = bw.SJ(tv_spk) # actually the same as win_len[i,trial], just in order to make it clear
        st_gauss = smth.gaussian(colMeans(tmp_st[trial_index, ]), window = 2*win, alpha = 1, tails = TRUE)*50
        st_gauss= st_gauss/max(st_gauss)
        
        tv_est[[trial]] = estspike.gaussian(dat, gam = gam, lam = lam_set[i], trial = trial, power = 1, st_gauss=st_gauss)
        tmp_st_tv[trial, tv_est[[trial]]$cp] =1 #matrix of indicators for change points
        #### evaluate the performance of the uniform varying L0 algorithm
        tv_vp[i] = tv_vp[i] + vp.dis(tv_est[[trial]]$cp, true_cp[[trial]], 0.01)
      }
      
      #####$$$$$$$$$$$$$$$$$$$$$ need to use dynamic estimator
      #empirical firing rate based on time-varying penalty from all trials 
      for(trial in 1:P)
      {
        cat('seed ', B, 'lam_set ', i, 'trial ', trial, 'tv smoothing\n')
        trial_index =c(1:P)[abs(1:P-trial)<5]
        # tv_fr[i, trial,] = smth.gaussian(colMeans(tmp_st_tv[trial_index,]), window = window_length, tails = TRUE)*50 #here means was used because tmp_st is binary
        tv_fr[i, trial,] = smth.gaussian(colMeans(tmp_st_tv[trial_index,]), window = 2*win_len[i, trial], alpha = 1, tails = TRUE)*50
      }
      
      #compute L2 error for time-varying penalty based on individual trials
      tv_fr2[i] = mean(na.omit((tv_fr[i,,] - true_fr)^2))
    }
    
    ##output: vp:vp distance, fr2:L2 distance, fr: estimated firing rate, true_cp: true spikes, est: estimated spikes
    comb.result = list(unif_vp = unif_vp, tv_singletrial_vp = tv_singletrial_vp, tv_vp = tv_vp, 
                       unif_fr2 = unif_fr2, tv_fr2 = tv_fr2, tv_singletrial_fr2 = tv_singletrial_fr2, 
                       unif_fr = unif_fr, tv_singletrial_fr = tv_singletrial_fr, tv_fr = tv_fr,
                       win_len = win_len,
                       lam_set = lam_set, true_fr = true_fr, dat = dat, true_cp = true_cp)
  }
close(pb)
stopCluster(cl)
save.image("dynamic_sim_pruned.RData")





