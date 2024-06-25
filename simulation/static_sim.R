setwd('/home/exx/Desktop/MTV-PAR/simulation/')
require(doSNOW)
require(doParallel)

# cores <- parallel::detectCores()
cores <- 30
cl <- makeSOCKcluster(cores, outfile = '')
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
    library(tictoc)
    library(Rcpp)
    source("function_l0spike.R")
    source('simulate_sptrain.R')
    sourceCpp('RcppFunctions.cpp')
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
        tic(paste0('trial', trial))
        cat('seed ', B, 'lam_set ', i, 'trial ', trial, 'constant penalty\n')
        unif_est[[trial]] = estspike_vanilla(dat, gam = gam, lam = lam_set[i], trial = trial, power = 0, st_gauss=rep(0,T))
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
        tv_singletrial_est[[trial]] = estspike_vanilla(dat, gam = gam, lam = lam_set[i], trial = trial, power = 1, st_gauss=st_gauss)
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
      for(trial in 1:P)
      {
        tic(paste0('trial', trial))
        cat('seed ', B, 'lam_set ', i, 'trial ', trial, 'TV-all\n')
        tv_est[[trial]] = estspike_vanilla(dat, gam = gam, lam = lam_set[i], trial = trial, power = 1, st_gauss=st_gauss)
        toc()
        tmp_st_tv[trial, tv_est[[trial]]$cp] =1 #matrix of indicators for change points
        #### evaluate the performance of the uniform varying L0 algorithm
        tv_vp[i] = tv_vp[i] + vp.dis(tv_est[[trial]]$cp, true_cp[[trial]], 0.01)
      }
      
      #empirical firing rate based on time-varying penalty from all trials
      tv_fr[i,] = smth.gaussian(colMeans(tmp_st_tv), window = 2*win_len[i], alpha = 1, tails = TRUE)*50 #here means was used because tmp_st is binary
      #compute L2 error for time-varying penalty based on individual trials
      tv_fr2[i] = mean(na.omit((tv_fr[i,] - true_fr)^2))
    }
    
    ##output: vp:vp distance, fr2:L2 distance, fr: estimated firing rate, true_cp: true spikes, est: estimated spikes
    comb.result = list(unif_vp = unif_vp, tv_singletrial_vp = tv_singletrial_vp, tv_vp = tv_vp, 
                       unif_fr2 = unif_fr2, tv_fr2 = tv_fr2, tv_singletrial_fr2 = tv_singletrial_fr2, 
                       unif_fr = unif_fr, tv_singletrial_fr = tv_singletrial_fr, tv_fr = tv_fr,
                       lam_set = lam_set, true_fr = true_fr, dat = dat, true_cp = true_cp, win_len = win_len)
  }
close(pb)
stopCluster(cl)
save.image("./result/static_sim_vanilla.RData")



### Make plots
setwd('/home/exx/Desktop/MTV-PAR/simulation/')
load('result/static_sim_pruned.RData')
lam_set = c(0.2, 0.3, 0.4, 0.5, 0.75, 1.0)
bseed = 100
P=50
T=1000
# vp distance
vp1 = array(NA, dim = c(bseed, length(lam_set)))
vp2 = array(NA, dim = c(bseed, length(lam_set)))
vp3 = array(NA, dim = c(bseed, length(lam_set)))
# L2-norm distance
fr1 = array(NA, dim = c(bseed, length(lam_set)))
fr2 = array(NA, dim = c(bseed, length(lam_set)))
fr3 = array(NA, dim = c(bseed, length(lam_set)))
# bandwidth
win1 = array(NA, dim = c(bseed, length(lam_set)))
# estimated firing rate
frr1 = array(NA, dim = c(bseed, T))
frr2 = array(NA, dim = c(bseed, T))
frr3 = array(NA, dim = c(bseed, T))

for (i in 1:bseed) {
  vp1[i,] = result[[i]]$unif_vp
  vp2[i,] = result[[i]]$tv_singletrial_vp
  vp3[i,] = result[[i]]$tv_vp
  fr1[i,] = result[[i]]$unif_fr2
  fr2[i,] = result[[i]]$tv_singletrial_fr2
  fr3[i,] = result[[i]]$tv_fr2
  win1[i,] = result[[i]]$win_len
  frr1[i,] = result[[i]]$unif_fr[2,]
  frr2[i,] = result[[i]]$tv_singletrial_fr[2,]
  frr3[i,] = result[[i]]$tv_fr[2,]
}

vp1 = vp1/P
vp2 = vp2/P
vp3 = vp3/P
fr1 = sqrt(fr1)
fr2 = sqrt(fr2)
fr3 = sqrt(fr3)
win1 = 2*win1

colMeans(vp1)
colMeans(vp2)
colMeans(vp3)

colMeans(fr1)
colMeans(fr2)
colMeans(fr3)

colMeans(win1)

library(ggplot2)
library(latex2exp)
library(gridExtra)

vp.data=data.frame(
  lambda=as.factor(c(rep(lam_set, each=bseed), rep(lam_set, each=bseed), rep(lam_set, each=bseed))),
  VP=c(c(vp1), c(vp2), c(vp3)),
  methods=c(rep("constant", bseed*length(lam_set)), rep("TV-1", bseed*length(lam_set)), rep("TV-all", bseed*length(lam_set)))
)
plot1=ggplot(vp.data, aes(x=lambda, y=VP, fill=methods)) + geom_boxplot(width=0.3) + 
  ggtitle("VP distance")+xlab(TeX("$\\lambda$"))+ylab("VP")+
  stat_summary(
    fun = median,
    geom = 'line',
    aes(group = methods, colour = methods),
    position = position_dodge(width = 0.3) #this has to be added
  )
plot1

fr2.data=data.frame(
  lambda=as.factor(c(rep(lam_set, each=bseed), rep(lam_set, each=bseed), rep(lam_set, each=bseed))),
  MSE=c(c(fr1), c(fr2), c(fr3)),
  methods=c(rep("constant", bseed*length(lam_set)), rep("TV-1", bseed*length(lam_set)), rep("TV-all", bseed*length(lam_set)))
)

plot2=ggplot(fr2.data, aes(x=lambda, y=MSE, fill=methods, main="L2 Norm")) + 
  ggtitle("L2 Norm") + xlab(TeX("$\\lambda$"))+ylab("L2 Norm")+
  geom_boxplot(width=0.3) + 
  stat_summary(
    fun = median,
    geom = 'line',
    aes(group = methods, colour = methods),
    position = position_dodge(width = 0.3) #this has to be added
  )
plot2

win.data = data.frame(
  lambda=as.factor(c(rep(lam_set, each=bseed))),
  win=c(c(win1)),
  methods=c(rep("constant", bseed*length(lam_set)))
)

plot3=ggplot(win.data, aes(x=lambda, y=win, fill=methods, main="Bandwidth")) + 
  ggtitle("Bandwidth") + xlab(TeX("$\\lambda$"))+ylab("Bandwidth")+
  geom_boxplot(width=0.3) + 
  stat_summary(
    fun = median,
    geom = 'line',
    aes(group = methods, colour = methods),
    position = position_dodge(width = 0.3) #this has to be added
  ) + 
  theme(legend.position = "none")
plot3

png("static-sim-bw.png")
grid.arrange(plot1, plot2, plot3, ncol=1, nrow=3)
dev.off()

png("bw-selection-fr8.png")
plot(result[[1]]$true_fr, type="l", xlab="Time Point", ylab="Firing Rate")
lines(colMeans(frr1), col=2)
lines(colMeans(frr2), col=3)
lines(colMeans(frr3), col=4)
legend("topleft", col = c(1:4), legend = c("true","constant","tv-1","tv-all"), 
       lty = c(1,1,1,1), text.col=1:4, bty = 'n', cex = 0.6)
dev.off()

# check one simulation (spike train)
## observed fluorescence trace and the underlying spike trains
i = 3
marker_true_cp = sim$true_cp[[i]]/50
marker_tv_cp = tv_est[[i]]$cp/50
marker_unif_cp = unif_est[[i]]$cp/50
marker_tv_singletrial_cp = tv_singletrial_est[[i]]$cp/50

marker_height = -1
plot((1:1000)/50, sim$y[i,], 'l', ylim = c(-4, 7), 
     xlab = 'Time (second)', ylab = 'Calcium Trace')
# abline(v = c(6, 14))
segments(x0 = marker_tv_cp, y0 = rep(marker_height, length(marker_tv_cp)),
         x1 = marker_tv_cp, y1 = rep(marker_height + 0.5, length(marker_tv_cp)),
         col = "blue", lwd = 0.5)
segments(x0 = marker_tv_singletrial_cp, y0 = rep(marker_height-1, length(marker_tv_singletrial_cp)),
         x1 = marker_tv_singletrial_cp, y1 = rep(marker_height + 0.5 -1, length(marker_tv_singletrial_cp)),
         col = "green", lwd = 0.5)
segments(x0 = marker_unif_cp, y0 = rep(marker_height-2, length(marker_unif_cp)),
         x1 = marker_unif_cp, y1 = rep(marker_height + 0.5 -2, length(marker_unif_cp)),
         col = "red", lwd = 0.5)
segments(x0 = marker_true_cp, y0 = rep(marker_height-3, length(marker_true_cp)),
         x1 = marker_true_cp, y1 = rep(marker_height + 0.5 -3, length(marker_true_cp)),
         col = "black", lwd = 0.5)
legend("topleft", col = c(1:4), legend = c("True","JW","TV-1","TV-all"), 
       lty = c(1,1,1,1), text.col=1:4, bty = 'n', cex = 0.8)






