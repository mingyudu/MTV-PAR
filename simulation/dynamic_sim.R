# module load R/4.2.1
# module load gcc/8.3.0
require(doSNOW)
require(doParallel)

# cores <- parallel::detectCores()
cores <- 24
cl <- makeSOCKcluster(cores, outfile='')
registerDoSNOW(cl)

pb <- txtProgressBar(min=1, max=100, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

t1 = Sys.time()
result <- 
  foreach(B=1:100,.options.snow=opts) %dopar% {
  tryCatch({
    library(smoother)
    library(MASS)
    library(zoo)
    library(tictoc)
    library(Rcpp)
    # source("function_l0spike.R")
    source("simulate_sptrain.R")
    sourceCpp("RcppFunctions.cpp")
    print("dynamic simulation vanilla!")
    
    #simulation parameters
    P=50 #the number of trials
    T=1000 #the number of time points
    ###window_length = 51; ### fixed window length. Not used!
    
    #AR1 parameter
    gam=0.96
    sd=0.3
    max_fr=0.2 #the maximum firing rate per time bin
    bump = c(300,700) #the peaks of a firing rate function
    lam_set = c(0.2, 0.3, 0.4, 0.5, 0.75, 1.0) # lambda values
    
    #true firing rate functions
    true_fr = matrix(0, P, T)
    pow_vec = exp(-(1:P-(P/2))^2/T) # varying peak height across trial
    for(trial in 1:P){
      pow = pow_vec[trial]
      obj=simulate.bump2(T, lam_star = max_fr,  bump = bump, pow=pow_vec[trial])
      true_fr[trial,]=obj$lam_s*10  
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
    sim = simulate.BP2(T, P, gam, sd, max_fr, bump=bump)
    
    dat = sim$y #p-by-matrix of calcium traces
    true_cp = sim$true_cp #a list of P vectors, with each vector records the spike times
    true_ct = sim$c #P-by-T matrix, storing the true traces
    true_st = sim$true_st #P-by-T matrix, storing the true spike status at each time point
    # cat('seed ', B, ' simulation finished...\n')
    
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
      
      for(trial in 1:P)
      {
        # incorporate between-trial information within neighboring B trials
        # B = 2*5 = 10
        trial_index =c(1:P)[abs(1:P-trial)<5] 
        # bandwidth selection
        unif_spk = unlist(lapply(trial_index, function(i){unif_est[[i]]$cp}))
        win_len[i, trial] = bw.SJ(unif_spk)
        #empirical firing rate based on uniform-penalty
        unif_fr[i, trial, ] = smth.gaussian(colMeans(tmp_st[trial_index, ]), window = 2*win_len[i, trial], alpha = 1, tails = TRUE)*50
        #empirical firing rate based on spikes estimates using individual trials
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
        tic(paste0('trial', trial))
        cat('seed ', B, 'lam_set ', i, 'trial ', trial, 'TV-10\n')
        trial_index =c(1:P)[abs(1:P-trial)<5]
        # tv_spk = unlist(lapply(trial_index, function(i){unif_est[[i]]$cp}))
        # win = bw.SJ(tv_spk) # actually the same as win_len[i,trial], just in order to make it clear
        
        ### estimate the firing rate using the neighbor trials only (TV-10)
        st_gauss = smth.gaussian(colMeans(tmp_st[trial_index, ]), window = 2*win_len[i,trial], alpha = 1, tails = TRUE)*50
        st_gauss= st_gauss/max(st_gauss)
        
        tv_est[[trial]] = estspike_vanilla(dat, gam = gam, lam = lam_set[i], trial = trial, power = 1, st_gauss=st_gauss)
        toc()
        tmp_st_tv[trial, tv_est[[trial]]$cp] =1 #matrix of indicators for change points
        #### evaluate the performance of the uniform varying L0 algorithm
        tv_vp[i] = tv_vp[i] + vp.dis(tv_est[[trial]]$cp, true_cp[[trial]], 0.01)
      }
      
      #####$$$$$$$$$$$$$$$$$$$$$ need to use dynamic estimator
      #empirical firing rate based on time-varying penalty from all trials 
      for(trial in 1:P)
      {
        trial_index =c(1:P)[abs(1:P-trial)<5]
        tv_fr[i, trial,] = smth.gaussian(colMeans(tmp_st_tv[trial_index,]), window = 2*win_len[i, trial], alpha = 1, tails = TRUE)*50
      }
      
      #compute L2 error for time-varying penalty based on individual trials
      tv_fr2[i] = mean(na.omit((tv_fr[i,,] - true_fr)^2))
    }
    
    ##output: vp:vp distance, fr2:L2 distance, fr: estimated firing rate, true_cp: true spikes, est: estimated spikes
    comb.result = list(unif_vp = unif_vp, tv_singletrial_vp = tv_singletrial_vp, tv_vp = tv_vp, 
                       unif_fr2 = unif_fr2, tv_fr2 = tv_fr2, tv_singletrial_fr2 = tv_singletrial_fr2, 
                       unif_fr = unif_fr, tv_singletrial_fr = tv_singletrial_fr, tv_fr = tv_fr,
                       lam_set = lam_set, true_fr = true_fr, dat = dat, true_cp = true_cp, win_len = win_len)
  }, error = function(e) return(paste0("The seed ", B, " caused the error: ", e)))
}
t2 = Sys.time()
print(t2-t1)

close(pb)
stopCluster(cl)
save.image("./result/20240709_dynamic_sim_vanilla_rcpp.RData")


load("./result/20240709_dynamic_sim_vanilla_rcpp.RData")
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
# estimated dynamic firing rate
frr1 = array(NA, dim = c(bseed, P, T))
frr2 = array(NA, dim = c(bseed, P, T))
frr3 = array(NA, dim = c(bseed, P, T))

for (i in 1:bseed) {
  vp1[i,] = result[[i]]$unif_vp
  vp2[i,] = result[[i]]$tv_singletrial_vp
  vp3[i,] = result[[i]]$tv_vp
  fr1[i,] = result[[i]]$unif_fr2
  fr2[i,] = result[[i]]$tv_singletrial_fr2
  fr3[i,] = result[[i]]$tv_fr2
  
  frr1[i,,] = result[[i]]$unif_fr[5,,]
  frr2[i,,] = result[[i]]$tv_singletrial_fr[5,,]
  frr3[i,,] = result[[i]]$tv_fr[5,,]
}

vp1 = vp1/P
vp2 = vp2/P
vp3 = vp3/P

colMeans(vp1)
colMeans(vp2)
colMeans(vp3)

colMeans(fr1)
colMeans(fr2)
colMeans(fr3)

library(ggplot2)
library(latex2exp)
library(patchwork)

vp.data=data.frame(
  lambda=as.factor(c(rep(lam_set, each=bseed), rep(lam_set, each=bseed), rep(lam_set, each=bseed))),
  VP=c(c(vp1), c(vp2), c(vp3)),
  methods=c(rep("JW", bseed*length(lam_set)), rep("TV-1", bseed*length(lam_set)), rep("TV-10", bseed*length(lam_set)))
)
plot1=
  ggplot(vp.data, aes(x=lambda, y=VP, col=methods)) + 
  geom_boxplot(width=0.4, position=position_dodge(width=0.6)) + 
  theme_bw() + 
  ggtitle("(b)")+xlab(TeX("$\\lambda$"))+ylab("VP")#+
# stat_summary(
#   fun = median,
#   geom = 'line',
#   aes(group = methods, colour = methods),
#   position = position_dodge(width = 0.3) #this has to be added
# )
plot1

fr2.data=data.frame(
  lambda=as.factor(c(rep(lam_set, each=bseed), rep(lam_set, each=bseed), rep(lam_set, each=bseed))),
  MSE=c(c(fr1), c(fr2), c(fr3)),
  methods=c(rep("JW", bseed*length(lam_set)), rep("TV-1", bseed*length(lam_set)), rep("TV-10", bseed*length(lam_set)))
)

plot2=
  ggplot(fr2.data, aes(x=lambda, y=MSE, col=methods)) + 
  geom_boxplot(width=0.4, position=position_dodge(width=0.6)) + 
  theme_bw() +
  ggtitle("(c)") + 
  xlab(TeX("$\\lambda$"))+ylab("L2 Norm")#+
# stat_summary(
#   fun = median,
#   geom = 'line',
#   aes(group = methods, colour = methods),
#   position = position_dodge(width = 0.3) #this has to be added
# )
plot2

library(dplyr)
library(reshape)

true_fr = result[[1]]$true_fr
unif_fr = colMeans(frr1) # average by seed number
tv_singletrial_fr = colMeans(frr2)
tv_fr = colMeans(frr3)
# calculate difference between true and estimated firing rates
diff_unif = (true_fr - unif_fr)
diff_tv_singletrial = (true_fr - tv_singletrial_fr)
diff_tv = (true_fr - tv_fr)
# create data frame for plot
df1 = melt(diff_unif)
df2 = melt(diff_tv_singletrial)
df3 = melt(diff_tv)
df1$model <- 'JW'
df2$model <- 'TV-1'
df3$model <- 'TV-10'
diff <- rbind(df1, df2, df3)
diff <- diff %>% mutate(model = factor(diff$model, 
                                       levels = c('JW', 'TV-1', 'TV-10')))
diff$X2 <- diff$X2/50 # divided by 50 Hz

# plot the heatmap for differences
plot3 = 
  ggplot(diff, aes(X2,X1,fill=value)) + 
  geom_tile() + 
  ggtitle('(d)') + 
  theme_minimal() + 
  facet_wrap(~model, nrow = 3) +
  scale_fill_viridis_b(limits=c(-0.31, 3.6), breaks=round(seq(-0.3,3.6,by=0.4),2)) +
  theme(legend.text = element_text(size = 7))+ # legend text font size
  xlab('Time (second)') +
  ylab("Trial Number") + 
  labs(fill = "Difference")
plot3

# heatmap: dynamic firing rate
library(reshape2)
df = reshape2::melt(true_fr)
colnames(df) = c('trial', 'time', 'fr')
df$time = df$time/50 # convert to seconds
plot4 = 
  df %>% 
  ggplot(aes(x = time, y = trial, fill = fr)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = '(a)', x = "Time (second)", y = "Trial Number", fill = 'Firing Rate') +
  theme_minimal()
plot4

fig = 
  plot4 / ((plot1 + plot2 + plot_layout(guides = 'collect', 
                                        nrow = 2) & 
              theme(legend.position = 'bottom')) | plot3 ) + 
  plot_layout(heights = c(1, 3))
fig

ggsave('20240719_dynamic_sim_accuracy_fig.png', plot = fig, 
       height = 7, width = 7, dpi = 1200)

# Ignore below
library(plotly)
# 3D plot for true firing rate
x = (1:1000)/50
y = 1:50
z = true_fr
fig1 <- plot_ly(type = 'surface', x = ~x, y = ~y, z = ~z, showscale=FALSE,
                colors = colorRamp(c("yellow","red", "black"))) %>% 
  plotly::layout(
    # title = 'True firing rate',
    scene = list(
      xaxis = list(title = 'Time (second)', titlefont = list(size=9), tickfont = list(size=10)),
      yaxis = list(title = 'Trial Number', titlefont = list(size=9), tickfont = list(size=10)),
      zaxis = list(title = 'Firing rate (Spikes/second)', titlefont = list(size=9), tickfont = list(size=10)),
      camera = list(eye = list(x = 1.25, y = 2, z = 1.25)),
      aspectratio = list(x=1.25,y=1,z=1)
    )
  )
fig1
orca(fig1,'true_fr.pdf')
# tried to stack but didn't work
fig2 <- ggplot(diff, aes(X2,X1,fill=value)) + 
  geom_tile() + 
  facet_wrap(~model, nrow = 3) +
  scale_fill_viridis_b(limits=c(0, 2.3), breaks=round(seq(0,2.3,by=0.2),2)) +
  theme(legend.text = element_text(size = 7))+ # legend text font size
  xlab('Time/sec') +
  ylab("Trial Number") + 
  ggtitle("Difference between the true and the estimated firing rate") + 
  labs(fill = "Difference")
fig2 <- ggplotly(fig2)
fig <- subplot(fig1, fig2, nrows=2)
# Warning messages: Can only have one: config 

