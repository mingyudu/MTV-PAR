
library(reticulate)
library(R.matlab)
library(smoother)
library(smoothie)
library(rlist)
library(plotly)
source('function_l0spike.R')
source("l0spike_extension.R")
source_python('est_gamma.py')

##estimate spikes without intercept using 1d smoothing
#dat: p(trial)-by-n matrix; lam: lambda;  iter:iteration of algorithm (iter=1 means constant L0)
#window: gaussian kernel 
###output: sp_train: estimated binary spike trains, fr_est: estimated firing rate, ct_est: estimated trace
##cp1: estimated spike locations; gam: estimated gamma
spike.ar <- function(dat,lam,iter,window){
  dat_avg <- colMeans(dat)
  #dat1 <- dat1 - dat_avg
  a = np_array(dat_avg, dtype = "float64")
  
  ##choose gamma value
  gam = c(estimate_parameters(a, p=1)[[1]])
  lam = lam
  p = dim(dat)[1]
  n = dim(dat)[2]
  ##do gaussian smoothing with iterations
  #lam = choose.lam(dat, gam, lamset, trial = trial)
  ct_est = matrix(0,p,n)
  for(it in 1:iter){
    sp0 = matrix(0,p,n)
    if(it==1){power=0;st_gauss=rep(0,n)}
    else power=1
    cp1 <- list()
    for(trial in 1:p){
      print(trial)
      rst = estspike.gaussian(dat, gam = gam, lam = lam, trial = trial, power = power,st_gauss = st_gauss)
      cp_init = rst$cp
      ct = rst$ct
      ct_est[trial,] = t(as.matrix(ct))
      cp1 <- list.append(cp1,cp_init)
      sp0[trial,cp_init] = sp0[trial,cp_init]+1}
    st_gauss = smth.gaussian(colSums(sp0), window = window,tails = TRUE)
    st_gauss = st_gauss/max(st_gauss)
  }
  fr_est = smth.gaussian(colSums(sp0), window = window,tails = TRUE)/p
  return(list(sp_train = sp0, fr_est=fr_est, cp1 = cp1,ct_est=ct_est,gam=gam))
}


##estimate spikes without intercept using 2d gaussian-boxcar smoothing
spike.2dar <- function(dat,lam,iter,window){
  dat_avg <- colMeans(dat)
  #dat1 <- dat1 - dat_avg
  a = np_array(dat_avg, dtype = "float64")
  
  ##choose gamma value
  gam = c(estimate_parameters(a, p=1)[[1]])
  lam = lam
  p = dim(dat)[1]
  n = dim(dat)[2]
  ##do gaussian smoothing with iterations
  ct_est = matrix(0,p,n)
  st_gauss = matrix(0,p,n)
  fr_est = matrix(0,p,n)
  for(it in 1:iter){
    sp0 = matrix(0,p,n)
    if(it==1){power=0;st_gauss=rep(0,n)}
    else power=1
    cp1 <- list()
    for(trial in 1:p){
      print(trial)
      rst = estspike.gaussian(dat, gam = gam, lam = lam, trial = trial, power = power,st_gauss = st_gauss[trial,])
      cp_init = rst$cp
      ct = rst$ct
      ct_est[trial,] = t(as.matrix(ct))
      cp1 <- list.append(cp1,cp_init)
      sp0[trial,cp_init] = sp0[trial,cp_init]+1}
    
    ##2d gaussian-boxcar
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
    fr_est = st_gauss
    st_gauss = st_gauss/max(st_gauss)
  }
  return(list(sp_train = sp0, fr_est=fr_est, cp1 = cp1,ct_est=ct_est,gam=gam))
}


##estimate spikes with intercept
#parameters similar to spike.ar
spike.inter <- function(dat,lam,iter,window){
  dat_avg <- colMeans(dat)
  #dat1 <- dat1 - dat_avg
  a = np_array(dat_avg, dtype = "float64")
  gam = c(estimate_parameters(a, p=1)[[1]])
  lam = lam
  p = dim(dat)[1]
  n = dim(dat)[2]
  ##do gaussian smoothing with iterations
  #lam = choose.lam(dat, gam, lamset, trial = trial)
  ct_est = matrix(0,p,n)
  for(it in 1:iter){
    sp0 = matrix(0,p,n)
    if(it==1){power=0;st_gauss=rep(0,n)}
    else power=1
    cp1 <- list()
    for(trial in 1:p){
      print(trial)
      rst = estspike.inter(dat, gam = gam, lam = lam, trial = trial, power = power,st_gauss = st_gauss)
      cp_init = rst$cp
      ct = rst$ct
      beta_t = rst$beta_t
      ct_est[trial,] = t(as.matrix(ct+beta_t))
      cp1 <- list.append(cp1,cp_init)
      sp0[trial,cp_init] = sp0[trial,cp_init]+1}
    st_gauss = smth.gaussian(colSums(sp0), window = window,tails = TRUE)
    st_gauss = st_gauss/max(st_gauss)
  }
  fr_est = smth.gaussian(colSums(sp0), window = window,tails = TRUE)/p
  return(list(sp_train = sp0, fr_est=fr_est, cp1 = cp1,ct_est=ct_est,gam=gam))
}


#### water lick data###############
cell_all <- readMat('allcell14.mat')
cell_all <- t(cell_all$dff.all)
cell_cr <- readMat('crcell14.mat')
cell_cr <- t(cell_cr$dff.cr)
cell_cl <- readMat('clcell14.mat')
cell_cl <- t(cell_cl$dff.cl)
cell_er <- readMat('ercell14.mat')
cell_er <- t(cell_er$dff.er)
cell_el <- readMat('elcell14.mat')
cell_el <- t(cell_el$dff.el)
timeonset <- readMat('watertime.mat')
timeonset <- timeonset$t

##plot firing rates
cr_fr = spike.ar(cell_cr, lam=1,iter=2,window = 10)$fr_est
cl_fr = spike.ar(cell_cl, lam=1,iter=2,window = 10)$fr_est
#er_fr = spike.ar(cell_er, lam=0.5,iter=2,window = 10)$fr_est
#el_fr = spike.ar(cell_el, lam=0.5,iter=2,window = 10)$fr_est
plot(timeonset,cr_fr*15,type = "l",col=4,xlab="Time from onset (s)",ylab = "Firing rate (Spikes/second)",main = "Estimated firing rates",
     ylim = c(0,4))
lines(timeonset,cl_fr*15,col=2,lty=1)
abline(h=0,lty=2)
abline(v=0,lty=2)
legend("topright", c("correct right", "correct left",'incorrect'), 
       col = c(4,2,1), lty = c(1,1,2),
       cex=0.8,  text.font = 2)



#####mouse C figures############
dat = readMat('neuron65_new.mat')
dat = (dat$dat.bytrial1)
p = dim(dat)[1]
n = dim(dat)[2]

load('mousec_ar099_200.RData')
sp0 = rst$sp_train
st_gauss0=sp0
###trace raster plot (figure 4) in Tong et al.(2021)
plot(1:n,dat[1,]/max(dat[1,])+1,type = "l",col=2,ylim = c(0,p+2),ylab = "Sessions",yaxt='n',cex.lab=1.5,
     xlab = "Time from shock start (s)",main = "MouseC neuron 65",xaxt='n')
axis(side=1,at=seq(1,1501,300),labels=seq(-40,60,20))
points(which(st_gauss0[1,]!=0),rep(0,length(which(st_gauss0[1,]!=0)))+1,pch=3,col=1)
for(i in 2:11){lines(1:n,dat[i,]/max(dat[i,])+i,col=2)
  points(which(st_gauss0[i,]!=0),rep(0,length(which(st_gauss0[i,]!=0)))+i,pch=3,col=1)
}
for(i in 12:p){lines(1:n,dat[i,]/max(dat[i,])+i,col=4)
  points(which(st_gauss0[i,]!=0),rep(0,length(which(st_gauss0[i,]!=0)))+i,pch=3,col=1)
}
abline(v=600,lty=2,col=2)


####2d firing rate plot 
n=dim(st_gauss0)[2];p=dim(st_gauss0)[1]

est_fr = matrix(0,p,n)
for(i in 1:p){
  est_fr[i,] = smth.gaussian(st_gauss0[i,], window = 100,tails = TRUE)*50
}

## between-trial smoothing using boxcar
for(i in 1:p){
  temp = est_fr[i,]
  for(j in 1:p){
    if(abs(i-j)<5&&i!=j){temp=rbind(temp,est_fr[j,])}
  }
  est_fr[i,] = colMeans(temp)
}


axx <- list(title = "Time from shockstart (s)")
axy <- list(title = "Trial number")
axz <- list(title = "Firing rate (Spikes/second)")
plot_ly(x = seq(-99/15,400/15,1/15), y = 1:p, z = est_fr[,501:1000]*15,colors = colorRamp(c("yellow","red", "black"))) %>%
  layout(title = 'Estimated firing rate (2D)', xaxis = list(
    ticktext = list("0", "10", "20"), 
    tickvals = seq(600,900,150),
    tickmode = "array"
  ), scene = list(xaxis=axx,yaxis=axy,zaxis=axz)) %>% add_surface()



