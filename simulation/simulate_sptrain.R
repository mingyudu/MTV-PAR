
##reference: https://www.r-bloggers.com/sampling-paths-from-a-gaussian-process/
# library(MASS)
#source('sample_fromGP.R')
#sigmoid function
# sigmoid <- function(z){
#   (1 + exp(-z))^(-1)
# }

# generate covariance matrix for points in `x` using given kernel function
# cov_matrix <- function(x, kernel_fn, ...) {
#   outer(x, x, function(a, b) kernel_fn(a, b, ...))
# }


# se_kernel <- function(x, y, sigma = 1, length = 1) {
#   sigma^2 * exp(- (x - y)^2 / (2 * length^2))
# }

# period_kernel <- function(x, y, p = 100, sigma = 1, length = 1) {
#   sigma^2 * exp(-2 * sin(pi * abs(x - y) / p)^2 / length^2)
# }

# rq_kernel <- function(x, y, alpha = 1, sigma = 1, length = 1) {
#   sigma^2 * (1 + (x - y)^2 / (2 * alpha * length^2))^(-alpha)
# }


# matern_kernel <- function(x, y, nu = 1.5, sigma = 1, l = 1) {
#   if (!(nu %in% c(0.5, 1.5, 2.5))) {
#     stop("p must be equal to 0.5, 1.5 or 2.5")
#   }
#   p <- nu - 0.5
#   d <- abs(x - y)
#   if (p == 0) {
#     sigma^2 * exp(- d / l)
#   } else if (p == 1) {
#     sigma^2 * (1 + sqrt(3)*d/l) * exp(- sqrt(3)*d/l)
#   } else {
#     sigma^2 * (1 + sqrt(5)*d/l + 5*d^2 / (3*l^2)) * exp(-sqrt(5)*d/l)
#   }
# }

# draw_samples <- function(x, N,  kernel_fn, ...) {
#   Y <- matrix(NA, nrow = length(x), ncol = N)
#   #set.seed(seed)
#   for (n in 1:N) {
#     K <- cov_matrix(x, kernel_fn, ...)
#     Y[, n] <- mvrnorm(1, mu = rep(0, times = length(x)), Sigma = K)
#   }
#   Y
# }


##sample from SCGP
# simulate.SGCP <- function(n, lam_star, rate_func = rate_fnunc){
#   #set.seed(1)
#   gs <- rate_func(1:n)
#   #gs <- c(rep(0.1,n/2),rep(10,n/2))
#   #gs <- draw_samples(1, 1, kernel_fn = kernel, length = n/10)
#   lam_s = sigmoid(gs)
#   #lam_s = c(rep(0.1,n/2), rep(1,n/2))
#   J = rpois(1, n*lam_star)
#   sj_set = sort(c(sample.int(n, J, replace = FALSE)))
#   eps_set = NULL
#   for(j in 1:J){
#     #set.seed(j)
#     rj = runif(1)
#     if(rj < lam_s[sj_set[j]]) eps_set = c(eps_set, sj_set[j])
#   }
#   return(list(eps_set = eps_set, sj_set = sj_set))
# }

##generate simulation data using gaussian kernel function
# simulate.GP <- function(n,p, gam, sd, lam_star,rate_func){
#   c = matrix(0,p,n)
#   y = matrix(0,p,n)
#   s = matrix(0,p,n)
#   eps = matrix(rnorm(p*n,0,sd),p,n)
#   cp = NULL
#   for(p0 in 1:p){
#     #set.seed(p0)
#     #eps = rnorm(n,0,sd)
#     true_cp = simulate.SGCP(n = n, lam_star = lam_star, rate_func = rate_func)$eps_set
#     s[p0, true_cp] = 1
#     s[p0, 1] = 0
#     cp[[p0]] = true_cp
#     for(i in 1:n){
#       if (i > 1) c[p0,i] = gam * c[p0,(i-1)] + s[p0,i]
#       else c[p0,i] = c[p0,i] + s[p0,i]
#       #y[p0,i] = c[p0,i]
#       #f[,i] = c[,i]
#     }
#   }
#   y = c + eps
#   return(list(y = y, c = c, true_cp = cp, true_st = s))
# }


###use bump functions as firing rate--------------------
# get.bump <- function(width){
#   y = exp(-width^2/(width^2-(x-width)^2))
#   return(y)
# }

## create bump in firing rate function
# 50Hz
# x = 1-1000 => t = 1/50, 2/50, ..., 1000/50
# d = 1
get.bump1 <- function(x, h0){
  exp(-abs(x-h0)^2/50^2)
}

# inv.sigmoid <- function(x){
#   log(x/(1-x))
# }

### Appendix B Algorithm:
### Simulate spike trains from an Inhomogeneous Poisson Process
# V: number of timesteps
# lam_star: upper bound of rate intensity
# bump: the locations where the peaks of firing rate occur
simulate.bump <- function(V, lam_star, bump = bump){
  lam_s = rep(0.05, V)
  for(i in 1:length(bump)){
    lam_s =  lam_s + 0.95*get.bump1(1:V, bump[i])
  }
  J = rpois(1, V*lam_star)
  sj_set = sort(c(sample.int(V, J, replace = FALSE)))
  eps_set = NULL
  for(j in 1:J){
    rj = runif(1)
    if(rj < lam_s[sj_set[j]]) eps_set = c(eps_set, sj_set[j])
  }
  return(list(eps_set = eps_set, sj_set = sj_set, lam_s = lam_s))
}


###generate simulation data under static firing rate
# n: number of timesteps
# p: number of trials
# gam: decay rate
# sd: standard deviation of the random noise
# lam_star: upper bound of rate intensity
# bump: the locations where the peaks of firing rate occur
simulate.BP <- function(n,p, gam, sd, lam_star, bump){
  c = matrix(0,p,n)
  y = matrix(0,p,n)
  s = matrix(0,p,n)
  eps = matrix(rnorm(p*n,0,sd),p,n)
  cp = NULL
  for(p0 in 1:p){
    true_cp = simulate.bump(V = n, lam_star = lam_star,  bump = bump)$eps_set
    s[p0, true_cp] = 1
    s[p0, 1] = 0
    cp[[p0]] = true_cp
    for(i in 1:n){
      if (i > 1) c[p0,i] = gam * c[p0,(i-1)] + s[p0,i]
      else c[p0,i] = c[p0,i] + s[p0,i]
    }
  }
  y = c + eps
  true_fr = simulate.bump(V = n, lam_star = lam_star,  bump = bump)$lam_s
  # y: observed fluorescence trace
  # c: underlying calcium concentration
  # true_cp: true spike locations (a vector of index)
  # true_st: true spike locations (a binary matrix)
  # true_fr: true firing rate
  return(list(y = y, c = c, true_cp = cp, true_st = s,true_fr=true_fr))
}

#--------------------------
# Simulation example
#--------------------------
# set.seed(1234)
# sim = simulate.BP(n = 1000, p = 50, gam = 0.96, sd = 0.3,
#                   lam_star = 0.2, bump = c(300, 700))
## observed fluorescence trace and the underlying spike occurrence
# i = 3
# marker_times = sim$true_cp[[i]]/50
# marker_height = -1
# plot((1:1000)/50, sim$y[i,], 'l', ylim = c(-1,8))
# abline(v = c(6, 14))
# segments(x0 = marker_times, y0 = rep(marker_height, length(marker_times)),
#          x1 = marker_times, y1 = rep(marker_height + 0.5, length(marker_times)),
#          col = "blue", lwd = 0.5)
## true firing rate
# obj = simulate.bump(V = 1000, lam_star = 0.2, bump = c(300, 700))
# true_fr = obj$lam_s*10
# plot((1:1000)/50, true_fr, 'l')

##generate simulation data with intercept
# simulate.BP.inter <- function(n,p, gam, sd, lam_star, bump,beta0=0){
#   c = matrix(0,p,n)
#   y = matrix(0,p,n)
#   s = matrix(0,p,n)
#   #beta0 = 1
#   eps = matrix(rnorm(p*n,0,sd),p,n)
#   cp = NULL
#   for(p0 in 1:p){
#     true_cp = simulate.bump(V = n, lam_star = lam_star,  bump = bump)$eps_set
#     s[p0, true_cp] = 1
#     s[p0, 1] = 0
#     cp[[p0]] = true_cp
#     for(i in 1:n){
#       if (i > 1) c[p0,i] = gam * c[p0,(i-1)] + s[p0,i]
#       else c[p0,i] = c[p0,i] + s[p0,i]
#     }
#   }
#   y = c + eps + beta0
#   true_fr = simulate.bump(V = n, lam_star = lam_star,  bump = bump)$lam_s
#   return(list(y = y, c = c, true_cp = cp, true_st = s,true_fr=true_fr))
# }



simulate.bump2 <- function(V, lam_star, bump = bump, pow = pow){
  lam_s = rep(0.05, V)
  lam_s = lam_s + 0.95*get.bump1(1:V, bump[1])*pow + 0.95*get.bump1(1:V, bump[2])*pow 
  J = rpois(1, V*lam_star)
  sj_set = sort(c(sample.int(V, J, replace = FALSE)))
  eps_set = NULL
  for(j in 1:J){
    rj = runif(1)
    if(rj < lam_s[sj_set[j]]) eps_set = c(eps_set, sj_set[j])
  }
  return(list(eps_set = eps_set, sj_set = sj_set))
}

####generate simulation data under dynamic firing rate
simulate.BP2 <- function(n,p, gam, sd, lam_star, bump){
  c = matrix(0,p,n)
  y = matrix(0,p,n)
  s = matrix(0,p,n)
  eps = matrix(rnorm(p*n,0,sd),p,n)
  cp = NULL

  pow_vec = exp(-(1:p-25)^2/1000)
  for(p0 in 1:p){
    pow = pow_vec[p0]
    true_cp = simulate.bump2(V = n, lam_star = lam_star,  bump = bump, pow = pow)$eps_set
    s[p0, true_cp] = 1
    s[p0, 1] = 0
    cp[[p0]] = true_cp
    for(i in 1:n){
      if (i > 1) c[p0,i] = gam * c[p0,(i-1)] + s[p0,i]
      else c[p0,i] = c[p0,i] + s[p0,i]
    }
  }
  y = c + eps
  return(list(y = y, c = c, true_cp = cp, true_st = s))
}








