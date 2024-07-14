BEF_simulate <- function(comm,
                         VCV_sp = NULL,
                         scale=1,
                         mean,
                         b1=0,
                         signals_X="sr",
                         # signals_intercept=T,
                         # signals_slope=F,
                         lambda_true=1,
                         sim=500,
                         seed=1000) {

  library(MASS)

  set.seed(seed)

  if (is.null(VCV_sp)) {
    pbtree <- pbtree(n=ncol(comm),scale=1,tip.label=colnames(comm),nsim=sim)
    VCV_sp <- lapply(1:sim,function(x) vcv(pbtree[[x]]))
  } else {
    VCV_sp <- lapply(1:sim,function(x) VCV_sp)
  }

  C <- lapply(1:sim,function(x) get_comm_pair_r(comm,VCV_sp[[x]],force.PD=F))

  if (signals_X == "sr") {
    x <- replicate(rowSums(comm),n=sim) #species richness
  }

  # if (signals_X == "phy_cor") {
  #   x <- replicate(t(chol(C_true)) %*% rnorm(n=nrow(comm),mean=ef_mean,sd=sd),n=sim) #signals in x
  # }

  # if (signals_X == "phy_cor") {
  #   x <- replicate(mvrnorm(1,rep(0,nrow(C),C)),n=sim) #signals in x
  # }
  #

  if (signals_X == "no_phy_cor") {
     x <- replicate(rnorm(nrow(comm),0,1),n=sim) #no signals in X
  }

  lambda_true <- matrix(lambda_true,nrow(VCV_sp[[1]]),ncol(VCV_sp[[1]]))
  diag(lambda_true) <- 1
  VCV_true <- lapply(1:sim,function(x) VCV_sp[[x]] * lambda_true)
  C_true <- lapply(1:sim, function(x) get_comm_pair_r(comm,VCV_true))

  # if (signals_intercept == T & signals_slope==T) {
  #   library(MASS)
  #   RE <- replicate(mvrnorm(length(x),
  #                           mu=c(0,0),
  #                           Sigma=rbind(c(1.0, 0),c(0, 1.0))),n = sim)#simulate random effect
  #   comm_int <- t(chol(C_true))%*%RE[,1]
  #   comm_beta <- t(chol(C_true))%*%RE[,2]
  #   sr_E_cor <- cor(x,comm_beta*x+comm_int)
  #   slope_int_cor <- cor(comm_beta,comm_int)
  #   y <- b1*x+comm_beta*x+comm_int#signals in intercept and slope
  # }
  #
  # if (signals_intercept == T & signals_slope == F) {
  #   comm_int <- t(chol(C_true)) %*% replicate(rnorm(n=nrow(comm),mean=ef_mean,sd=sd),n=sim)
  #   y <- b1*x+comm_int#signals in intercept but not slope
  #   sr_E_cor <- NA
  #   slope_int_cor <- NA
  # }
  #
  # if (signals_intercept==F & signals_slope == F) {
  #   y <- b1*x+replicate(rnorm(nrow(comm)),n=sim)#no signals in error
  #   sr_E_cor <- NA
  #   slope_int_cor <- NA
  # }

  y <- lapply(1:sim, function(x) 0+b1*x+mvrnorm(1,rep(0,nrow(C_true[[x]])),C_true[[x]]))

  sim_all <- lapply(1:sim, function(z) data.frame(y=y[[z]],x=x[,z]))
  sim_C <- lapply(1:sim, function(z) data.frame)
  sim_dat <- list(sim_dat = sim_all, true_phy = VCV_true)

  return(sim_dat)
}
