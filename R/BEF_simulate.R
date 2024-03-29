BEF_simulate <- function(comm,
                         V,
                         ef_mean,
                         sd,
                         b1=0,
                         signals_X="sr",
                         signals_intercept=T,
                         signals_slope=T,
                         lambda_true=1,
                         sim=500) {

  C <- get_comm_pair_r(comm,V)

  V_true <- V*lambda_true
  diag(V_true) <- diag(V)
  C_true <- get_comm_pair_r(comm,V_true)

  if (signals_X == "sr") {
    x <- replicate(rowSums(comm),n=sim) #species richness
  }

  if (signals_X == "phy_cor") {
    x <- replicate(t(chol(C_true)) %*% rnorm(n=nrow(comm),mean=ef_mean,sd=sd),n=sim) #signals in x
  }

  if (signals_X == "no_phy_cor") {
    x <- replicate(rnorm(nrow(comm),0,1),n=sim) #no signals in X
  }

  if (signals_intercept == T & signals_slope==T) {
    library(MASS)
    RE <- replicate(mvrnorm(length(x),
                            mu=c(0,0),
                            Sigma=rbind(c(1.0, 0),c(0, 1.0))),n = sim)#simulate random effect
    comm_int <- t(chol(C_true))%*%RE[,1]
    comm_beta <- t(chol(C_true))%*%RE[,2]
    sr_E_cor <- cor(x,comm_beta*x+comm_int)
    slope_int_cor <- cor(comm_beta,comm_int)
    y <- b1*x+comm_beta*x+comm_int #signals in intercept and slope
  }

  if (signals_intercept == T & signals_slope == F) {
    comm_int <- t(chol(C_true)) %*% replicate(rnorm(n=nrow(comm),mean=ef_mean,sd=sd),n=sim)
    y <- b1*x+comm_int #signals in intercept but not slope
    sr_E_cor <- NA
    slope_int_cor <- NA
  }

  if (signals_intercept==F & signals_slope == F) {
    y <- b1*x+replicate(rnorm(nrow(comm)),n=sim)#no signals in error
    sr_E_cor <- NA
    slope_int_cor <- NA
  }

  sim_all <- lapply(1:sim, function(z) data.frame(y=y[,z],x=x[,z]))
  sim_dat <- list(sim_dat = sim_all, true_phy = C_true)

  return(sim_dat)
}
