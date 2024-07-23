BEF_simulate <- function(comm,
                         VCV_sp = NULL,
                         scale=1,
                         b1=0,
                         signals_X="phy_cor",
                         signals_Y = T,
                         intercept = 0,
                         y_mean = 0,
                         y_sd = 1,
                         x_mean=0,
                         x_sd = 1,
                         noise_mean = 0,
                         noise_sd = 0,
                         lambda_true=1,
                         sim=500,
                         seed=1000) {

  library(MASS)
  library(vegan)
  set.seed(seed)

  sim_all <- VCV_sp_list <- VCV_true_list <- list()

  for (i in 1:sim) {
    if(is.null(VCV_sp)){
      tree <- pbtree(n=ncol(comm))
      tree$tip.label <- colnames(comm)

      vcv <- vcv(tree)
      VCV_sp_list[[i]] <- vcv
      vcv_true <- vcv*lambda_true
      diag(vcv_true) <- diag(vcv)
      } else {
      VCV_sp_list[[i]] <- VCV_sp
      vcv_true <- VCV_sp*lambda_true
      diag(vcv_true) <- diag(VCV_sp)
    }


    VCV_true_list[[i]] <- vcv_true
    C_true <- get_comm_pair_r(comm,vcv_true)

    if (signals_X == "phy_cor") {
    #x <- t(chol(C_true)) %*% rnorm(nrow(comm),x_mean,x_sd)
      x <- mvrnorm(1,rep(0,nrow(C_true)),C_true)
    }

    if (signals_X == "sr") {
      x <- specnumber(comm)
    }

    if (signals_X == "no_phy_cor") {
      x <- rnorm(nrow(comm),x_mean,x_sd)
    }

    if (signals_Y == T) {
    #y <- t(chol(C_true)) %*% rnorm(nrow(comm),y_mean,y_sd)
      y <- intercept+b1*x+mvrnorm(1,rep(0,nrow(C_true)),C_true)+rnorm(nrow(comm),noise_mean,noise_sd)
    } else {
      y <- intercept+b1*x+rnorm(nrow(comm),y_mean,y_sd)+rnorm(nrow(comm),noise_mean,noise_sd)
    }

    sim_all[[i]] <- data.frame(y=y,x=x)
  }


  sim_dat <- list(sim_dat = sim_all,
                  sim_phy = VCV_sp_list,
                  true_phy = VCV_true_list)

  return(sim_dat)
}
