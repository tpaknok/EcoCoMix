BEF_simulate <- function(comm,
                         VCV_sp = NULL,
                         scale=1,
                         nspp=15,
                         nsite=50,
                         min_richness=1,
                         max_richness=4,
                         spaMM_formula,
                         b1=0,
                         signals_X="phy_cor",
                         signals_Y = T,
                         intercept = 0,
                         y_mean = 0,
                         y_sd = 1,
                         x_mean=0,
                         x_sd = 1,
                         noise_mean = 0,
                         noise_sd = 1,
                         lambda_true=1) {
                         #seed=1000) {

  # library(MASS)
  # library(vegan)
  # set.seed(seed)
  #
  # sim_all <- VCV_sp_list <- VCV_true_list <- list()
  #
  # for (i in 1:sim) {
  #   if(is.null(VCV_sp)){
  #     tree <- pbtree(n=ncol(comm))
  #     tree$tip.label <- colnames(comm)
  #
  #     vcv <- vcv(tree)
  #     VCV_sp_list[[i]] <- vcv
  #     vcv_true <- vcv*lambda_true
  #     diag(vcv_true) <- diag(vcv)
  #     } else {
  #     VCV_sp_list[[i]] <- VCV_sp
  #     vcv_true <- VCV_sp*lambda_true
  #     diag(vcv_true) <- diag(VCV_sp)
  #   }
  #
  #
  #   VCV_true_list[[i]] <- vcv_true
  #   C_true <- get_comm_pair_r(comm,vcv_true)
  #
  #   if (signals_X == "phy_cor") {
  #   #x <- t(chol(C_true)) %*% rnorm(nrow(comm),x_mean,x_sd)
  #     x <- mvrnorm(1,rep(0,nrow(C_true)),C_true)
  #   }
  #
  #   if (signals_X == "sr") {
  #     x <- specnumber(comm)
  #   }
  #
  #   if (signals_X == "no_phy_cor") {
  #     x <- rnorm(nrow(comm),x_mean,x_sd)
  #   }
  #
  #   if (signals_Y == T) {
  #   #y <- t(chol(C_true)) %*% rnorm(nrow(comm),y_mean,y_sd)
  #     y <- intercept+b1*x+mvrnorm(1,rep(0,nrow(C_true)),C_true)+rnorm(nrow(comm),noise_mean,noise_sd)
  #   } else {
  #     y <- intercept+b1*x+rnorm(nrow(comm),y_mean,y_sd)+rnorm(nrow(comm),noise_mean,noise_sd)
  #   }
  #
  #   sim_all[[i]] <- data.frame(y=y,x=x)
  # }
  # require(parallel)
  #
  # .simulate_data <- function(comm,
  #                            VCV_sp,
  #                            scale,
  #                            b1,
  #                            signals_X,
  #                            signals_Y,
  #                            intercept,
  #                            y_mean,
  #                            y_sd,
  #                            x_mean,
  #                            x_sd,
  #                            noise_mean,
  #                            noise_sd,
  #                            lambda_true,
  #                            force.PD = F) {
    library(MASS)
    library(vegan)
    library(CPR)
    library(ape)
    library(phytools)

    conv <- c(0,0,0,0)
    count <- 0

    comp <- list()

    if(is.null(comm)) {
    comm <- matrix(0,nrow=nspp,ncol=nsite)

    while (any(colSums(comm) == 0) | max(rowSums(comm)) < max_richness) {
      for (j in 1:(nsite-2)) {
        pos <- sample(1:nspp, sample(min_richness:max_richness,1))
        temp_comp <- rep(0,nspp)
        temp_comp[pos] <- 1
        comp[[j]] <- temp_comp
      }
      comm <- do.call(rbind,comp)
      comm <- rbind(comm,c(rep(1,max_richness),rep(0,nspp-max_richness)))
      comm <- rbind(comm,c(1,rep(0,nspp-1)))
      colnames(comm) <- paste0("sp",1:nspp)
      }
    }

    while (sum(grepl("0",conv)) > 0) {

    if(is.null(VCV_sp)){
      tree <- pbtree(n=ncol(comm))
      tree$tip.label <- colnames(comm)

      vcv <- vcv(tree)
      vcv_true <- vcv*lambda_true
      diag(vcv_true) <- diag(vcv)
    } else {
      vcv_true <- VCV_sp*lambda_true
      diag(vcv_true) <- diag(VCV_sp)
    }

    C_true <- get_comm_pair_r(comm,vcv_true,force.PD=F)$covM
    x <- mvrnorm(1,rep(0,nrow(C_true)),C_true)
    y <- b1*x+mvrnorm(1,rep(0,nrow(C_true)),C_true)+rnorm(nrow(comm),noise_mean,noise_sd)
    data <- data.frame(y=y,x=x,comp_id=as.character(1:nrow(comm)))

    sim_data <- list(data=data,
                     sim_phy=vcv,
                     true_phy=vcv_true)

    spaMM_formula <- as.formula(spaMM_formula)

    models <- CPR_spaMM(spaMM_formula,
                        data=sim_data$data,
                        VCV_sp = sim_data$sim_phy,
                        true_VCV = sim_data$true_phy,
                        comm=comm,
                        family= "gaussian",
                        optim.lambda = T,
                        control.HLfit = list(max.iter=10000),
                        control.optim=list(factr=1e15))

    conv <- models$conv

    count=count+1
    }

    # result <- list(models=models,
                   # true_lambda = lambda_true,
                   # count = count)

    #extract significance?
    m_optim_sig <- ifelse(models$optimized_lambda_model_satt$`Pr(>F)` < 0.05,1,0)
    m_true_sig <- ifelse(models$true_model_satt$`Pr(>F)` < 0.05,1,0)
    m_original_sig <- ifelse(models$original_VCV_m_satt$`Pr(>F)` < 0.05,1,0)
    m_best_sig <- ifelse(models$best_model_satt$`Pr(>F)` < 0.05,1,0)
    m_without_comp_sig <- ifelse(summary(models$without_comp_model,details=T,verbose=F)$beta_table[2,4] < 0.05,1,0)
    sig <- c(m_optim_sig=m_optim_sig,
             m_true_sig=m_true_sig,
             m_original_sig=m_original_sig,
             m_best_sig=m_best_sig,
             m_without_comp_sig=m_without_comp_sig)

    slope <- c(m_optim_slope=summary(models$optimized_lambda_model,verbose=F)$beta_table[2,1],
               m_true_slope = summary(models$true_model,verbose=F)$beta_table[2,1],
               m_original_slope = summary(models$original_VCV_model,verbose=F)$beta_table[2,1],
               m_best_slope = summary(models$best_model,verbose=F)$beta_table[2,1],
               m_without_comp_slope = summary(models$without_comp_model,verbose=F)$beta_table[2,1])

    optim_lambda <- c(models$optimized_lambda)

    AIC <- c(models$AIC)

    result <- c(sig,slope,
                optim_lambda,
                AIC,
                min_richness=models$min_richness,
                max_richness=models$max_richness,
                nspp=models$nspp,
                true_lambda = lambda_true,
                count=count)

    return(result)
  # }
}


