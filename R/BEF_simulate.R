BEF_simulate <- function(comm,
                         VCV_sp = NULL,
                         nspp=15,
                         nsite=50,
                         min_richness=1,
                         max_richness=4,
                         spaMM_formula,
                         b0=0,
                         b1=0,
                         signals_X="phy_cor",
                         noise_mean = 0,
                         noise_sd = 1,
                         non_phy_cor_mean = 0,
                         non_phy_cor_sd = 1,
                         lambda_true=1,
                         scale_all=FALSE,
                         conv_fail_drop = TRUE,
                         ...) {
                         #seed=1000) {
    library(MASS)
    library(vegan)
    library(CPR)
    library(ape)
    library(phytools)

    conv <- c(0,0,0,0)
    count <- 0
    dummy <- 1

    while (sum(grepl("0",conv)) > 0) {
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
    } else {
      nspp <- ncol(comm)
      nsite <- nrow(comm)
      min_richness <- min(rowSums(comm))
      max_richness <- max(rowSums(comm))
    }

    if(is.null(VCV_sp)){
      tree <- pbtree(n=ncol(comm))
      tree$tip.label <- colnames(comm)

      vcv <- vcv(tree)
      vcv_true <- vcv * lambda_true
      diag(vcv_true) <- diag(vcv)
    } else {

      vcv <- VCV_sp
      vcv_true <- VCV_sp * lambda_true
      diag(vcv_true) <- diag(VCV_sp)
      }

    comm_pair <- get_comm_pair_r(comm,vcv_true,force.PD=F)
    C_true <- comm_pair$corM
    comm_kronecker <- comm_pair$comm_kronecker

    if (signals_X == "phy_cor") x1 <- mvrnorm(1,rep(0,nrow(C_true)),C_true)
    if (signals_X == "non_phy_cor") x1 <- rnorm(nrow(comm),non_phy_cor_mean,non_phy_cor_sd)
    if (signals_X == "sr") x1 <- rowSums(comm)

    x2 <- mvrnorm(1,rep(0,nrow(C_true)),C_true)

    x3 <- rnorm(nrow(comm),noise_mean,noise_sd)

    data <- data.frame(x1=x1,x2=x2,x3=x3,comp_id=as.character(1:nrow(comm)))

    if (scale_all) data[,1:3] <- scale(data[,1:3])

    y <- b0+b1*x1+x2+x3

    data <- data.frame(y=y,data)

    sim_data <- list(data=data,
                     sim_phy=vcv,
                     true_phy=vcv_true)

    spaMM_formula <- as.formula(spaMM_formula)

    models <- CPR_spaMM(spaMM_formula,
                        data=sim_data$data,
                        VCV_sp = sim_data$sim_phy,
                        true_VCV = sim_data$true_phy,
                        comm=comm,
                        comm_kronecker = comm_kronecker,
                        ...)

    conv <- models$conv

    count=count+1

    if (conv_fail_drop == F) conv <- c(1,1,1,1)

    }


    # result <- list(models=models,
                   # true_lambda = lambda_true,
                   # count = count)

    #extract significance

    m_optim_sig <- NA

    if (all(!is.na(unlist(models$optimized_lambda_model_satt))))
    m_optim_sig <- ifelse(models$optimized_lambda_model_satt$`Pr(>F)` < 0.05,1,0)

    m_true_sig <- ifelse(models$true_model_satt$`Pr(>F)` < 0.05,1,0)
    m_original_sig <- ifelse(models$original_VCV_m_satt$`Pr(>F)` < 0.05,1,0)

    m_best_sig <- NA

    if (all(!is.na(unlist(models$best_model_satt))))
    m_best_sig <- ifelse(models$best_model_satt$`Pr(>F)`[[1]] < 0.05,1,0)

    m_without_comp_sig <- ifelse(anova(models$without_comp_model)$`Pr(>F)`[[1]] < 0.05,1,0)

    sig <- c(m_optim_sig=m_optim_sig,
             m_true_sig=m_true_sig,
             m_original_sig=m_original_sig,
             m_best_sig=m_best_sig,
             m_without_comp_sig=m_without_comp_sig)

    slope <- c(m_optim_slope=ifelse(!is.na(m_optim_sig),summary(models$optimized_lambda_model,verbose=F)$beta_table[2,1],NA),
               m_true_slope = summary(models$true_model,verbose=F)$beta_table[2,1],
               m_original_slope = summary(models$original_VCV_model,verbose=F)$beta_table[2,1],
               m_best_slope = ifelse(!is.na(m_best_sig),summary(models$best_model,verbose=F)$beta_table[2,1],NA),
               m_without_comp_slope = summary(models$without_comp_model,verbose=F)$beta_table[2,1])

    optim_lambda <- ifelse(!is.na(m_optim_sig),c(optim_lambda=models$optimized_lambda),NA)

    AIC <- c(models$AIC)

    result <- c(sig,
                slope,
                optim_lambda = optim_lambda,
                optim_lambda_int = models$optimized_lambda_int,
                AIC,
                min_richness=models$min_richness,
                max_richness=models$max_richness,
                nspp=models$nspp,
                true_lambda = lambda_true,
                r_x1x2 = cor(x1,x2),
                b1 = b1,
                count=count,
                optim_r2m = ifelse(!is.na(m_optim_sig),get_R2(models$optimized_lambda_model)[[1]],NA),
                optim_r2c = ifelse(!is.na(m_optim_sig),get_R2(models$optimized_lambda_model)[[2]],NA),
                NumDF = ifelse(!is.na(m_optim_sig),models$optimized_lambda_model_satt$NumDF,NA),
                DenDF = ifelse(!is.na(m_optim_sig),models$optimized_lambda_model_satt$DenDF,NA)
                )
    result <- data.frame(t(result))

    result$signals_X <- NA
    result$phylo_structure <- NA
    result$signals_X <- signals_X
    result$phylo_structure <- ifelse(is.null(VCV_sp),"Simulated","Provided")

    return(result)
  # }
}


