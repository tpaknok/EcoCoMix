CPR_spaMM <- function(formula,
                    data,
                    VCV_sp,
                    comm,
                    optim.lambda=T,
                    original.VCV=T,
                    AIC_threshold = -4,
                    init=list(lambda=NaN,phi=NaN),
                    ...) {

  require(nlme)
  require(NMOF)
  require(spaMM)

  message("modelling")
  data$comp_id <- 1:nrow(data)

  formula_text <- Reduce(paste,deparse(formula))
  formula_elements <- unlist(strsplit(formula_text,split="\\+"))
  formula_elements <- formula_elements[!grepl("comp_id",formula_elements)]
  no_phylo_formula <- as.formula(paste(formula_elements,collapse="+"))

  formula <- as.formula(formula)
  no_phylo_formula <- as.formula(no_phylo_formula)

  m_without_comp <- fitme(no_phylo_formula,data=data)
  AIC_without_phylo <- AIC(m_without_comp,verbose=F)[[1]]

  lambda_spaMM <- NA
  AIC_optim <- NA
  optimized_model_result <- NA
  best_m <- m_optim <- NA

  C.lambda.spaMM <- get_comm_pair_r(comm,VCV_sp)
  rownames(C.lambda.spaMM) <- data$comp_id
  m_original_VCV <- fitme(formula,
                          corrMatrix=C.lambda.spaMM,
                          data=data,
                          init=init,
                          ...)

  AIC_original_VCV <- AIC(m_original_VCV,verbose=F)[[1]]
  VCV_sp_lambda0 <- VCV_sp*0
  diag(VCV_sp_lambda0) <- diag(VCV_sp)

  C.lambda0.spaMM <- get_comm_pair_r(comm,VCV_sp_lambda0)
  rownames(C.lambda0.spaMM) <- data$comp_id
  m_lambda0 <- fitme(formula,
                   corrMatrix=C.lambda0.spaMM,
                   data=data,
                   init=init,
                   ...)
  AIC_star <- AIC(m_lambda0,verbose=F)[[1]]

  if (optim.lambda == T) {

    grid_result <- gridSearch(fun=likelihood.lambda.spaMM,
                              levels=list(lambda=c(0.2,0.4,0.6,0.8)),
                              lower=0,
                              upper=1,
                              formula = formula,
                              data = data,
                              VCV_sp = VCV_sp,
                              comm = comm,
                              printDetail = F,
                              init=init,
                              ...)

    ML.opt<-optim(grid_result$minlevels,
                  likelihood.lambda.spaMM,
                  formula=formula,
                  data=data,
                  VCV_sp=VCV_sp,
                  method = "L-BFGS-B",
                  comm=comm,
                  lower=0.0,
                  upper=1,
                  init=init,
                  ...)

    lambda_spaMM<-ML.opt$par
    logL<--ML.opt$value
    VCV_sp_optim <- VCV_sp*lambda_spaMM
    diag(VCV_sp_optim) <- diag(VCV_sp)
    C.lambda.spaMM<- get_comm_pair_r(comm,VCV_sp_optim)
    rownames(C.lambda.spaMM) <- data$comp_id

    m_optim <- fitme(formula,
                     corrMatrix=C.lambda.spaMM,
                     data=data,
                     ...)

    AIC_optim <- AIC(m_optim,verbose=F)[[1]]

    if (AIC_star < AIC_optim & AIC_star < AIC_original_VCV) {
      m_optim <- m_lambda0
      lambda_spaMM <- 0
      AIC_optim <- AIC_star
    }

    if (AIC_original_VCV < AIC_optim & AIC_star > AIC_original_VCV) {
      m_optim <- m_original_VCV
      lambda_spaMM = 1
      AIC_optim <- AIC_original_VCV
    }

    if (AIC_optim - AIC_without_phylo < AIC_threshold) {
      best_m <- m_optim
    } else {
      best_m <- m_without_comp
    }
  }

  output <- list(best_m = best_m,
                 optimized_lambda_model = m_optim,
                 without_comp_model = m_without_comp,
                 original_VCV_model = m_original_VCV,
                 star_model = m_lambda0,
                 AIC = c(AIC_without_phylo =  AIC_without_phylo,
                         AIC_original_VCV =  AIC_original_VCV,
                         AIC_optim = AIC_optim,
                         AIC_star = AIC_star),
                 optimized_lambda = lambda_spaMM)
}
