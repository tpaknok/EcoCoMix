CPR_spaMM <- function(formula,
                      data,
                      VCV_sp,
                      comm,
                      force.PD = FALSE,
                      optim.lambda=TRUE,
                      original.VCV=TRUE,
                      AIC_threshold = -4,
                      method.spaMM = "ML",
                      true_VCV = NULL,
                      control.optim=NULL,
                      comm_kronecker = NULL,
                      int_model = TRUE,
                      ...) {

  require(nlme)
  require(NMOF)
  require(spaMM)
  require(admisc)

  time1<-Sys.time()
  data$comp_id <- 1:nrow(data)

  formula_text <- Reduce(paste,deparse(formula))
  formula_elements <- unlist(strsplit(formula_text,split="\\+"))
  formula_elements <- formula_elements[!grepl("comp_id",formula_elements)]
  no_phylo_formula <- as.formula(paste(formula_elements,collapse="+"))

  formula <- as.formula(formula)
  no_phylo_formula <- as.formula(no_phylo_formula)


  m_without_comp <- fitme(no_phylo_formula,
                          data=data,
                          method=method.spaMM,
                          ...)


  AIC_without_comp <- AIC(m_without_comp,verbose=F)[[2]]

  lambda_spaMM <- NA
  AIC_optim <- AIC_true <- AIC_optim_int <- AIC_int_only <- NA
  optimized_model_result <- NA
  best_m <- m_optim <- NA
  lambda_spaMM_int <- lambda_spaMM <- NA
  comm_cov <- get_comm_pair_r(comm,
                              VCV_sp,
                              comm_kronecker=comm_kronecker,
                              force.PD=force.PD)


  C.lambda.spaMM <- comm_cov$corM
  comm_kronecker <- comm_cov$comm_kronecker

  rownames(C.lambda.spaMM) <- as.character(data$comp_id)

  m_original_VCV <- fitme(formula,
                          corrMatrix=as_precision(C.lambda.spaMM),
                          data=data,
                          method=method.spaMM,
                          ...)


  .drop1_spamm <- function(model,C) {
    conv <- 0

    if (length(fixef(model)) == 1 & names(fixef(model))[[1]] == ("(Intercept)")) {
      msg <- "drop1 not conducted"
      result_satt <- "Intercept-only model. drop1 was not conducted"
    } else {
      model$call$corrMatrix <- as_precision(C)
      model$call$method <- ifelse(method.spaMM == "REML","REML","ML")
      model$call$data <- model$data
      model$family$family <- as.character(paste0(model$family$family))
      model$family$link <- as.character(paste0(model$family$link))

      result_satt <- drop1(model,verbose=F)
      msg <- tryCatchWEM(drop1(model,verbose=F),capture=F)

      if (is.null(msg$warning)) {
        msg$warning <- "OK"
      }

      if (!grepl("converge",msg$warning)) {
        conv <- 1
      }
    }

      result <- list(result=result_satt,
                     msg=msg,
                     conv=conv)
    return(result)
  }

  original_VCV_model_satt <- .drop1_spamm(m_original_VCV,C.lambda.spaMM)

  AIC_original_VCV <- AIC(m_original_VCV,verbose=F)[[2]]
  VCV_sp_lambda0 <- VCV_sp*0
  diag(VCV_sp_lambda0) <- diag(VCV_sp)

  C.lambda0.spaMM <- get_comm_pair_r(comm,
                                     VCV_sp_lambda0,
                                     comm_kronecker=comm_kronecker,
                                     force.PD=force.PD)$corM

  rownames(C.lambda0.spaMM) <- as.character(data$comp_id)

  m_lambda0 <- fitme(formula,
                   corrMatrix=as_precision(C.lambda0.spaMM),
                   data=data,
                   method=method.spaMM,
                   ...)

  AIC_star <- AIC(m_lambda0,verbose=F)[[2]]

  re <- names(ranef(m_lambda0))
  re <- paste0(re,collapse="+")
  response <- as.character(m_lambda0$predictor)[[2]]
  f_null <- as.formula(paste0(response,"~1+",re))

  C.true <- true_model_satt <- m_true<- NULL
  if (!is.null(true_VCV)) {
    C.true <- get_comm_pair_r(comm,
                              true_VCV,
                              comm_kronecker=comm_kronecker,
                              force.PD=force.PD)$corM
    rownames(C.true) <- as.character(data$comp_id)

    m_true <- fitme(formula,
                    corrMatrix=as_precision(C.true),
                    data=data,
                    method=method.spaMM,
                    ...)
    AIC_true <- AIC(m_true,verbose=F)[[2]]
    true_model_satt <- .drop1_spamm(m_true,C.true)
  }

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
                              method.spaMM=method.spaMM,
                              comm_kronecker=comm_kronecker,
                              ...)

    ML.opt<-optim(grid_result$minlevels,
                  #runif(1),
                  likelihood.lambda.spaMM,
                  formula=formula,
                  data=data,
                  VCV_sp=VCV_sp,
                  method = "L-BFGS-B",
                  comm=comm,
                  lower=0.0,
                  upper=1,
                  method.spaMM = method.spaMM,
                  comm_kronecker=comm_kronecker,
                  control=control.optim,
                  ...)


    lambda_spaMM <- ML.opt$par
    logL<--ML.opt$value
    VCV_sp_optim <- VCV_sp*lambda_spaMM
    diag(VCV_sp_optim) <- diag(VCV_sp)
    C.lambda.optim.spaMM <- get_comm_pair_r(comm,
                                            VCV_sp_optim,
                                            comm_kronecker = comm_kronecker,
                                            force.PD=force.PD)$corM
    rownames(C.lambda.optim.spaMM) <- as.character(data$comp_id)

    m_optim <- fitme(formula,
                     corrMatrix=as_precision(C.lambda.optim.spaMM),
                     data=data,
                     method=method.spaMM,
                     ...)

    AIC_optim <- AIC(m_optim,verbose=F)[[2]]

    if (AIC_star < AIC_optim & AIC_star < AIC_original_VCV) {
      m_optim <- m_lambda0
      lambda_spaMM <- 0
      AIC_optim <- AIC_star
    }

    if (AIC_original_VCV <= AIC_optim & AIC_original_VCV <= AIC_star) {
      m_optim <- m_original_VCV
      lambda_spaMM <- 1
      AIC_optim <- AIC_original_VCV
    }

    if (AIC_optim - AIC_without_comp < AIC_threshold) {
      best_m <- m_optim
    } else {
      best_m <- m_without_comp
    }

    optim_model_satt <- .drop1_spamm(m_optim,get(gsub("as_precision(*)","",m_optim$call$corrMatrix)[[2]]))

    if (!is.null(best_m$call$corrMatrix)) {
      best_model_satt <- .drop1_spamm(best_m,get(gsub("as_precision(*)","",best_m$call$corrMatrix)[[2]]))
    } else {
      best_model_satt <- list(result=anova(m_without_comp))
    }

    m_optim_int <- NULL
    if (int_model & optim.lambda) {
      grid_result_int <- gridSearch(fun=likelihood.lambda.spaMM,
                                    levels=list(lambda=c(0.2,0.4,0.6,0.8)),
                                    lower=0,
                                    upper=1,
                                    formula = f_null,
                                    data = data,
                                    VCV_sp = VCV_sp,
                                    comm = comm,
                                    printDetail = F,
                                    method.spaMM=method.spaMM,
                                    comm_kronecker=comm_kronecker,
                                    ...)

      ML.opt_int<-optim(grid_result_int$minlevels,
                        #runif(1),
                        likelihood.lambda.spaMM,
                        formula=f_null,
                        data=data,
                        VCV_sp=VCV_sp,
                        method = "L-BFGS-B",
                        comm=comm,
                        lower=0.0,
                        upper=1,
                        method.spaMM = method.spaMM,
                        comm_kronecker=comm_kronecker,
                        control=control.optim,
                        ...)

      lambda_spaMM_int <- ML.opt_int$par
      VCV_sp_optim_int <- VCV_sp*lambda_spaMM_int
      diag(VCV_sp_optim_int) <- diag(VCV_sp)
      C.lambda.optim.spaMM.int <- get_comm_pair_r(comm,
                                                  VCV_sp_optim_int,
                                                  comm_kronecker = comm_kronecker,
                                                  force.PD=force.PD)$corM
      rownames(C.lambda.optim.spaMM.int) <- as.character(data$comp_id)

      m_optim_int <- fitme(f_null,
                           corrMatrix=as_precision(C.lambda.optim.spaMM.int),
                           data=data,
                           method=method.spaMM,
                           ...)

      m_optim_int_star <- fitme(f_null,
                           corrMatrix=as_precision(C.lambda0.spaMM),
                           data=data,
                           method=method.spaMM,
                           ...)

      m_optim_int_BM <- fitme(f_null,
                                corrMatrix=as_precision(C.lambda.spaMM),
                                data=data,
                                method=method.spaMM,
                                ...)

      m_int <- fitme(as.formula(paste0(response,"~1")),
                     data=data,
                     ...)

      AIC_int <- c(AIC(m_optim_int,verbose=F)[[2]], AIC(m_optim_int_star,verbose=F)[[2]],AIC(m_optim_int_BM,verbose=F)[[2]])
      AIC_int_pos <- which.min(AIC_int)
      AIC_optim_int <- AIC_int[AIC_int_pos]

      if (AIC_int_pos == 2) {
        m_optim_int <- m_optim_int_star
        optimized_lambda_int <- 0
        }
      if (AIC_int_pos == 3) {
        m_optim_int <- m_optim_int_BM
        optimized_lambda_int <- 1
      }

      AIC_int_only <- AIC(m_int,verbose=F)[[2]]

    }
  }

  ##########drop 1
   # best_model_satt <- drop1(best_m)
  # optim_model_satt <- drop1(m_optim)
  # original_VCV_model_satt <- drop1(m_original_VCV)
  # true_model_satt <- drop1(m_true)

  if (optim.lambda) {
    conv_best_m <- best_model_satt$conv
    conv_optim_m <- optim_model_satt$conv
    best_model_satt <- best_model_satt$result
    optim_model_satt <- optim_model_satt$result
    optimized_lambda <- lambda_spaMM
    optimized_lambda_int <- lambda_spaMM_int
  } else {
    best_model_satt <- optim_model_satt <- optimized_lambda <- optimized_lambda_int <- NA
    conv_best_m <- conv_optim_m <- 1
  }

  output <- list(best_model = best_m,
                 best_model_satt = best_model_satt,
                 optimized_lambda_model = m_optim,
                 optimized_lambda_model_satt = optim_model_satt,
                 original_VCV_model = m_original_VCV,
                 original_VCV_m_satt = original_VCV_model_satt$result,
                 true_model = m_true,
                 true_model_satt = true_model_satt$result,
                 without_comp_anova = anova(m_without_comp),
                 without_comp_model = m_without_comp,
                 intercept_only_model = m_optim_int,
                 star_model = m_lambda0,
                 AIC = c(AIC_without_comp =  AIC_without_comp,
                         AIC_original_VCV =  AIC_original_VCV,
                         AIC_optim = AIC_optim,
                         AIC_star = AIC_star,
                         AIC_true = AIC_true,
                         AIC_optim_int = AIC_optim_int,
                         AIC_int_only = AIC_int_only),
                 optimized_lambda = optimized_lambda ,
                 optimized_lambda_int = optimized_lambda_int,
                 conv = c(best_m=conv_best_m,
                          optim_m=conv_optim_m,
                          orig_m=original_VCV_model_satt$conv,
                          true_m=true_model_satt$conv),
                 min_richness=min(rowSums(comm)),
                 max_richness=max(rowSums(comm)),
                 nspp = ncol(comm))

}
