CPR_GLS <- function(formula,
                    df,
                    VCV_sp,
                    comm,
                    optim.lambda=T,
                    original.VCV=T,
                    ...) {

  require(nlme)

  formula <- as.formula(formula)
  m_without_phylo <- gls(formula,data=df)

  lambda_GLS <- NA
  AIC_optim <- NA
  optimized_model_result <- NA

  C.lambda.GLS <- get_comm_pair_r(comm,VCV_sp)
  m_original_VCV <- gls(formula,
                        data=df,
                        correlation=corSymm(C.lambda.GLS[lower.tri(C.lambda.GLS)], fixed = T))
  AIC_original_VCV <- AIC(m_original_VCV)

  VCV_sp_lambda0 <- VCV_sp*0
  diag(VCV_sp_lambda0) <- diag(VCV_sp)

  C.lambda0.GLS <- get_comm_pair_r(comm,VCV_sp_lambda0)
  m_lambda0_VCV <- gls(formula,
                        data=df,
                        correlation=corSymm(C.lambda0.GLS[lower.tri(C.lambda0.GLS)], fixed = T))
  AIC_star <- AIC(m_lambda0_VCV)

  AIC_optim <-NA
  if (optim.lambda == T) {

    grid_result <- gridSearch(fun=likelihood.lambda,
                              levels=list(lambda=c(0.2,0.4,0.6,0.8)),
                              lower=0,
                              upper=1,
                              formula = formula,
                              df = df,
                              VCV_sp = VCV_sp,
                              comm = comm,
                              printDetail = F,
                              ...)

    ML.opt<-optim(grid_result$minlevels,
                  #runif(1,0.2,0.8),
                  likelihood.lambda,
                  formula=formula,
                  df=df,
                  VCV_sp=VCV_sp,
                  comm=comm,
                  method="L-BFGS-B",
                  lower=0.0,
                  upper=1,
                  ...)

    lambda_GLS<-ML.opt$par
    logL<--ML.opt$value
    VCV_sp_optim <- VCV_sp*lambda_GLS
    diag(VCV_sp_optim) <- diag(VCV_sp)
    C.lambda.GLS<- get_comm_pair_r(comm,VCV_sp_optim)

  m_optim <- gls(formula,
                 data=df,
                 correlation=corSymm(C.lambda.GLS[lower.tri(C.lambda.GLS)], fixed = T))

  AIC_optim <- AIC(m_optim)

  if (AIC_star < AIC_optim & AIC_star < AIC_original_VCV) {
    m_optim <- m_lambda0_VCV
    lambda_GLS = 0
  }

  if (AIC_original_VCV < AIC_optim & AIC_star > AIC_original_VCV) {
    m_optim <- m_original
    lambda_GLS = 1
  }

  optimized_model_result <- summary(m_optim)$tTable
  optimized_model_result[,4] <- 2*pt(abs(summary(m_optim)$tTable[1:2,3]),m_optim$dims$N-m_optim$dims$p-1,lower.tail=F)
  }

  output <- list(optimized_model = optimized_model_result,
                 without_phylo_model = summary(m_without_phylo)$tTable,
                 original_VCV_model = summary(m_original_VCV)$tTable,
                 star_model = summary(m_lambda0_VCV)$tTable,
                 AIC = c(AIC_without_phylo = AIC(m_without_phylo),
                         AIC_original_VCV = AIC(m_original_VCV),
                         AIC_optim = AIC_optim,
                         AIC_star = AIC_star),
                 optimized_lambda = lambda_GLS)
  return(output)
}
