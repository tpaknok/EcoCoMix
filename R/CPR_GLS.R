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

  if (optim.lambda == T) {
    ML.opt<-optim(runif(1,0.2,0.8),
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

  optimized_model_result <- summary(m_optim)$tTable
  optimized_model_result[,4] <- 2*pt(abs(summary(m_optim)$tTable[1:2,3]),m_optim$dims$N-m_optim$dims$p-1,lower.tail=F)
  AIC_optim <- AIC(m_optim)
  }

  output <- list(without_phylo_model = summary(m_without_phylo)$tTable,
                 optimized_model = optimized_model_result,
                 original_VCV_model = summary(m_original_VCV)$tTable,
                 AIC = c(AIC_without_phylo = AIC(m_without_phylo),
                         AIC_original_VCV = AIC(m_original_VCV),
                         AIC_optim = AIC_optim),
                 optimized_lambda = lambda_GLS)
  return(output)
}
