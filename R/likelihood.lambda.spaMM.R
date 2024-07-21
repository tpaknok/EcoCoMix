likelihood.lambda.spaMM<-function(lambda,
                                  formula,
                                  data,
                                  VCV_sp,
                                  comm,
                                  init=list(lambda=NaN,phi=NaN),
                                  method.spaMM = "REML",
                                  ...){

  message(lambda)

  formula <- as.formula(formula)
  VCV_sp_lambda <- VCV_sp*lambda
  diag(VCV_sp_lambda) <- diag(VCV_sp)

  C.lambda<- get_comm_pair_r(comm,VCV_sp_lambda)
  rownames(C.lambda) <- 1:nrow(data)

  n <- ncol(C.lambda)

  m_spamm <- fitme(formula,
              corrMatrix=as_precision(C.lambda),
              data=data,
              init=init,
              method=method.spaMM,
              ...)

  AIC <- AIC(m_spamm,verbose=F)[[1]]
  return(AIC)
}

