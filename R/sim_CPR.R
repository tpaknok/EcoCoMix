sim_CPR <- function(comm,V,ef_mean,sd,
                     b1=0,
                     signals_X="sr",
                     signals_intercept=T,signals_slope=T,
                     lambda_true=1,true_model=F,optim_model=F,
                     pc.prior = "pc.prior.auto",
                     INLA_formula,
                     hyper = NULL) {
  require(INLA)
  require(nlme)

  C <- get_comm_pair_r(comm,V)

  V_true <- V*lambda_true
  diag(V_true) <- diag(V)
  C_true <- get_comm_pair_r(comm,V_true)

  if (signals_X == "sr") {
    x <- rowSums(comm) #species richness
  }

  if (signals_X == "phy_cor") {
    x <- c(t(chol(C_true)) %*% rnorm(n=nrow(comm),mean=ef_mean,sd=sd)) #signals in x
  }

  if (signals_X == "no_phy_cor") {
    x <- rnorm(nrow(comm),0,1) #no signals in X
  }

  if (signals_intercept == T & signals_slope==T) {
    library(MASS)
    RE <- mvrnorm(length(x), mu=c(0,0), Sigma=rbind(c(1.0, 0),
                                                   c(0, 1.0))) #simulate random effect
    comm_int <- c(t(chol(C_true))%*%RE[,1])
    comm_beta <- c(t(chol(C_true))%*%RE[,2])
    sr_E_cor <- cor(x,comm_beta*x+comm_int)
    slope_int_cor <- cor(comm_beta,comm_int)
    y <- b1*x+comm_beta*x+comm_int #signals in intercept and slope
  }

  if (signals_intercept == T & signals_slope == F) {
    comm_int <- c(t(chol(C_true)) %*% rnorm(n=nrow(comm),mean=ef_mean,sd=sd))
    y <- b1*x+comm_int #signals in intercept but not slope
    sr_E_cor <- NA
    slope_int_cor <- NA
  }

  if (signals_intercept==F & signals_slope == F) {
    y <- b1*x+rnorm(nrow(comm))#no signals in error
    sr_E_cor <- NA
    slope_int_cor <- NA
  }

  sim_data <- data.frame(y=y,x=x,comm=1:nrow(comm),comm2=1:nrow(comm))

  m <- lm(y~x,data=sim_data)
  p <- summary(m)$coefficients[1:2,4]
  effect_lm <- summary(m)$coefficients[1:2,1]

  ### INLA
  sdres <- sd(residuals(m))
  param <- c(3*sdres,0.01)

  m2 <- gls(y~x,data=sim_data,correlation=corSymm(C[lower.tri(C)], fixed = T))
  p_gls <- summary(m2)$tTable[1:2,4]
  effect_gls <- summary(m2)$tTable[1:2,1]

  if (true_model == T) {
  m2_true <- gls(y~x,data=sim_data,correlation=corSymm(C_true[lower.tri(C_true)], fixed = T))
  p_gls_true <- summary(m2_true)$tTable[1:2,4]
  effect_gls_true <- summary(m2_true)$tTable[1:2,1]
  } else {
    p_gls_true <- effect_gls_true <- c(NA,NA)
  }

  if (optim_model == T) {
  message("searching best lambda")
  require(phylosignal)
  #init_lambda <- round(lambdaTest(residuals(m),C)$Lambda,3)

  ###GLS-optim
  start_GLS <- Sys.time()
  ML.opt<-optim(runif(1,0.2,0.8),likelihood.lambda,y=y,X=x,V=V,comm=comm,method="L-BFGS-B",
                 lower=0.0,upper=1)
  lambda_GLS<-ML.opt$par
  logL<--ML.opt$value
  V_GLS <- V*lambda_GLS
  diag(V_GLS) <- diag(V)
  C.lambda.GLS<- get_comm_pair_r(comm,V_GLS)
  end_GLS <- Sys.time()
  time_GLS <- end_GLS-start_GLS

  ###INLA - optim
  start_INLA <- Sys.time()

  prior1 <- list(prec=list(prior="pc.prec"),param=param)

  ML.opt2<-optim(runif(1,0.2,0.8),
                 likelihood.lambda.INLA,
                 inla_formula= INLA_formula,
                 data=sim_data,
                 phyV=V,
                 comm=comm,
                 prior = list(prior1 = prior1),
                 control.compute = list(waic=T),
                 method="L-BFGS-B",
                 lower=0.0,upper=1.0)

  V2<-diag(diag(C))
  C_temp2<-C-V2
  n<-nrow(C_temp2)
  lambda_INLA<-ML.opt2$par
  logL2 <--ML.opt2$value
  V_INLA <- V*lambda_INLA
  diag(V_INLA) <- diag(V)
  C.lambda.INLA<- get_comm_pair_r(comm,V_INLA)
  end_INLA <- Sys.time()
  time_INLA <- end_INLA-start_INLA

  m2_optim <- gls(y~x,data=sim_data,correlation=corSymm(C.lambda.GLS[lower.tri(C.lambda.GLS)], fixed = T))
  p_gls_optim <- 2*pt(abs(summary(m2_optim)$tTable[1:2,3]),m2_optim$dims$N-m2_optim$dims$p-1,lower.tail=F) #estimated lambda, thus df = n-3.
  effect_gls_optim <- summary(m2_optim)$tTable[1:2,1]

  phylo_prec_mat_lambda <- solve(C.lambda.INLA)

  if (signals_slope == T) {
  m_INLA_lambda <- inla(y~x+f(comm,model="generic0",Cmatrix=phylo_prec_mat_lambda,prior = "pc.prec",hyper=param)+
                          f(comm2,x,model="generic0",Cmatrix=phylo_prec_mat_lambda,prior = "pc.prec",hyper=param),data=sim_data)
  } else {
    m_INLA_lambda <- inla(y~x+f(comm,model="generic0",Cmatrix=phylo_prec_mat_lambda,prior = "pc.prec",hyper=param),data=sim_data)

  }
  INLA_C_lambda<- m_INLA_lambda$summary.fixed
  INLA_C_lambda_sig <- sign(INLA_C_lambda[,3]) * sign(INLA_C_lambda[,5])*-1 #same sign (not overlap with zero) = negative number. easier for counting
  } else {
    lambda_INLA <- lambda_GLS <- time_GLS <- time_INLA <- NA
    INLA_C_lambda_sig <- p_gls_optim <- effect_gls_optim <-  c(NA,NA)
    INLA_C_lambda <- data.frame(c(NA))
  }

  if (true_model == T) {
  phylo_prec_mat_true <- solve(C_true)
  if (signals_slope == T) {
  m_INLA_true <- inla(y~x+f(comm,model="generic0",Cmatrix=phylo_prec_mat_true,prior = "pc.prec",hyper=param)+
                        f(comm2,x,model="generic0",Cmatrix=phylo_prec_mat_true,prior = "pc.prec",hyper=param),data=sim_data)
  } else {
    m_INLA_true <- inla(y~x+f(comm,model="generic0",Cmatrix=phylo_prec_mat_true,prior = "pc.prec",hyper=param),data=sim_data)

  }


  INLA_C_true <- m_INLA_true$summary.fixed
  INLA_C_true_sig <- sign(INLA_C_true[,3]) * sign(INLA_C_true[,5])*-1
  } else {
    INLA_C_true_sig <-  c(NA,NA)
    INLA_C_true <- data.frame(c(NA))

  }

  phylo_prec_mat <- solve(C)

  if (signals_slope == T) {
    m_INLA <- inla(y~x+f(comm,model="generic0",Cmatrix=phylo_prec_mat,prior = "pc.prec",hyper=param)+
                          f(comm2,x,model="generic0",Cmatrix=phylo_prec_mat,prior = "pc.prec",hyper=param),data=sim_data)
  } else {
  m_INLA <- inla(y~x+f(comm,model="generic0",Cmatrix=phylo_prec_mat,prior = "pc.prec",hyper=param),data=sim_data)
  }

  INLA_C <- m_INLA$summary.fixed
  INLA_C_sig <- sign(INLA_C[,3]) * sign(INLA_C[,5])*-1

  r <- cor(sim_data$x,sim_data$y)

  end <- Sys.time()
  message("best lambda = ",lambda_INLA)
  return(data.frame(p[[2]],p_gls[[2]],p_gls_true[[2]],p_gls_optim[[2]],INLA_C_sig[[2]],INLA_C_true_sig[[2]],INLA_C_lambda_sig[[2]],
                    effect_lm[[2]],effect_gls[[2]],effect_gls_true[[2]],effect_gls_optim[[2]],INLA_C[2,1],INLA_C_true[2,1],INLA_C_lambda[2,1],
                    p[[1]],p_gls[[1]],p_gls_true[[1]],INLA_C_sig[[1]],INLA_C_true_sig[[1]],INLA_C_lambda_sig[[1]],
                    effect_lm[[1]],effect_gls[[1]],effect_gls_true[[1]],INLA_C[1,1],INLA_C_true[1,1],INLA_C_lambda[1,1],
                    sr_E_cor,slope_int_cor,
                    lambda_GLS=lambda_GLS,lambda_INLA = lambda_INLA,lambda_true=lambda_true, x_y_r = r,
                    time_GLS=as.numeric(time_GLS),time_INLA = as.numeric(time_INLA)))
}
