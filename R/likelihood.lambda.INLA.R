likelihood.lambda.INLA<-function(inla_formula,data,phyV=NULL,comm=NULL,lambda,prior = NULL,...){
  require(INLA)
  message(lambda, "-INLA")
  V_lambda <- phyV*lambda
  diag(V_lambda) <- diag(phyV)
  C.lambda<- get_comm_pair_r(comm,V_lambda)

  P.lambda <- solve(C.lambda)

  unlist(unname(prior[1]),recursive=F)

  # create prior list in local environment
  for (i in 1:length(prior)) {
    if (lengths(prior[i]) > 1) {
    assign(names(prior[i]),unlist(unname(prior[i]),recursive=F))
    }

    if (length(prior[i]) == 1) {
      assign(names(prior[i]),unname(prior[i])[[1]])
    }
  }
  #res_f <- gsub("\\bf\\([^)]*\\)", "",Reduce(paste, deparse(formula)))
  #res_f <- gsub("\\+\\s*", "", res_f)

  #if (prior == "pc.prior.auto") {
  #m <- lm(res_f,data=data)
  #sdres <- sd(residuals(m))
  #hyper_param <- c(3*sdres,0.01)
  #} else {
    #hyper_param <- hyper
  #}

  #argus <- c(list(formula = as.formula(formula),
                  #data = data))

  #m_INLA_lambda <- do.call(inla,argus)

  #formula <- as.formula(gsub("hyper_param","hyper_param",Reduce(paste, deparse(formula)))) #so replace it with a new hyper_param works??

  inla_formula <- as.formula(Reduce(paste, deparse(inla_formula)))

  m_INLA_lambda <- inla(inla_formula,
                        data=data,...)

  m_INLA_lambda <- inla.rerun(m_INLA_lambda)
  #summary(m_INLA_lambda)

  wAIC <- m_INLA_lambda$waic$waic
  return(wAIC)
}

#Example
#optim(runif(1),likelihood.lambda.INLA,inla_formula=f,data=data1,phyV=V_sp,
#comm=comm, prior=list(prior1=prior1),
#control.compute = list(waic=T),
#lower=0,upper=1)
