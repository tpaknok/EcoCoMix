likelihood.lambda.INLA<-function(inla_formula,
                                 data,
                                 family="gaussian",
                                 VCV_sp=NULL,
                                 comm=NULL,
                                 lambda,
                                 priors = NULL,
                                 ...){
  require(INLA)
  V_lambda <- VCV_sp*lambda
  diag(V_lambda) <- diag(VCV_sp)
  C.lambda<- get_comm_pair_r(comm,V_lambda)

  Phylo <- solve(C.lambda)

  # create prior list in local environment
  if (!is.null(priors)) {
    for (i in 1:length(priors)) {
      if (lengths(priors[i]) > 1) {
      assign(names(priors[i]),unlist(unname(priors[i]),recursive=F))
      }

    if (length(priors[i]) == 1) {
      assign(names(priors[i]),unname(priors[i])[[1]])
      }
    }
  }

  inla_formula <- as.formula(Reduce(paste, deparse(inla_formula)))

  m_INLA_lambda <- inla(inla_formula,
                        family=family,
                        data=data,...)

  m_INLA_lambda <- inla.rerun(m_INLA_lambda)

  wAIC <- m_INLA_lambda$waic$waic
  return(wAIC)
}
