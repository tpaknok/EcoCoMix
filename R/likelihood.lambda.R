likelihood.lambda<-function(lambda,
                            y,
                            X,
                            V,
                            comm){
  message(lambda)
  V_lambda <- V*lambda
  diag(V_lambda) <- diag(V)

  C.lambda<- get_comm_pair_r(comm,V_lambda)
  n <- ncol(C.lambda)

  sim_data <- data.frame(y=y,X=X)
  gls_m <- gls(y~X,data=sim_data,correlation=corSymm(C.lambda[lower.tri(C.lambda)], fixed = T))
  AIC <- AIC(gls_m)
  logL <- logLik(gls_m)
  #beta<-solve(t(X)%*%solve(C.lambda)%*%X)%*%(t(X)%*%solve(C.lambda)%*%y)
  #sig2e<-as.double((1/n)*(t(y-X%*%beta)%*%solve(C.lambda)%*%(y-X%*%beta)))
  #logL<--(1/2)*t(y-X%*%beta)%*%solve(sig2e*C.lambda)%*%(y-X%*%beta)-(1/2)*
    #determinant(sig2e*C.lambda,logarithm=TRUE)$modulus-(n/2)*log(2*pi)

  #AIC <- 2*(ncol(as.matrix(X))+2) - 2*logL
  return(AIC)
}
