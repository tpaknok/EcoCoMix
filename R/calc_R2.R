calc_R2 <- function(INLA_m) {
  re <- 1/INLA_m$summary.hyperpar$mean
  var_re <- sum(re)

  fe_names <- INLA_m$names.fixed[-grep("(Intercept)",INLA_m$names.fixed)]
  data_fe <- as.matrix(INLA_m$model.matrix[,fe_names])
  var_predictor <- apply(data_fe,2,var)
  var_fe <- sum(var_predictor * INLA_m$summary.fixed[fe_names,]$mean^2)

  R2m <- 1-(var_re)/(var_fe+var_re)
  R2c <- 1-(re[[1]])/(var_fe+var_re)

  R2_vec <- c(R2m=R2m,R2c=R2c)
  return(R2_vec)
}


