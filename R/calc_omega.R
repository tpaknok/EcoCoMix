calc_omega <- function(INLA_m) {
  re <- 1/INLA_m$summary.hyperpar$mean
  
  INLA_formula <- INLA_m$.args$formula
  INLA_formula <-Reduce(paste,deparse(INLA_formula))
  
  formula_var <- unlist(strsplit(Reduce(paste,deparse(INLA_lambda_1[[1]]$best_model$.args$formula)),"\\+"))
  pos <- grep("prec.mat.INLA",formula_var)
  omega <- sum(re[pos])/(sum(re))
  
  return(omega)
}