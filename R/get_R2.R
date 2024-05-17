get_R2 <- function(variance_componenet) {
  var_comp <- variance_componenet

  total_var <- var_comp["var_fixed"]+var_comp["var_random"]+var_comp["var_residual"]
  r2m <- var_comp["var_fixed"]/total_var
  r2c <- (var_comp["var_fixed"]+var_comp["var_random"])/total_var
  r2Phylo <- var_comp["var_phylo"]/total_var
  r2Phylo_latent <- var_comp["var_phylo"]/(var_comp["var_phylo"]+var_comp["var_residual"])

  c(R2marginal=unname(r2m),R2conditional=unname(r2c),R2Phylo=unname(r2Phylo),R2Phylo_latent=unname(r2Phylo_latent))
}
