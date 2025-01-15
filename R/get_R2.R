get_R2 <- function(spaMM_m,spaMM_m_int_only = NULL) {

  get_variance_spaMM <- function(spaMM_m, spaMM_m_int_only){

    ###checked with glmmTMB output based on performance::r2
    require(vctrs)
    var_fixed <- NULL
    var_random <- NULL
    var_residual <- NULL
    var_phylo <- NULL

    ###fixed effect calculation
    beta <- fixef(spaMM_m)
    X <- get_matrix(spaMM_m)

    var_fixed <- var(as.vector(beta%*%t(as.matrix(X))))
    ###random effect calculation

    var_random <- .get_variance_random_spaMM(spaMM_m)

    ###var_residual
    var_residual <- .get_variance_residual(spaMM_m,spaMM_m_int_only)
    # ###var_phylo
    # var_phylo <- .get_variance_random_phylo_INLA(spaMM_m,phylo_Cmatrix_names = "prec.mat.INLA")
    #
    ###output
    list(var_fixed=var_fixed,
         var_random=var_random,
         var_phylo=var_phylo,
         var_residual=var_residual
    )
  }

  ####
  .get_variance_random_spaMM <- function (spaMM_m) {

    vc <- VarCorr(spaMM_m)

    vc$Group <- sub("\\..*", "", vc$Group)
    terms <- unique(vc$Group)

    vc_list <- list()
    for (i in 1:length(terms)) {
      vc_terms <- subset(vc, Group == terms[i])
      vc_list[[i]] <- vc_terms$Variance
      names(vc_list[[i]]) <- vc_terms$Term
    }

    names(vc_list) <- terms

    vc_list <- vc_list[names(vc_list) != "Residual"]
    re_var <- NULL
    for (j in 1:(length(vc_list))) {
      Sigma <- vc_list[[j]]
      rn <- names(Sigma)

      X <- spaMM_m$data
      X$`(Intercept)` <- 1

      if (!is.null(rn)) {
        valid <- names(Sigma) %in% colnames(X)
        if (!all(valid)) {
          rn <- rn[valid]
          Sigma <- Sigma[valid, valid, drop = FALSE]
        }
      }

      if (length(Sigma) > 1) {
      Sigma<-diag(Sigma)
      }

      Z <- as.matrix(X[, rn, drop = FALSE])
      Z.m <- as.matrix(Z) %*% Sigma

      re_var <- c(re_var,sum(diag(crossprod(Z.m, Z))) / nrow(spaMM_m$data))
    }

    re_var <- sum(re_var)
  }

  ####
  .get_variance_residual <- function (spaMM_m,
                                      spaMM_m_int_only) {

    if (family(spaMM_m)$family == "gaussian") {
      mu <- 1
    }

    if (family(spaMM_m)$family == "negbin2") { #experimental

      mu <- exp(as.vector(fixef(spaMM_m_int_only)) + 0.5 * .get_variance_random_spaMM(spaMM_m_int_only))
      sig <- get_inits_from_fit(spaMM_m)$init$NB_shape
    }

    cvsquared <- tryCatch(
      {
        vv<-switch(family(spaMM_m)$family,
                   gaussian = VarCorr(spaMM_m)[VarCorr(spaMM_m)$Group == "Residual","Variance"],
                   negbin2 = log(1+1/mu+1/sig)
        )
      }
    )
  }

  if (spaMM_m$family$family != "gaussian" & spaMM_m$family$family != "negbin2") {
    stop("Currently the function only supports gaussian distribution.")
  }

  if (length(ranef(spaMM_m)) == 0) {
    obs_R2 <- 1-sum((spaMM_m$y-predict(spaMM_m))^2)/sum((spaMM_m$y-mean(spaMM_m$y))^2)
    p <- ncol(get_matrix(spaMM_m)[,colnames(get_matrix(spaMM_m)) != "(Intercept)",drop=F])
    N <- length(spaMM_m$y)
    adj_R2 <- 1-(1-obs_R2)*(N-1)/(N-p-1)
    output <- c(`R-squared` = obs_R2, `Adjusted R-squared` = adj_R2)

  } else {
  m_var <- get_variance_spaMM(spaMM_m,spaMM_m_int_only)

  R2m <- m_var$var_fixed/(m_var$var_fixed+m_var$var_random+m_var$var_residual)
  R2c <- (m_var$var_fixed+m_var$var_random)/(m_var$var_fixed+m_var$var_random+m_var$var_residual)
  output <- c(R2m=R2m,R2c=R2c)
  }
    return(output)
}
