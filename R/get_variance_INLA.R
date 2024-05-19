get_variance_INLA <- function(INLA_m,
                              INLA_m_int_only,
                              phylo_Cmatrix_names = "prec.mat.INLA"){

  require(vctrs)
  var_fixed <- NULL
  var_random <- NULL
  var_residual <- NULL
  var_phylo <- NULL

  ###fixed effect calculation
  message("get_fixed")
  beta <- INLA_m$summary.fixed$`0.5quant`
  X <- INLA_m$model.matrix

  var_fixed <- var(as.vector(beta%*%t(as.matrix(X))))

  ###random effect calculation
  message("get_random")
  var_random <- .get_variance_random_INLA(INLA_m)

  ###var_residual
  var_residual <- .get_variance_residual(INLA_m,INLA_m_int_only)

  ###var_phylo
  var_phylo <- .get_variance_random_phylo_INLA(INLA_m,phylo_Cmatrix_names = "prec.mat.INLA")

  ###output
  list(var_fixed=var_fixed,
       var_random=var_random,
       var_phylo=var_phylo,
       var_residual=var_residual
       )
}

####
.get_variance_random_INLA <- function (INLA_m) {
  require(stringr)
  all_vars <- labels(terms(INLA_m$.args$formula))
  re <- all_vars[grep("f\\(",all_vars)] #random effect

  generic_0_re <- re[grepl("Cmatrix",re)]

  if(length(generic_0_re) > 0) {
    generic_0_re <- gsub(" ","",generic_0_re)
    split_re <- strsplit(generic_0_re,",")
    generic_0_mat_names <- sapply(split_re, function(x) unlist(x)[grep("Cmatrix=",x)])
    generic_0_mat_names <- unname(sapply(generic_0_mat_names,function(x) gsub("Cmatrix=|\\)","",x)))
    for (i in 1:length(generic_0_mat_names)) {
      assign(generic_0_mat_names[[i]],diag(nrow(INLA_m$.args$data)))
    }
  }

  if (length(re) > 0) {
  ###constructure fake C matrix such that the f() function runs....
    re_list  <- lapply(re, function(x) eval(parse(text=x))[c("label","weights")])
    re_list <- purrr::transpose(re_list)
    re_list$var <- lapply(INLA_m$summary.hyperpar[match(paste0("Precision for ",re_list$label),rownames(INLA_m$summary.hyperpar)),"0.5quant"],
                          function(x) 1/x)
    re_list$weights <- lapply(re_list$weights, function(x) ifelse(is.null(x),"(Intercept)",x))

    Z <- INLA_m$.args$data[,unlist(re_list$label),drop=F]

    re_list$label <- lapply(colnames(Z)[vec_duplicate_id(t(Z))],function(x) x)
    unique.terms <- unique(unlist(re_list$label))

    var_random <- sum(sapply(unique.terms,function(x) {
      pos <- which(re_list$label == x)
      re_subset<- lapply(re_list,function(x) x[pos])

      Z <- as.matrix(INLA_m$model.matrix[,unlist(unique(re_subset$weights))])
      Sigma <- matrix(0,nrow=ncol(Z),ncol=ncol(Z))
      diag(Sigma) <- unlist(re_subset$var)
      Z.m <- Z %*% Sigma

      sum(diag(crossprod(Z.m, Z))) / nrow(Z)
    }
    )
    )
  } else {
    message("No random effect detected...")
    var_random <- 0
  }

  return(var_random)
}

.get_variance_random_phylo_INLA <- function (INLA_m,phylo_Cmatrix_names = "prec.mat.INLA") {
  all_vars <- labels(terms(INLA_m$.args$formula))
  re <- all_vars[grep("f\\(",all_vars)] #random effect
  phylo_re <- all_vars[grep(phylo_Cmatrix_names,all_vars)]

  generic_0_re <- re[grepl("Cmatrix",phylo_re)]

  if(length(generic_0_re) > 0) {
    generic_0_re <- gsub(" ","",generic_0_re)
    split_re <- strsplit(generic_0_re,",")
    generic_0_mat_names <- sapply(split_re, function(x) unlist(x)[grep("Cmatrix=",x)])
    generic_0_mat_names <- unname(sapply(generic_0_mat_names,function(x) gsub("Cmatrix=|\\)","",x)))
    for (i in 1:length(generic_0_mat_names)) {
      assign(generic_0_mat_names[[i]],diag(nrow(INLA_m$.args$data)))
    }
  }

  if (length(phylo_re) > 0) {
    re_list  <- lapply(phylo_re, function(x) eval(parse(text=x))[c("label","weights")])
    re_list <- purrr::transpose(re_list)
    re_list$var <- lapply(INLA_m$summary.hyperpar[match(paste0("Precision for ",re_list$label),rownames(INLA_m$summary.hyperpar)),"0.5quant"],
                          function(x) 1/x)
    re_list$weights <- lapply(re_list$weights, function(x) ifelse(is.null(x),"(Intercept)",x))

    Z <- INLA_m$.args$data[,unlist(re_list$label),drop=F]

    re_list$label <- lapply(colnames(Z)[vec_duplicate_id(t(Z))],function(x) x)
    unique.terms <- unique(unlist(re_list$label))

    var_phylo_random <- sum(sapply(unique.terms,function(x) {
      pos <- which(re_list$label == x)
      re_subset<- lapply(re_list,function(x) x[pos])

      Z <- as.matrix(INLA_m$model.matrix[,unlist(unique(re_subset$weights))])
      Sigma <- matrix(0,nrow=ncol(Z),ncol=ncol(Z))
      diag(Sigma) <- unlist(re_subset$var)
      Z.m <- Z %*% Sigma

      sum(diag(crossprod(Z.m, Z))) / nrow(Z)
    }
    )
    )
  } else {
    message("No phylogenetic random effect detected...")
    var_phylo_random <- 0
  }

  return(var_phylo_random)
}

.get_variance_residual <- function (INLA_m,
                                    INLA_m_int_only) {
  supported_family = c("gaussian","tweedie","poisson","nbinomial","binomial","zeroinflatednbinomial1")

  if (INLA_m$.args$family %in% supported_family) {
  mu <- exp(INLA_m_int_only$summary.fixed$`0.5quant`)
  } else {
    stop("Currently the package doesn't support ",INLA_m$.args$family," distribution.")
  }

  if (mu < 6) {
    warning("mu is close to zero. Estimates of distributional residuals (and hence other R2) may not be reliable?")
  }
  cvsquared <- tryCatch(
    {
      vv<-switch(INLA_m$.args$family,
             gaussian = (1/INLA_m$summary.hyperpar["Precision for the Gaussian observations","0.5quant"]),
             nbinomial = switch(INLA_m$misc$linkfunctions$names,
                                log=.get_variance_residual_nbinomial(INLA_m,mu),
                                sqrt=0.25,
                                "Only supporting log and sqrt link for negative binomial distribution."),
             poisson = switch(INLA_m$misc$linkfunctions$names,
                              log=mu,
                              sqrt=0.25,
                              "Only supporting log and sqrt link for poisson distribution."),
             tweedie = switch(INLA_m$misc$linkfunctions$names,
                              log= .get_variance_residual_tweedie(INLA_m,mu),
                              "Only supporting log link for tweedie distribution."),
             binomial = switch(INLA_m$misc$linkfunctions$names,
                               logit=pi^2/3,
                               cloglog=,
                               clogloglink=pi^2/6,
                               probit=1,
                               "Only supporting logit, clogloglink, and probit link for tweedie distribution."),
             zeroinflatednbinomial1 = switch(INLA_m$misc$linkfunctions$names,
                                             log=.get_variance_residual_zinfnb(INLA_m),
                                             "Only supporting log link for zero-inflated negative binomial distribution 1."
             )
      )

      if (INLA_m$.args$family == "gaussian") {
        vv <- vv
      } else {
        vv <- vv / mu^2
      }
    }
  )

  if (INLA_m$.args$family == "gaussian") {
    cvsquared
  } else {
    log1p(cvsquared)
  }
}

  .get_variance_residual_nbinomial <- function (INLA_m,mu) {
    theta <- INLA_m$summary.hyperpar["size for the nbinomial observations (1/overdispersion)","0.5quant"]
    log1p(mu * (1 + mu/theta)/(mu^2))
  }

  .get_variance_residual_tweedie <- function (INLA_m,mu) {
    phi <- INLA_m$summary.hyperpar["Dispersion parameter for Tweedie","0.5quant"]
    power <- INLA_m$summary.hyperpar["p parameter for Tweedie","0.5quant"]
    phi * mu^power
  }

  .get_variance_residual_zinfnb <- function(INLA_m,mu) {
    # zi probability
    p <- INLA_m$summary.hyperpar["zero-probability parameter for zero-inflated nbinomial_1","0.5quant"]

    # median of conditional distribution
    conditional <- INLA_m$summary.fitted.values[,"0.5quant"]/(1-p)

    # sigma
    k <- INLA_m$summary.hyperpar["size for nbinomial_1 zero-inflated observations","0.5quant"]

    # variance
    var <- conditional * (1 + conditional/k)

    pvar <- (1 - p) * var + conditional^2 * (p^2 + p)
    mean(pvar)
  }



