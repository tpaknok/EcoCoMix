.get_variance_information <- function(model,
                                      faminfo,
                                      name_fun = "get_variances",
                                      verbose = TRUE) {
  # sanity check
  if (is.null(model)) {
    return(NULL)
  }


  check_if_installed("lme4", reason = "to compute variances for mixed models")

  if (inherits(model, "lme")) {
    check_if_installed("nlme", reason = "to compute variances for mixed models")
  }

  if (inherits(model, "rstanarm")) {
    check_if_installed("rstanarm", reason = "to compute variances for mixed models")
  }

  # stanreg
  # ---------------------------
  if (inherits(model, "stanreg")) {
    mixed_effects_info <- list(
      beta = lme4::fixef(model),
      X = rstanarm::get_x(model),
      vc = lme4::VarCorr(model),
      re = lme4::ranef(model)
    )

    # GLMMapdative
    # ---------------------------
  } else if (inherits(model, "MixMod")) {
    vc1 <- vc2 <- NULL
    re_names <- find_random(model)

    vc_cond <- !startsWith(colnames(model$D), "zi_")
    if (any(vc_cond)) {
      vc1 <- model$D[vc_cond, vc_cond, drop = FALSE]
      attr(vc1, "stddev") <- sqrt(diag(vc1))
      attr(vc1, "correlation") <- stats::cov2cor(model$D[vc_cond, vc_cond, drop = FALSE])
    }

    vc_zi <- startsWith(colnames(model$D), "zi_")
    if (any(vc_zi)) {
      colnames(model$D) <- gsub("^zi_(.*)", "\\1", colnames(model$D))
      rownames(model$D) <- colnames(model$D)
      vc2 <- model$D[vc_zi, vc_zi, drop = FALSE]
      attr(vc2, "stddev") <- sqrt(diag(vc2))
      attr(vc2, "correlation") <- stats::cov2cor(model$D[vc_zi, vc_zi, drop = FALSE])
    }

    model_deviance <- get_deviance(model, verbose = FALSE)
    residual_df <- get_df(model, type = "residual", verbose = FALSE)

    vc1 <- list(vc1)
    names(vc1) <- re_names[[1]]
    attr(vc1, "sc") <- sqrt(abs(model_deviance) / residual_df)

    if (!is.null(vc2)) {
      vc2 <- list(vc2)
      names(vc2) <- re_names[[2]]
      attr(vc2, "sc") <- sqrt(abs(model_deviance) / residual_df)
    }

    vcorr <- compact_list(list(vc1, vc2))
    names(vcorr) <- c("cond", "zi")[seq_along(vcorr)]

    mixed_effects_info <- list(
      beta = lme4::fixef(model),
      X = get_modelmatrix(model),
      vc = vcorr,
      re = list(lme4::ranef(model))
    )
    names(mixed_effects_info$re) <- model$id_name

    # joineRML
    # ---------------------------
  } else if (inherits(model, "mjoint")) {
    re_names <- find_random(model, flatten = TRUE)
    vcorr <- summary(model)$D
    attr(vcorr, "stddev") <- sqrt(diag(vcorr))
    attr(vcorr, "correlation") <- stats::cov2cor(vcorr)
    vcorr <- list(vcorr)
    names(vcorr) <- re_names[1]
    attr(vcorr, "sc") <- model$coef$sigma2[[1]]

    mixed_effects_info <- list(
      beta = lme4::fixef(model),
      X = matrix(1, nrow = n_obs(model), dimnames = list(NULL, "(Intercept)_1")),
      vc = vcorr,
      re = list(lme4::ranef(model))
    )
    names(mixed_effects_info$re) <- re_names[seq_along(mixed_effects_info$re)]

    # nlme / glmmPQL
    # ---------------------------
  } else if (inherits(model, "lme")) {
    re_names <- find_random(model, split_nested = TRUE, flatten = TRUE)
    comp_x <- get_modelmatrix(model)
    rownames(comp_x) <- seq_len(nrow(comp_x))
    if (.is_nested_lme(model)) {
      vals_vc <- .get_nested_lme_varcorr(model)
      vals_re <- lme4::ranef(model)
    } else {
      vals_vc <- list(nlme::getVarCov(model))
      vals_re <- list(lme4::ranef(model))
    }
    mixed_effects_info <- list(
      beta = lme4::fixef(model),
      X = comp_x,
      vc = vals_vc,
      re = vals_re
    )
    names(mixed_effects_info$re) <- re_names
    names(mixed_effects_info$vc) <- re_names

    # ordinal
    # ---------------------------
  } else if (inherits(model, "clmm")) {
    if (requireNamespace("ordinal", quietly = TRUE)) {
      mm <- get_modelmatrix(model)
      mixed_effects_info <- list(
        beta = c("(Intercept)" = 1, stats::coef(model)[intersect(names(stats::coef(model)), colnames(mm))]),
        X = mm,
        vc = ordinal::VarCorr(model),
        re = ordinal::ranef(model)
      )
    }

    # glmmadmb
    # ---------------------------
  } else if (inherits(model, "glmmadmb")) {
    mixed_effects_info <- list(
      beta = lme4::fixef(model),
      X = get_modelmatrix(model),
      vc = lme4::VarCorr(model),
      re = lme4::ranef(model)
    )

    # brms
    # ---------------------------
  } else if (inherits(model, "brmsfit")) {
    comp_x <- get_modelmatrix(model)
    rownames(comp_x) <- seq_len(nrow(comp_x))
    vc <- lapply(names(lme4::VarCorr(model)), function(i) {
      element <- lme4::VarCorr(model)[[i]]
      if (i != "residual__") {
        if (is.null(element$cov)) {
          out <- as.matrix(drop(element$sd[, 1])^2)
          colnames(out) <- rownames(out) <- gsub("Intercept", "(Intercept)", rownames(element$sd), fixed = TRUE)
        } else {
          out <- as.matrix(drop(element$cov[, 1, ]))
          colnames(out) <- rownames(out) <- gsub("Intercept", "(Intercept)", rownames(element$cov), fixed = TRUE)
        }
        attr(out, "sttdev") <- element$sd[, 1]
      } else {
        out <- NULL
      }
      out
    })
    vc <- compact_list(vc)
    names(vc) <- setdiff(names(lme4::VarCorr(model)), "residual__")
    attr(vc, "sc") <- lme4::VarCorr(model)$residual__$sd[1, 1]
    mixed_effects_info <- list(
      beta = lme4::fixef(model)[, 1],
      X = comp_x,
      vc = vc,
      re = lapply(lme4::ranef(model), function(re) {
        reval <- as.data.frame(drop(re[, 1, ]))
        colnames(reval) <- gsub("Intercept", "(Intercept)", dimnames(re)[[3]], fixed = TRUE)
        reval
      })
    )
    names(mixed_effects_info$beta) <- gsub("Intercept", "(Intercept)", names(mixed_effects_info$beta), fixed = TRUE)

    # cpglmm
    # ---------------------------
  } else if (inherits(model, "cpglmm")) {
    check_if_installed("cplm")

    mixed_effects_info <- list(
      beta = cplm::fixef(model),
      X = cplm::model.matrix(model),
      vc = cplm::VarCorr(model),
      re = cplm::ranef(model)
    )

    # lme4 / glmmTMB
    # ---------------------------
  } else {
    mixed_effects_info <- list(
      beta = lme4::fixef(model),
      X = lme4::getME(model, "X"),
      vc = lme4::VarCorr(model),
      re = lme4::ranef(model)
    )
  }


  # for models with zero-inflation, we only want the conditional part
  if (inherits(model, c("glmmTMB", "MixMod"))) {
    mixed_effects_info <- lapply(mixed_effects_info, .collapse_cond)
  }

  # currently, we don't support calculating all variance components
  # for the zero-inflated part of the model only. This is not fully implemented
  # mixed_effects_info <- lapply(mixed_effects_info, .collapse_zi)

  # for glmmTMB, tell user that dispersion model is ignored
  if (!is.null(find_formula(model)[["dispersion"]]) && verbose) {
    format_warning(sprintf("%s ignores effects of dispersion model.", name_fun))
  }

  # fix rank deficiency
  rankdef <- is.na(mixed_effects_info$beta)
  if (any(rankdef)) {
    rankdef_names <- names(mixed_effects_info$beta)[rankdef]
    mixed_effects_info$beta <- mixed_effects_info$beta[setdiff(names(mixed_effects_info$beta), rankdef_names)]
  }

  mixed_effects_info
}

.collapse_cond <- function(x) {
  if (is.list(x) && "cond" %in% names(x)) {
    x[["cond"]]
  } else {
    x
  }
}

.compute_variance_random <- function(model, terms, mixed_effects_info) {
  if (is.null(terms)) {
    return(NULL)
  }
  .sigma_sum <- function(Sigma) {
    rn <- rownames(Sigma)

    # fix for models w/o intercept
    if (!any(startsWith(colnames(mixed_effects_info$X), "(Intercept)"))) {
      mixed_effects_info$X <- cbind("(Intercept)" = 1, mixed_effects_info$X)
    }

    if (!is.null(rn)) {
      valid <- rownames(Sigma) %in% colnames(mixed_effects_info$X)
      if (!all(valid)) {
        rn <- rn[valid]
        Sigma <- Sigma[valid, valid, drop = FALSE]
      }
    }

    Z <- mixed_effects_info$X[, rn, drop = FALSE]
    Z.m <- Z %*% Sigma
    sum(diag(crossprod(Z.m, Z))) / n_obs(model)
  }

  # if (inherits(x, "MixMod")) {
  #   .sigma_sum(mixed_effects_info$vc)
  # } else {
  #   sum(sapply(mixed_effects_info$vc[terms], .sigma_sum))
  # }
  sum(sapply(mixed_effects_info$vc[terms], .sigma_sum))
}

.get_variance_information(m1)$vc


nr <- vapply(.get_variance_information(m1)$re, nrow, numeric(1))
not_obs_terms <- names(nr[nr != n_obs(m1)])
obs_terms <- names(nr[nr == n_obs(m1)])

var.random <- .compute_variance_random(m1, not_obs_terms, .get_variance_information(m1))

###
Sigma <- .get_variance_information(m1)$vc[not_obs_terms]
mixed_effects_info <- .get_variance_information(m1)
.get_variance_information(m1)$vc[not_obs_terms] %*% .get_variance_information(m1)$X

###
mu <- exp(c(fixef(null_model(m1)))$cond+0.5*0.5337962)
sig <- sigma(m1)

disp <- mu * (1 + sig)

log(1+1/mu+1/sig)
