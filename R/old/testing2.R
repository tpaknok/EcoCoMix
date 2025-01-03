dat <- data.frame(y =c(rZINBI(100, mu = 10, sigma = .6, nu=0.1),rZINBI(100, mu = 5, sigma = .3, nu=.5)),sites =c(rep("a", 100), rep("b", 100)),year = rep(1:4, each = 10, times = 5),trans = rep(1:40, each = 5, times = 1), area=rNO(200,20))
head(dat)

(m2 <- glmmTMB(count ~ scale(DOY) + (scale(DOY)||site),
               family=nbinom2, data=Salamanders))

performance::r2(m2)

insight::get_variance(m2)


library(spaMM)
library(lmerTest)
?fitme
data("blackcap")
(fullfit <- fitme(migStatus ~ means+ Matern(1|longitude+latitude),data=blackcap) )
(fullfit_int <- fitme(migStatus ~ 1+ Matern(1|longitude+latitude),data=blackcap) )

get_variance_spaMM(fullfit,fullfit_int)

data("scotlip")
m <- fitme(cases~I(prop.ag/10)+offset(log(expec)),family=negbin(), data=scotlip)
c(m$family)

str(m)

m$family$resid.model


library(lme4)

data(Orthodont,package="nlme")
Orthodont$nsex <- as.numeric(Orthodont$Sex=="Male")
Orthodont$nsexage <- with(Orthodont, nsex*age)
m <- lmer(distance ~ age + (age||Subject), data=Orthodont)
m_TMB <- glmmTMB(distance ~ age + nsex+ (age+nsex||Subject), data=Orthodont,REML=T)

m_spaMM <- fitme(distance ~ age + nsex +(age||Subject), data=Orthodont,method="REML")

performance::r2(m)
performance::r2(m_TMB)
get_variance_spaMM(m_spaMM)

Salamanders$mined <- ifelse(Salamanders$mined == "yes", 1, 0)

(m1 <- glmmTMB(count ~ mined + (mined||site)+(mined||spp),
               family=nbinom2, data=Salamanders))
summary(m1)
get_variance(m1)
performance::r2(m1)

m_spaMM2 <- fitme(count ~ mined + (mined||site) + (mined||spp), data=Salamanders,family=negbin,method="ML")
summary(m_spaMM2)

m_spaMM3 <- fitme(count ~ 1+ (mined||site) + (mined||spp), data=Salamanders,family=negbin,method="ML")
summary(m_spaMM3)

get_variance_spaMM(m_spaMM2,m_spaMM3)

get_variance(m1)

###
mixed_effects_info <- list(
  beta = lme4::fixef(model),
  X = lme4::getME(model, "X"),
  vc = lme4::VarCorr(model),
  re = lme4::ranef(model)
)

.collapse_cond <- function(x) {
  if (is.list(x) && "cond" %in% names(x)) {
    x[["cond"]]
  } else {
    x
  }
}

if (inherits(model, c("glmmTMB", "MixMod"))) {
  mixed_effects_info <- lapply(mixed_effects_info, .collapse_cond)
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
