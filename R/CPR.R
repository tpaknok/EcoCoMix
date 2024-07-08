CPR <- function(formula,
                priors = NULL,
                df,
                VCV_sp,
                comm,
                family,
                optim.lambda=T,
                original.VCV=T,
                optim.control=NULL,
                control.compute=list(waic=T,config=T),
                prediction.df = "auto",
                predictedfittedresponse = "best_m",
                inla.rerun = 1,
                wAIC_threshold = -2,
                ...
                ) {

  require(INLA)
  require(stringr)
  require(tidyverse)

  formula_text <- Reduce(paste,deparse(formula))
  response_name <- gsub(" ~.*","",formula_text)
  all_vars <- labels(terms(formula))

  community_names <- word(all_vars[grep("Phylo",all_vars)],sep=",")
  community_names <- gsub("f\\(","",community_names)

  fe <- all_vars[-grep("f\\(",all_vars)] #fixed effect
  re <- all_vars[grep("f\\(",all_vars)] #random effect

  if (length(re) == 0) {
    stop("need to specify at least one phylogenetic random effect (i.e., a generic 0 model with the Cmatrix argument sets as Phylo....such as f(comm,model='generic0',Cmatrix=Phylo)")
  }

  if (!grepl("Cmatrix=Phylo",gsub(" ","",re))){
    stop("need to specify at least one phylogenetic random effect (i.e., a generic 0 model with the Cmatrix argument sets as Phylo....such as f(comm,model='generic0',Cmatrix=Phylo)")
  }

  formula_elements <- unlist(strsplit(formula_text,split="\\+"))
  formula_elements <- formula_elements[!grepl("Phylo",formula_elements)]
  no_phylo_formula <- as.formula(paste(formula_elements,collapse="+"))
  message(response_name)

  df$predictor <- "observed"
  df$response <- response_name

  missing_community_names <- community_names[!community_names %in% colnames(df)]
  if (length(missing_community_names > 0)) { #create a community ID column if missing
  comm_id <- replicate(length(missing_community_names),1:nrow(df))
  colnames(comm_id) <- community_names
  df <- cbind(df,comm_id)
  }

  if (length(grep("Phylo", names(priors))) > 0) {
    stop('The names of your prior should not contain the word "Phylo". Please change the names.')
  }

  ###Create variable
  wAIC_GLM <- wAIC_optim <- wAIC_original <- NA
  lambda_INLA <- NA
  m_optim_result <- m_original_result <- NA
  prec.mat.INLA <- NA
  prediction <- NA
  prediction_phylo <- NA
  best_m <- NULL
  best_model <- "Starting"

  build_prediction <- function(x) {
    if (is.numeric(x) | is.integer(x)) {
      output <- mean(x)
    }

    if (is.factor(x) | is.character(x) ) {
      output <- names(sort(table(x),decreasing=T))[[1]]
    }
    return(output)
  }

  if (prediction.df == "auto") {
    predict_df <- lapply(fe,function(x) {
        marginal_predict <- as.data.frame(lapply(df, function(x) build_prediction(x)))
        unique_value <- unique(df[,x])
        marginal_predict <- marginal_predict[rep(seq_len(nrow(marginal_predict)), each = length(unique_value)), ]
        marginal_predict[,x] <- unique_value
        marginal_predict[,!colnames(marginal_predict) %in% fe] <- NA
        marginal_predict$predictor <- as.character(x)
        marginal_predict$response <- response_name
        return(marginal_predict)
      })

    predict_df <- do.call(rbind,predict_df)
    df <- rbind(df,predict_df)
  }

  if (is.data.frame(prediction.df)) {
    prediction.df[,response_name] <- NA
    df <- cbind(df,prediction.df)
  }

  pos <- which(is.na(df[,response_name])) #the position of the supplied df

  m1_INLA_GLM <- inla(no_phylo_formula,
                      data=df,
                      family=family,
                      control.predictor=list(compute=TRUE),
                      control.compute = control.compute,
                      ...)

  wAIC_GLM <- m1_INLA_GLM$waic$waic
  # m1_INLA_GLM_r2<- cor(m1_INLA_GLM$summary.fitted.values[-pos,"mean"],df[-pos,response_name])^2
  f <- get(paste0("inla.link.",m1_INLA_GLM$misc$linkfunctions$names)) #link function!!

  phylo_formula <- Reduce(paste,deparse(formula))
  phylo_formula <- gsub("Phylo","prec.mat.INLA",phylo_formula)
  phylo_formula <- as.formula(phylo_formula)
  int_only_formula <- paste(response_name,"~ 1+",re,collapse="+")
  int_only_formula <- gsub("Phylo","prec.mat.INLA",int_only_formula)
  int_only_formula <- as.formula(int_only_formula)

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

  if (original.VCV == T) {
    C.lambda.INLA<- get_comm_pair_r(comm,VCV_sp)
    prec.mat.INLA <- solve(C.lambda.INLA)

    m1_INLA_original <- inla(phylo_formula,
                              data = df,
                              family = family,
                              control.predictor=list(compute=TRUE),
                              control.compute = control.compute,
                              ...)

    m1_INLA_original_int_only <- inla(int_only_formula,
                                       data = df,
                                       family = family,
                                       control.predictor=list(compute=TRUE),
                                       control.compute = control.compute,
                                       ...)

    wAIC_original <- m1_INLA_original$waic$waic
    m1_INLA_original_r2<- cor(m1_INLA_original$summary.fitted.values[-pos,"mean"],df[-pos,response_name])^2
    m_original_result <- m1_INLA_original$summary.fixed
  }

  # if (priors == "pc.prec") {
  #   message("setup penalizing complexity priors")
  #   # m <- glm(no_phylo_formula,
  #   #          family = family,
  #   #          data = df)
  #   #
  #   # summary(m)
  #   # m_sig <- summary(m)$coefficients[2,4]
  #   # sdres <- sd(residuals(m))
  #   sdres <- sd(m1_INLA_GLM$summary.fitted.values[1:nrow(df),"mean"]-df[,response_name])
  #   param <- c(3*sdres,0.01)
  #   priors <- list(priors=list(prec=list(prior="pc.prec"),param=param))
  # }

  if (optim.lambda==T) {
    require(NMOF)

    grid_result <- gridSearch(fun=likelihood.lambda.INLA,
                              levels=list(lambda=c(0.05,0.25,0.5,0.75,0.95)),
                              lower=0,
                              upper=1,
                              inla_formula = formula,
                              family = family,
                              data = df,
                              VCV_sp = VCV_sp,
                              comm = comm,
                              prior = priors,
                              control.compute = control.compute,
                              printDetail = F,
                              ...) #Determine best starting value among 0.2,0.4,0.6,0.8

    ML.opt<-optim(grid_result$minlevels,
                  likelihood.lambda.INLA,
                  inla_formula = formula,
                  family = family,
                  data = df,
                  VCV_sp = VCV_sp,
                  comm = comm,
                  prior = priors,
                  control.compute = control.compute,
                  method = "L-BFGS-B",
                  lower = 0.0,
                  upper = 1.0,
                  control = optim.control,
                  ...) #advanced search using L-BFGS-B

    lambda_INLA<-ML.opt$par
    wAIC <- ML.opt$value
    VCV_sp_INLA <- VCV_sp*lambda_INLA
    diag(VCV_sp_INLA) <- diag(VCV_sp)

    C.lambda.INLA<- get_comm_pair_r(comm,VCV_sp_INLA)
    prec.mat.INLA <- solve(C.lambda.INLA)

    m1_INLA_optim <- inla(phylo_formula,
                          data = df,
                          family = family,
                          control.predictor=list(compute=TRUE),
                          control.compute = control.compute,
                          ...) #optimized_model

    m1_INLA_int_only <- inla(int_only_formula,
                          data = df,
                          family = family,
                          control.predictor=list(compute=TRUE),
                          control.compute = control.compute,
                          ...)

    wAIC_optim <- m1_INLA_optim$waic$waic
    # m1_INLA_optim_r2 <- cor(m1_INLA_optim$summary.fitted.values[-pos,"mean"],df[-pos,response_name])^2

    m_optim_result <- m1_INLA_optim$summary.fixed

    prediction_phylo <- cbind(df,
                              f(m1_INLA_optim$summary.fitted.values, inverse=T),
                              model="With phylogeny",
                              better_fit = ifelse(best_model == "With phylogeny","Yes","No"))
    prediction_phylo <- prediction_phylo[pos,]

    fe_sig <- ifelse(sign(m1_INLA_optim$summary.fixed[,3])* sign(m1_INLA_optim$summary.fixed[,5]) == 1, "Sig.","Insig.")
    fe_sig <- data.frame(predictor = rownames(m1_INLA_optim$summary.fixed),Sig=fe_sig)
    prediction_phylo <- prediction_phylo %>% inner_join(fe_sig,"predictor")
  }

  wAIC_all <- c(wAIC_optim=wAIC_optim, wAIC_no_phylo=wAIC_GLM, wAIC_original_VCV=wAIC_original)

  best_model_name <- "Phylogeny without optimization"
  return_lambda <- NA

  if (optim.lambda == T) {
    return_lambda <- lambda_INLA
    if (wAIC_optim - wAIC_GLM >= wAIC_threshold) {
      best_m <- m1_INLA_GLM
      best_model_name <- "Without phylogeny"
    } else { #will ignore NA
      best_model_name <- "Optimized phylogeny"
      }
    }

  min_wAIC_m <- names(which.min(wAIC_all))
  best_m <- switch(min_wAIC_m,
                   wAIC_optim = m1_INLA_optim,
                   wAIC_no_phylo = m1_INLA_GLM,
                   wAIC_original_VCV = m1_INLA_original)

  rerun <- 1
  if(inla.rerun != 0) {
    while (rerun <= inla.rerun) {
    message("rerunning INLA...",rerun,"/",inla.rerun) # https://groups.google.com/g/r-inla-discussion-group/c/qkoV9ZtA1Wo
    m1_INLA_GLM <- inla.rerun(m1_INLA_GLM)
    if (optim.lambda == T) {
      m1_INLA_optim <- inla.rerun(m1_INLA_optim)
    }
    if (original.VCV == T) {
      m1_INLA_original <- inla.rerun(m1_INLA_original)
    }
    rerun <- rerun+1
    }
  }

  prediction_no_phylo <- cbind(df,
                               f(m1_INLA_GLM$summary.fitted.values, inverse=T),
                               model = "Without phylogeny",
                               better_fit = ifelse(best_model == "With phylogeny","No","Yes"))
  prediction_no_phylo <- prediction_no_phylo[pos,]

  fe_sig <- ifelse(sign(m1_INLA_GLM$summary.fixed[,3])* sign(m1_INLA_GLM$summary.fixed[,5]) == 1, "Sig.","Insig.")
  fe_sig <- data.frame(predictor = rownames(m1_INLA_GLM$summary.fixed),Sig=fe_sig)
  prediction_no_phylo <- prediction_no_phylo %>% inner_join(fe_sig,"predictor")

  prediction <- rbind(prediction_no_phylo,prediction_phylo)
  prediction <- split(prediction,prediction$predictor)

  pfr <- switch(predictedfittedresponse,
           best_m = ifelse(is.null(best_m),NA,list(best_m$summary.fitted.values[-pos,"0.5quant"])),
           original_VCV = ifelse(is.null(m1_INLA_original),NA,list(m1_INLA_original$summary.fitted.values[-pos,"0.5quant"])),
           no_phylo = ifelse(is.null(m1_INLA_GLM),NA,list(m1_INLA_GLMl$summary.fitted.values[-pos,"0.5quant"])),
           NA)

  pfr <- unlist(pfr)

  if (optim.lambda == T) {
    best_model <- m1_INLA_optim
    var_comp <- unlist(get_variance_INLA(best_model,m1_INLA_int_only))
  } else {
    best_model <- m1_INLA_original
    var_comp <- unlist(get_variance_INLA(best_model,m1_INLA_original_int_only))
  }

  R2 <- get_R2(var_comp)

  output <- list(best_model = best_model,
                 best_model_name = best_model_name,
                 predictedfittedresponse = pfr,
                 wAIC = wAIC_all,
                 optimized_lambda = return_lambda,
                 initial_formula = formula,
                 without_phylo_model = m1_INLA_GLM$summary.fixed,
                 optimized_model = m_optim_result,
                 original_VCV_model = m_original_result,
                 optimized_phylo_matrix = prec.mat.INLA,
                 prediction = prediction,
                 variance_component = var_comp,
                 R2 = R2,
                 init_lambda = ifelse(optim.lambda==T,grid_result$minlevels,NA)
                 )

  return(output)
}
