BEF_simulate_eval <- function(models, type="INLA") {
    without_phylo_model <- lapply(1:length(models),function(x) models[[x]]$without_phylo_model[2,])
    original_VCV_model <- lapply(1:length(models),function(x) models[[x]]$original_VCV_model[2,])
    optimized_model <- lapply(1:length(models),function(x) models[[x]]$optimized_model)

    NA_pos <- which(is.na(optimized_model))

    if (length(NA_pos) >= 1) {
    fill <- models[[1]]$without_phylo_model
    fill[!is.na(fill)] <- NA
    optimized_model[NA_pos] <- list(fill)
    }

    optimized_model <- lapply(1:length(models),function(x) optimized_model[[x]][2,])

    without_phylo_model <- as.data.frame(do.call(rbind,without_phylo_model))
    original_VCV_model <- as.data.frame(do.call(rbind,original_VCV_model))
    optimized_model <- as.data.frame(do.call(rbind,optimized_model))
    original_VCV_model$i <- 1:length(models)
    without_phylo_model$i <-  1:length(models)
    optimized_model$i <- 1:length(models)

    original_VCV_model$model <- "User supplied VCV GLS"
    without_phylo_model$model <-  "Without phylogeny"
    optimized_model$model <- "Optimized phylogeny"

    optim_lambda <- lapply(1:length(models),function(x) models[[x]]$optim_lambda)
    optim_lambda <- do.call(rbind,optim_lambda)

    original_VCV_model$optim_lambda <- NA
    without_phylo_model$optim_lambda <-  NA
    optimized_model$optim_lambda <- optim_lambda

    all_models <- rbind(original_VCV_model,without_phylo_model,optimized_model)

    summary_stat <- switch(type,
                           GLS = {all_models %>%
                                    group_by(model) %>%
                                    summarize(sig_result = sum(`p-value` < 0.05,na.rm=T),
                                              proportion = sum(`p-value` < 0.05,na.rm=T)/sum(`p-value` > 0,na.rm=T),
                                              min_slope = min(Value,na.rm=T),
                                              max_slope = max(Value,na.rm=T),
                                              mean_slope = mean(Value,na.rm=T)
                                    )
                           },
                           INLA = {all_models %>%
                               group_by(model) %>%
                               summarize(sig_result = sum(sign(`0.025quant`)*sign(`0.975quant`) == 1,na.rm=T),
                                         proportion = sum(sign(`0.025quant`)*sign(`0.975quant`) == 1,na.rm=T)/sum(`sd` > 0,na.rm=T),
                                         min_slope = min(`0.5quant`,na.rm=T),
                                         max_slope = max(`0.5quant`,na.rm=T),
                                         mean_slope = mean(`0.5quant`,na.rm=T))
                           }
                      )

    AIC <- switch(type,
                  GLS = do.call(rbind,lapply(1:length(models),function(x) models[[x]]$AIC)),
                  INLA = do.call(rbind,lapply(1:length(models),function(x) models[[x]]$wAIC))
    )

    optim_lambda <- data.frame(min_lambda = min(optimized_model$optim_lambda,na.rm=T),
                               max_lambda = max(optimized_model$optim_lambda,na.rm=T),
                               mean_lambda = mean(optimized_model$optim_lambda,na.rm=T),
                               n = length(na.omit(optimized_model$optim_lambda)))

    output <- list(summary_stat = summary_stat,
                   result = all_models,
                   AIC = AIC,
                   optim_lambda=optim_lambda,
                   type = type)

    return(output)

  }
