library(phytools)
library(tidyverse)
library(CPR)
library(ape)
library(DHARMa)

data(KSR)
data(KSR_MLtree)
data(KSR_EF)

m_LAI <- CPR_spaMM(LAI~Real.rich+corrMatrix(1|comp_id),
                   data=KSR_EF,
                   comm=KSR,
                   VCV_sp = vcv(KSR_MLtree),
                   method.spaMM = "REML")

get_R2(m_LAI$best_model) #get R2 for best model
get_R2(m_LAI$original_VCV_model) #get R2 for model based on the provided phylogenetic covariance matrix between species
####

m_bugs <- CPR_spaMM(bugs~Real.rich+corrMatrix(1|comp_id),
                   data=KSR_EF,
                   comm=KSR,
                   VCV_sp = vcv(KSR_MLtree),
                   family="negbin2")

# for the full model with lambda = 1
phy_comm_M <- get_comm_pair_r(KSR,vcv(KSR_MLtree))$corM

m_bugs_int <- fitme(bugs~1+corrMatrix(1|comp_id),corrMatrix=as_precision(phy_comm_M),data=KSR_EF,family="negbin2")

get_R2(m_bugs$original_VCV_model,m_bugs_int)

# for the full model with optimized lambda (also the best model in this case)
optim_VCV_sp <- vcv(KSR_MLtree)*m_LAI$optimized_lambda
diag(optim_VCV_sp) <- diag(vcv(KSR_MLtree))
phy_comm_M <- get_comm_pair_r(KSR,optim_VCV_sp)$corM

m_bugs_int <- fitme(bugs~1+corrMatrix(1|comp_id),corrMatrix=as_precision(phy_comm_M),data=KSR_EF,family="negbin2")

get_R2(m_bugs$best_model,m_bugs_int)
