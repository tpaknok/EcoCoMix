get_comm_pair_r <- function (comm,
                             VCV_sp,
                             comm_kronecker = NULL,
                             force.PD = F) {
  require(Matrix)
  require(matrixcalc)

  if (sum(colnames(comm) %in% rownames(VCV_sp)) != ncol(comm)) {
    stop("Inconsistent species name between species covariance matrix and community data matrix")
  }

  comm <- comm[,match(colnames(comm),colnames(VCV_sp))] #rearrange

  if (is.null(comm_kronecker)) {
  product <- as.data.frame(kronecker(as.matrix(comm),as.matrix(comm)))
  } else {
    product <- comm_kronecker
  }

  time1 <- Sys.time()
  row_op <- cbind(expand.grid(1:nrow(comm),1:nrow(comm)))
  row_op$cov <- as.matrix(product) %*% c(VCV_sp) #covariance
  covM <- matrix(row_op$cov,nrow(comm),nrow(comm))
  corM <- cov2cor(covM)
  time2 <- Sys.time()
  time2-time1
  #comm_var <- product[row_op$Var1==row_op$Var2,]
  #comm_var <- as.matrix(comm_var) %*% c(VCV_sp)
  #rownames(comm_var) <- 1:nrow(comm_var)
  #row_op$numerator1 <- comm_var[match(row_op$Var1,rownames(comm_var))]
  #row_op$numerator2 <- comm_var[match(row_op$Var2,rownames(comm_var))]

  #row_op$cor <- row_op$cov/sqrt(row_op$numerator1*row_op$numerator2)
  #corM <- matrix(row_op$cor,nrow(comm),nrow(comm))

  ### Calculate species turnover component
  # VCV_sp0 <- VCV_sp*0
  # diag(VCV_sp0) <- diag(VCV_sp)
  #
  # row_op0 <- cbind(expand.grid(1:nrow(comm),1:nrow(comm)))
  # row_op0$cov <- as.matrix(product) %*% c(VCV_sp0) #covariance
  # comm_var <- product[row_op$Var1==row_op$Var2,]
  # comm_var <- as.matrix(comm_var) %*% c(VCV_sp0)
  # rownames(comm_var) <- 1:nrow(comm_var)
  # #row_op$numerator1 <- comm_var[match(row_op$Var1,rownames(comm_var))]
  # #row_op$numerator2 <- comm_var[match(row_op$Var2,rownames(comm_var))]
  #
  # #row_op$cor <- row_op$cov/sqrt(row_op$numerator1*row_op$numerator2)
  # #corM <- matrix(row_op$cor,nrow(comm),nrow(comm))
  # covM0 <- matrix(row_op0$cov,nrow(comm),nrow(comm))
  # covMphylo <- covM-covM0

  if (force.PD == T & is.positive.definite(round(corM,5)) == F) { #loss of signifiance. need to round the matrix
    corM <-as.matrix(nearPD(corM,corr=T,keepDiag=T,maxit=100000)$mat)
  }

  result <- list(corM = corM,
                 covM = covM,
                 comm_kronecker = product)
  return(result)
}
