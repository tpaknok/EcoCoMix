get_comm_pair_r <- function (comm,
                             VCV_sp,
                             comm_kronecker = NULL,
                             force.PD = FALSE) {
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
  corM <- cov2cor(covM) #correlation matrix

  if (force.PD == T & is.positive.definite(round(corM,5)) == F) { #could be due to loss of significance. need to round the matrix
    corM <-as.matrix(nearPD(corM,corr=T,keepDiag=T,maxit=100000)$mat)
  }

  result <- list(corM = corM,
                 covM = covM,
                 comm_kronecker = product)
  return(result)
}
