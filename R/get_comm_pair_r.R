get_comm_pair_r <- function (comm,V,force.PD = T) {
  require(Matrix)
  require(matrixcalc)
  if (sum(colnames(comm) %in% rownames(V)) != ncol(comm)) {
    stop("Inconsistent species name between species covariance matrix and community data matrix")
  }

  comm <- comm[,match(colnames(comm),colnames(V))] #rearrange
  product <- as.data.frame(kronecker(as.matrix(comm),as.matrix(comm)))
  row_op <- cbind(expand.grid(1:nrow(comm),1:nrow(comm)))

  row_op$cov <- as.matrix(product) %*% c(V) #covariance
  comm_var <- product[row_op$Var1==row_op$Var2,]
  comm_var <- as.matrix(comm_var) %*% c(V)
  rownames(comm_var) <- 1:nrow(comm_var)
  row_op$numerator1 <- comm_var[match(row_op$Var1,rownames(comm_var))]
  row_op$numerator2 <- comm_var[match(row_op$Var2,rownames(comm_var))]

  row_op$cor <- row_op$cov/sqrt(row_op$numerator1*row_op$numerator2)
  corM <- matrix(row_op$cor,nrow(comm),nrow(comm))

  if (force.PD == T & is.positive.definite(round(corM,5)) == F) { #loss of signifiance. need to round the matrix
    corM <-as.matrix(nearPD(corM,corr=T,keepDiag=T,maxit=100000)$mat)
  }

  return(corM)
}
