estimate_spatial_multi_N <- function(Y, X, W, GAMMA_esti = F) {
  
  # initialization
  L <- ncol(Y) # the dimension of Y
  n <- nrow(Y) # number of observations
  p <- ncol(X) # the size of X
  res_beta <- matrix(0, p, L)
  hat_Uk <- matrix(0, n, L)
  RHO <- matrix(0, L, L)
  GAMMA <- matrix(0, L, L)
  
  # verification
  stopifnot(nrow(X) == n)
  stopifnot(nrow(W) == n, ncol(W) == n)
  # check if the 1st column is the constant
  if (p > 1) {
    if (! all(X[, 1] == 1)) {
      X <- cbind(1, X)
    }
  }
  
  # initialisation
  W_X <- W %*% X[, -1]
  W2_X <- W %*% W_X
  H_n <- cbind(X, W_X, W2_X)
  W_Y <- W %*% Y
  
  P_H <- H_n %*% chol2inv(qr(H_n)$qr) %*% t(H_n)
  Z_tilde <- cbind(X, P_H %*% W %*% Y)
  
  if (GAMMA_esti) {
    Z_tilde <- cbind(Z_tilde, P_H %*% Y)
  } else {
    cste <- chol2inv(qr(Z_tilde)$qr) %*% t(Z_tilde)
  }
  
  for (k in 1:L) {
    if (GAMMA_esti) {
      cste <- chol2inv(qr(Z_tilde[, -(L + p + k)])$qr) %*% t(Z_tilde[, -(L + p + k)])  
    }
    res_k <- cste %*% Y[, k] 
    res_beta[, k] <- res_k[1:p]
    RHO[k, ] <- res_k[(p + 1):(p + L)]   
    hat_Uk[, k] <-  Y[, k] - X %*% res_beta[, k] - W_Y %*% RHO[k, ]
    if (GAMMA_esti) { 
      GAMMA[k, (1:L)[-k]] <- res_k[(p + L + 1):(p + 2 * L - 1)]
      hat_Uk[, k] <-  hat_Uk[, k] - Y %*% GAMMA[k, ] 
    }
  }
  
  SIGMA <- t(hat_Uk) %*% hat_Uk / n
  
  return(list(res_beta = res_beta,
              GAMMA = GAMMA,
              RHO = RHO,
              SIGMA = SIGMA))
}
