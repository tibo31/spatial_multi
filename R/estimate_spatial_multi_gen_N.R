estimate_spatial_multi_gen_N <- function(Y, X, W, 
                                         ind_beta = matrix(T, ncol(X), ncol(Y)),
                                         ind_RHO = matrix(T, ncol(Y), ncol(Y)),
                                         ind_GAMMA = matrix(F, ncol(Y), ncol(Y))) {
  
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
  Z_tilde <- cbind(X, P_H %*% W %*% Y, P_H %*% Y)

  for (k in 1:L) {
    ind_beta_l <- ind_beta[, k]
    p_l <- length(which(ind_beta_l))
    ind_rho_l <- ind_RHO[k, ]
    rho_l <- length(which(ind_rho_l))
    ind_GAMMA_l <- ind_GAMMA[, k]
    GAMMA_l <- length(which(ind_GAMMA_l))
    ind_l <- c(ind_beta_l, ind_rho_l, ind_GAMMA_l) 
    Z_tilde_l <- Z_tilde[, ind_l]
    cste <- chol2inv(qr(Z_tilde_l)$qr) %*% t(Z_tilde_l)
    res_k <- cste %*% Y[, k] 
    res_beta[ind_beta_l, k] <- res_k[1:p_l]
    RHO[k, ind_rho_l] <- res_k[(p_l + 1):(p_l + rho_l)]
    GAMMA[k, ind_GAMMA_l] <- res_k[(p_l + rho_l + 1):(p_l + rho_l + GAMMA_l)]
    hat_Uk[, k] <-  Y[, k] - Y %*% GAMMA[k, ] - X %*% res_beta[, k] - W_Y %*% RHO[k, ]
  }
  
  SIGMA <- t(hat_Uk) %*% hat_Uk / n
  
  return(list(res_beta = res_beta,
              GAMMA = GAMMA,
              RHO = RHO,
              SIGMA = SIGMA))
}