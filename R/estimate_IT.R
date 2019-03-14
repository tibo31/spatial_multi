estimate_IT <- function(Y, X, nu, eps_precision = 10^(-8)) {
  
  # initialization
  L <- ncol(Y) # the dimension of Y
  n <- nrow(Y) # number of observations
  p <- ncol(X) # the size of X
  
  # verification
  stopifnot(nrow(X) == n)
  
  # initialization
  Sig2SE <- (nu - 2)/nu
  S2SigE <- nu/(nu - 2)
  
  # Transform X into a kronecker form
  XX <- matrix(0, L*n, L*p)
  # Transform X into a kronecker form
  ind_col1 <- 1
  ind_col2 <- p
  for (k in 1:L) {
    XX[, ind_col1:ind_col2] <- kronecker(X, diag(L)[, k])
    ind_col1 <- ind_col1 + p
    ind_col2 <- ind_col2 + p
  }
  
  # Transform Y into kronecker form
  YY <- as.vector(t(Y))
  
  # Initial estimates (k = 0)
  hat_beta0 <- solve(crossprod(XX), crossprod(XX, YY))
  hat_U0 <- Y - X%*% matrix(hat_beta0, p, L)
  hat_sigma0 <- t(hat_U0) %*% hat_U0 / n
  
  # Estimates at step k 
  hat_betak <- hat_beta0
  hat_Uk <- hat_U0
  hat_sigmak <- hat_sigma0
  
  k <- 1
  repeat {
    
    # hat_sigmak_inv <- solve(hat_sigmak)
    hat_sigmak_inv <- chol2inv(chol(hat_sigmak))
    temp_uk <- hat_Uk %*% hat_sigmak_inv
    temp_uk <- temp_uk * hat_Uk
    temp_uk <- apply(temp_uk, 1, sum) 
    wk1 <- (nu + L) / (nu - 2 + temp_uk)
    
    # use formula with Kronecker on X
    temp_k1 <- t(XX) %*% kronecker(diag(wk1), hat_sigmak_inv)
    # hat_betak1 <- solve(temp_k1 %*% XX, temp_k1 %*% YY)
    hat_betak1 <- chol2inv(chol(temp_k1 %*% XX)) %*% temp_k1 %*% YY    
    hat_Uk1 <-  Y - X %*% matrix(hat_betak1, p, L)
    hat_sigmak1 <- t(hat_Uk1) %*% diag(wk1) %*% hat_Uk1 / n
    
    # Check the stopping criterion
    if (sum((hat_betak - hat_betak1)^2) <= eps_precision*crossprod(hat_beta0) || 
        k > 1000) {
      # Save the results (estimates and tests)
      res <- list(
        hat_beta = hat_betak,
        hat_Sigma = hat_sigmak,
        residuals = hat_Uk
      )
      return(res)
    } else {
      k <- k + 1
      hat_betak <- hat_betak1
      hat_Uk <- hat_Uk1
      hat_sigmak <- hat_sigmak1
    }
  }
  print("The algorithm did not converge towards a solution")
}