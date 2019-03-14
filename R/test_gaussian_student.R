mahalo_gaussian_student <- function(Y, mu = NULL, Sigma = NULL,
                                  type = "Gaussian", 
                                  nu = NULL) {
  
  # initialisation 
  eps_precision <- 10^(-8)
  
  # verification
  stopifnot(type %in% c("Gaussian", "Student"))
  
  # verification 
  if (type == "Student" & is.null(nu)) {
    stop("nu must be given when type is Student")
  }
  
  # initialization
  L <- ncol(Y) # the dimension of Y
  n <- nrow(Y) # number of observations
  
  # Case Gaussian 
  if (type == "Gaussian") {
    if (is.null(mu)) {
      mu <- apply(Y, 2, mean)
      mu <- matrix(rep(mu, n), n, L, byrow = T)
    }
    
    if (is.null(Sigma)) {
      Sigma <- cov(Y)
    }
    
    dM <- crossprod(t(Y - mu), solve(Sigma, t(Y - mu)))
    # TN <- ks.test(diag(dM), "pchisq", L)
  }
  
  # Case Student   
  if ((is.null(mu) | is.null(Sigma)) & type == "Student") {
    X <- matrix(1, nrow = n, ncol = 1)
    S2SigE <- nu/(nu - 2)
    
    res_esti <- estimate_IT(Y, X, nu = nu, eps_precision = 10^(-8))
    
    hat_Sigma <- res_esti$hat_Sigma
    hat_Uk <- res_esti$residuals 
    dM <- crossprod(t(hat_Uk), solve(hat_Sigma, t(hat_Uk)))
    # TN <- ks.test(diag(1/L * S2SigE * dM), "pf", L, nu)
  } else {
    if (type == "Student") { # user gives the parameters
      dM <- crossprod(t(Y - mu), solve(Sigma, t(Y - mu)))
      # TN <- ks.test(diag(1/L * S2SigE * dM), "pf", L, nu)
    }
  }
  return(dM)
}