simu_spatial_multi_y <- function(X, beta_true, method_simulate = "N",
                         Sigma, GAMMA, RHO, W, nu = NULL) {
  
  # packages needed
  require("Matrix")
  
  # verification
  stopifnot(method_simulate %in% c("N", "IT"))
  
  # verification
  if (method_simulate != "N" & is.null(nu))
    stop("nu_D must be given when method_simulate is IT or UT")
  
  # initialization
  n <- nrow(X) # number of observations
  L <- ncol(beta_true) # the dimension of Y
  p <- nrow(beta_true) # the size of X
  
  # verification
  stopifnot(ncol(X) == p)
  stopifnot(ncol(GAMMA) == nrow(GAMMA) && nrow(GAMMA) == L)
  stopifnot(ncol(RHO) == nrow(RHO) && nrow(RHO) == L)
  
  # If method is Student, calculate the scatter matrix
  if (method_simulate != "N") {
    Sig2SD <- (nu - 2)/nu
    Scatter <- Sig2SD * Sigma
  }
  
  # the mean of the espilon
  Mu <- rep(0, nrow(Sigma))
  
  # Generate the residuals
  if (method_simulate == "IT") { # DGP = IT
    epsi <- mvnfast::rmvt(n, mu = Mu, sigma = Scatter, df = nu)
  } else { # DGP = N
      epsi <- MASS::mvrnorm(n, Mu, Sigma =  Sigma)
    }
  
  # rewrite the data in vectors 
  x_beta_epsi <- as.vector(X %*% beta_true + epsi)
  
  # Define the Y variable
  I_nL <- Diagonal(n * L)
  kro_GAMMA <- kronecker(GAMMA, Diagonal(n)) 
  kro_RHO <- kronecker(RHO, as(W, "Matrix"))    
  Y <- solve(I_nL - kro_GAMMA - kro_RHO, x_beta_epsi)
  
  # sortie
  return(matrix(Y, n, L))
}
