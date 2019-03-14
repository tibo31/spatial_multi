simu_jrfm <- function(core, X, beta_true, Sigma, B, eps_precision) {
  
  # initialization 
  n <- nrow(X)
  nb_beta <- prod(dim(beta_true))
  nb_sigma <- dim(Sigma)[1] + 1
  crossprod_X <- crossprod(X)
  
  # table of estimate 
  table_simu_beta <- matrix(0, nb_beta * 6, 7) 
  table_simu_sigma <- matrix(0, nb_sigma * 6, 7) 
  UT_UT_beta <- numeric(nb_beta)
  UT_UT_sigma <- numeric(nb_sigma)
  table_simu_beta_square <- matrix(0, nb_beta * 6, 7) 
  table_simu_sigma_square <- matrix(0, nb_sigma * 6, 7) 
  UT_UT_beta_square <- numeric(nb_beta)
  UT_UT_sigma_square <- numeric(nb_sigma)
  
  for (b in 1:B) {
    # simulate the data
    for (i in 1:7) {
      # if i = 1,  method_simulate = "N"
      # if i = 2, method_simulate = "UT", nu_D = 3
      # if i = 3, method_simulate = "UT", nu_D = 4
      # if i = 4, method_simulate = "UT", nu_D = 5
      # if i = 5, method_simulate = "IT", nu_D = 3
      # if i = 6, method_simulate = "IT", nu_D = 4
      # if i = 7, method_simulate = "IT", nu_D = 5
      if (i == 1) {
        method_simulate <- "N"
        nu_D = NULL
      } else {
        if (i %in% c(2, 3, 4)) {
          method_simulate <- "UT"
          nu_D <- i + 1
        } else {
          method_simulate <- "IT"
          nu_D <- i - 2
        }
      }
      
      # simulation
      set.seed((core - 1) * B + b) 
      y <- simu_multi_y(X = X, beta_true = beta_true, 
                        method_simulate = method_simulate, Sigma = Sigma, 
                        nu = nu_D) 
      
      # estimate with type = "N"
      hat_beta_N <- solve(crossprod_X, t(X) %*% y)
      residuals_N <- (y - X %*% hat_beta_N) 
      hat_Sigma_N <- t(residuals_N) %*% residuals_N / n
      s1_N <- hat_Sigma_N[1, 1]
      s2_N <- hat_Sigma_N[2, 2]
      rho_N <- hat_Sigma_N[1, 2]/(sqrt(s1_N * s2_N))
      sigma_N <- c(rho_N, s1_N, s2_N)
      # storage
      table_simu_beta[1:nb_beta, i] <- table_simu_beta[1:nb_beta, i] + 
        hat_beta_N
      table_simu_sigma[1:nb_sigma, i] <- table_simu_sigma[1:nb_sigma, i] + 
        sigma_N   
      table_simu_beta_square[1:nb_beta, i] <- table_simu_beta_square[1:nb_beta, i] + 
        hat_beta_N^2
      table_simu_sigma_square[1:nb_sigma, i] <- table_simu_sigma_square[1:nb_sigma, i] + 
        sigma_N^2   
      
      # estimate with type = "U" only when i = 2 
      if (i == 2) {
        hat_beta_U <- hat_beta_N 
        residuals_U <- residuals_N
        hat_Sigma_U <- hat_Sigma_N  
        hat_Sigma_U <- nu_D/(nu_D - 2) * hat_Sigma_U
        s1_U <- hat_Sigma_U[1, 1]
        s2_U <- hat_Sigma_U[2, 2]
        rho_U <- hat_Sigma_U[1, 2]/(sqrt(s1_U * s2_U))
        sigma_U <- c(rho_U, s1_U, s2_U)
        # storage
        UT_UT_beta <- UT_UT_beta + as.vector(hat_beta_U)
        UT_UT_sigma <- UT_UT_sigma + sigma_U
        UT_UT_beta_square <- UT_UT_beta_square + as.vector(hat_beta_U^2)
        UT_UT_sigma_square <- UT_UT_sigma_square + sigma_U^2
      }
      
      # estimate with type = "IT"
      for (k in 1:5) {
        nu_D_estimate <- switch(k, "1" = 3,
                                "2" = 4,
                                "3" = 5,
                                "4" = 10,
                                "5" = 20)
        res_IT <- tryCatch(estimate_IT(y, X, nu = nu_D_estimate, 
                              eps_precision = eps_precision), error = function(e) "error")
        if (length(res_IT) == 1) {
          # if problem of inversion of the matrix Sigma, replace NA by average mean
          if (b != 1) {
            hat_beta_IT <- table_simu_beta[(k * nb_beta + 1):((k + 1) * nb_beta), i] / (b - 1)
            sigma_IT <- table_simu_sigma[(k * nb_sigma + 1):((k + 1) * nb_sigma), i] / (b - 1)
          } else {
            hat_beta_IT <- 0
            sigma_IT <- 0
          }
        } else {
          hat_beta_IT <- res_IT$hat_beta
          hat_Sigma_IT <- res_IT$hat_Sigma
          s1_IT <- hat_Sigma_IT[1, 1]
          s2_IT <- hat_Sigma_IT[2, 2]
          rho_IT <- hat_Sigma_IT[1, 2]/(sqrt(s1_IT * s2_IT))
          sigma_IT <- c(rho_IT, s1_IT, s2_IT)
        }
        
        # storage
        table_simu_beta[(k * nb_beta + 1):((k + 1) * nb_beta), i] <- 
          table_simu_beta[(k * nb_beta + 1):((k + 1) * nb_beta), i] + hat_beta_IT
        table_simu_sigma[(k * nb_sigma + 1):((k + 1) * nb_sigma), i] <- 
          table_simu_sigma[(k * nb_sigma + 1):((k + 1) * nb_sigma), i] + sigma_IT  
        table_simu_beta_square[(k * nb_beta + 1):((k + 1) * nb_beta), i] <- 
          table_simu_beta_square[(k * nb_beta + 1):((k + 1) * nb_beta), i] + hat_beta_IT^2
        table_simu_sigma_square[(k * nb_sigma + 1):((k + 1) * nb_sigma), i] <- 
          table_simu_sigma_square[(k * nb_sigma + 1):((k + 1) * nb_sigma), i] + sigma_IT^2  
      }
    }
  }
  return(list(
    table_simu_beta = table_simu_beta, 
    table_simu_sigma = table_simu_sigma,
    UT_UT_beta = UT_UT_beta,
    UT_UT_sigma = UT_UT_sigma,
    table_simu_beta_square = table_simu_beta_square, 
    table_simu_sigma_square = table_simu_sigma_square,
    UT_UT_beta_square = UT_UT_beta_square,
    UT_UT_sigma_square = UT_UT_sigma_square
    ))
}
