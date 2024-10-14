## Function Details

### **Maximum Likelihood Estimation (MLE) and Expectation-Maximization (EM)**
# The function uses the **M-step** and **EM algorithm** to iteratively estimate model parameters, including SNP effect sizes (`alpha`) and noise variances (`sigma1`, `sigma2`, `sigmau`).

# **M-step**: Optimizes the parameters for gene expression (`y`) and phenotype (`z`) given the current estimates of the SNP effects and noise levels.
# **EM algorithm**: Iteratively updates the parameter estimates using the M-step until convergence.

### **Model Estimation and Cross-validation**
# The function supports k-fold cross-validation to tune regularization parameters (`lambda`) for different models, such as Lasso and Elastic Net. It calculates the mean squared error (MSE) to evaluate model performance.

# - `k`: SNP loadings, defining how many SNPs affect each gene.
tips <- function(simulated_data, k, sigma1_init = 0.1, sigma2_init = 0.1, sigmau_init = 0.1, truealpha_init) {
  
  # Extract variables from the simulated data
  wg1 <- simulated_data$wg1
  wg2 <- simulated_data$wg2
  u_init <- simulated_data$u_init
  y <- simulated_data$y
  z <- simulated_data$z
  
  # Derive dimensions from the data
  n1 <- dim(wg1)[1]  # Number of samples in wg1
  n2 <- dim(wg2)[1]  # Number of samples in wg2
  m1 <- dim(y)[2]    # Number of genes
  p1 <- dim(u_init)[1]  # Number of SNPs
  
  # Load required libraries
  library(foreach)
  library(glmnet)
  library(MASS)
  library(parallel)
  
  # Define helper functions (M-step and EM algorithm)
  
  M_step <- function(old, lambda, a_para, wg1, y, wg2, u, truealpha, e2, m1, n1, n2, p1, k) {
    # Extract parameters
    sigma1 <- old[1]
    sigma2 <- old[2]
    sigmau <- old[3]
    alpha_g <- old[-c(1:3)]
    a <- a_para
    abc <-  (1 - a) * lambda * k
    ab <- a * lambda
    wg1_product <- t(wg1) %*% wg1
    wg2_product <- t(wg2) %*% wg2
    sigma_ui_constant <- diag(1, p1, p1)
    
    # Parallel computation over m1 genes
    results <- mclapply(1:m1, function(i) {
      sigma_ui_i <- solve(1 / sigma1 * wg1_product + (alpha_g[i]^2) / sigma2 * wg2_product + 1 / sigmau * sigma_ui_constant)
      mu_ui_i <- sigma_ui_i %*% (1 / sigma1 * t(wg1) %*% y[, i] + alpha_g[i] / sigma2 * t(wg2) %*% (wg2 %*% u[, i] * truealpha[i] + e2))
      
      si <- 2 * t(wg2 %*% u[, i] * truealpha[i] + e2) %*% wg2 %*% mu_ui_i - a * lambda * sign(alpha_g[i])
      abc <-  (1 - a) * lambda * k[i]
      ab <- a * lambda
      alphaest <- ifelse(abs(si) <= abc | abs(si) <= ab, 0, (1 - abc / (abs(si))) * 
                           (2 * t(mu_ui_i) %*% t(wg2) %*% wg2 %*% mu_ui_i + 2 * sum(diag(wg2 %*% sigma_ui_i %*% t(wg2))))^(-1) * si)
      newz <- wg2 %*% u[, i] * truealpha[i] + e2
      E1 <- t(y[, i]) %*% y[, i] - 2 * t(y[, i]) %*% wg1 %*% mu_ui_i + t(mu_ui_i) %*% t(wg1) %*% wg1 %*% mu_ui_i + sum(diag(t(wg1) %*% wg1 %*% sigma_ui_i))
      E2 <- t(newz) %*% newz - 2 * alphaest * t(newz) %*% wg2 %*% mu_ui_i + t(as.numeric(alphaest) * wg2 %*% mu_ui_i) %*% (as.numeric(alphaest) * wg2 %*% mu_ui_i) + sum(diag(as.numeric(alphaest)^2 * wg2 %*% sigma_ui_i %*% t(wg2)))
      E3 <- t(mu_ui_i) %*% mu_ui_i + sum(diag(sigma_ui_i))
      list(E1 = E1, E2 = E2, E3 = E3, alphaest = alphaest, mu_ui_i = mu_ui_i)
    }, mc.cores = parallel::detectCores())
    
    sigma1est <- 1 / (m1 * n1) * sum(sapply(results, function(x) x$E1))
    sigma2est <- 1 / (m1 * n2) * sum(sapply(results, function(x) x$E2))
    sigmauest <- 1 / (m1 * p1) * sum(sapply(results, function(x) x$E3))
    log_likelihood_y <- -0.5 * n1 * m1 * log(2 * pi * sigma1est) - 0.5 * sum(sapply(results, function(x) x$E1)) / sigma1est
    log_likelihood_z <- -0.5 * n2 * log(2 * pi * sigma2est) - 0.5 * sum(sapply(results, function(x) x$E2)) / sigma2est
    log_likelihood_u <- -0.5 * p1 * m1 * log(2 * pi * sigmauest) - 0.5 * sum(sapply(results, function(x) x$E3)) / sigmauest
    log_likelihood <- log_likelihood_y + log_likelihood_z + log_likelihood_u
    # Print estimated parameters for debugging
    # print(c(sigma1est, sigma2est, sigmauest))
    return(list(est = c(sigma1est, sigma2est, sigmauest, unlist(sapply(results, function(x) x$alphaest))),
                lik = log_likelihood))
  }
  
  EM <- function(init, lambda, a_para, wg1, y, wg2, u, truealpha, e2, m1, n1, n2, 
                 p1, k, max_iter = 10, tol = 0.001) {
    # Initialize the parameters
    theta_old <- init - 0.0001
    iter <- 1
    # Initialize the previous difference to an arbitrarily large value
    prev_diff <- Inf
    # Loop over the maximum number of iterations
    for (i in 1:max_iter) {
      step_results <- M_step(old = theta_old, lambda = lambda, a_para = a_para, 
                             wg1 = wg1, y = y, wg2 = wg2, u = u, 
                             truealpha = truealpha, e2 = e2, m1 = m1, n1 = n1, n2 = n2, 
                             p1 = p1, k = k)
      # Extract the new parameter estimates and log-likelihood
      theta_new <- step_results$est
      log_likelihood_new <- step_results$lik
      # Compute the current difference
      current_diff <- mean(abs(theta_new - theta_old))
      # Print progress
      # cat("Iteration:", iter, "parameter difference:", current_diff, "\n")
      # Check for convergence using the difference between the current and previous differences
      if (abs(current_diff - prev_diff) < tol) break
      # Update for the next iteration
      theta_old <- theta_new
      prev_diff <- current_diff
      iter <- iter + 1
      # print(prev_diff)
    }
    # Return the estimated parameters
    return(list(final_parameters = c(theta_new, lambda, log_likelihood_new)))
  }
  
  # Function to calculate MSE
  msecal <- function(wg2, mle, lam_values) {
    mse <- rep(0, length(lam_values))
    for (i in 1:length(lam_values)) {
      alphaest0 <- mle[i, 4:(m1 + 3)]
      z_est0 <- wg2 %*% u %*% alphaest0  + e2
      z_true <- wg2 %*% u %*% truealpha + e2
      mse[i] <- mean((z_true - z_est0)^2)
    }
    return(mse)
  }
  
  # Set lambda and a_para values based on size
  if (size == "small") {
    a_mle <-  c(0, 0.25, 0.75, rep(1, 12))
    lam_values_l0 <- c(0, 10, 15, 20, seq(30, 200, length = 9), 400, 700)
    lam_values_l1 <- c(0, 5, 10, 14, seq(17, 200, length = 11))
    lam_values_l2 <- c(0, 15, 30, seq(40, 250, length = 10), 350, 550)
  } else if (size == "large") {
    a_mle <-  c(0, 0.5, rep(1, 13))
    lam_values_l0 <- c(0, 5, 15, seq(20, 300, length = 10), 400, 500)
    lam_values_l1 <- c(0, 3, 7, seq(15, 200, length = 10), 300, 400)
    lam_values_l2 <- c(0, 5, 10, seq(20, 300, length = 10), 400, 500)
  } else if (size == "hcsmall") {
    a_mle <- c(0.00, 0.25, 0.75, rep(1, 12)) 
    lam_values_l0 <- c(0, 1, 2, 3, seq(3.5, 40, length = 10), 50)
    lam_values_l1 <- c(0, 1, 2, 3, 5, seq(6, 25, length = 9), 50)
    lam_values_l2 <- c(0, 3, 6, 10, seq(12, 50, length = 10), 100)
  } else if (size == "htsmall") {
    a_mle <- c(0.00, 0.50, 0.75, rep(1, 12))
    lam_values_l0 <- c(0, 3, 7, 10, seq(15, 60, length = 10), 100)
    lam_values_l1 <- c(seq(0, 10, length = 11), 12, 15, 20, 200)
    lam_values_l2 <-  c(0, 3, 6, seq(8, 33, length = 10), 50, 100)
  }
  
  # Run EM algorithm for different lambda and a_para values
  mle_l0 <- foreach(i = seq_along(lam_values_l0), .combine = rbind) %do% {
    lambda <- lam_values_l0[i]
    a_para <- a_mle[i]
    EM(init = c(sigma1, sigma2, sigmau, truealpha), 
       lambda = lambda, a_para = a_para, 
       wg1, y, wg2, u, truealpha, e2, m1, n1, n2, p1, k)$final_parameters
  }
  
  mle_l1 <- foreach(i = lam_values_l1, .combine = rbind) %do% 
    EM(init = c(sigma1, sigma2, sigmau, truealpha), 
       lambda = i, a_para = 0,  wg1, y, wg2, u, truealpha, e2, m1, n1, n2, p1, k)$final_parameters
  
  mle_l2 <-  foreach(i = lam_values_l2, .combine = rbind) %do% 
    EM(init = c(sigma1, sigma2, sigmau, truealpha), 
       lambda = i, a_para = 1,  wg1, y, wg2, u, truealpha, e2, m1, n1, n2, p1, k)$final_parameters
  
  # Calculate MSE for each method
  mse_l0 <- msecal(wg2, mle_l0, lam_values_l0)
  mse_l1 <- msecal(wg2, mle_l1, lam_values_l1)
  mse_l2 <- msecal(wg2, mle_l2, lam_values_l2)
  
  # Compile results into a list
  result_list <- list(
    mle_l0 = mle_l0,
    mle_l1 = mle_l1,
    mle_l2 = mle_l2,
    mse_l0 = mse_l0,
    mse_l1 = mse_l1,
    mse_l2 = mse_l2)
  return(result_list)
}
