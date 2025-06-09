## Function Details

### **Maximum Likelihood Estimation (MLE) and Expectation-Maximization (EM)**
# The function uses the **M-step** and **EM algorithm** to iteratively estimate model parameters, including SNP effect sizes (`alpha`) and noise variances (`sigma1`, `sigma2`, `sigmau`).

# **M-step**: Optimizes the parameters for gene expression (`y`) and phenotype (`z`) given the current estimates of the SNP effects and noise levels.
# **EM algorithm**: Iteratively updates the parameter estimates using the M-step until convergence.

### **Model Estimation and Cross-validation**
# The function supports k-fold cross-validation to tune regularization parameters (`lambda`) for different models, such as Lasso and Elastic Net. It calculates the mean squared error (MSE) to evaluate model performance.

# - `k`: SNP loadings, defining how many SNPs affect each gene.
tips <- function(simulated_data, k, sigma1_init = 0.1, sigma2_init = 0.1, sigmau_init = 0.1, truealpha_init) {
  library(foreach)
  library(glmnet)
  library(MASS)
  library(parallel)

  # Extract variables
  wg1 <- simulated_data$wg1
  wg2 <- simulated_data$wg2
  u <- simulated_data$u_init
  y <- simulated_data$y
  z <- simulated_data$z
  e2 <- z - wg2 %*% u %*% truealpha_init

  n1 <- nrow(wg1)
  n2 <- nrow(wg2)
  m1 <- ncol(y)
  p1 <- nrow(u)

  M_step <- function(old, lambda, a, wg1, y, wg2, u, truealpha, e2, m1, n1, n2, p1, k) {
    sigma1 <- old[1]
    sigma2 <- old[2]
    sigmau <- old[3]
    alpha <- old[-(1:3)]

    results <- mclapply(1:m1, function(i) {
      wg1t_wg1 <- t(wg1) %*% wg1
      wg2t_wg2 <- t(wg2) %*% wg2
      Su <- diag(1, p1)
      alpha_i <- alpha[i]

      Sigma_u <- solve((1/sigma1) * wg1t_wg1 + (alpha_i^2 / sigma2) * wg2t_wg2 + (1/sigmau) * Su)
      mu_u <- Sigma_u %*% ((1/sigma1) * t(wg1) %*% y[,i] + (alpha_i / sigma2) * t(wg2) %*% (wg2 %*% u[,i] * truealpha[i] + e2))

      si <- 2 * t(wg2 %*% u[,i] * truealpha[i] + e2) %*% wg2 %*% mu_u - a * lambda * sign(alpha_i)
      thresh <- (1 - a) * lambda * k[i]

      alpha_est <- ifelse(abs(si) <= thresh | abs(si) <= a * lambda, 0,
                          (1 - thresh / abs(si)) *
                          (2 * t(mu_u) %*% t(wg2) %*% wg2 %*% mu_u +
                           2 * sum(diag(wg2 %*% Sigma_u %*% t(wg2))))^(-1) * si)

      Ez <- wg2 %*% u[,i] * truealpha[i] + e2

      list(
        E1 = sum((y[,i] - wg1 %*% mu_u)^2) + sum(diag(wg1 %*% Sigma_u %*% t(wg1))),
        E2 = sum((Ez - alpha_est * wg2 %*% mu_u)^2) + sum(diag(alpha_est^2 * wg2 %*% Sigma_u %*% t(wg2))),
        E3 = sum(mu_u^2) + sum(diag(Sigma_u)),
        alpha = alpha_est
      )
    }, mc.cores = parallel::detectCores())

    sigma1 <- mean(sapply(results, function(x) x$E1)) / n1
    sigma2 <- mean(sapply(results, function(x) x$E2)) / n2
    sigmau <- mean(sapply(results, function(x) x$E3)) / p1

    loglik <- sum(sapply(results, function(x) -0.5 * (log(2 * pi * sigma1) + x$E1 / sigma1))) +
              sum(sapply(results, function(x) -0.5 * (log(2 * pi * sigma2) + x$E2 / sigma2))) +
              sum(sapply(results, function(x) -0.5 * (log(2 * pi * sigmau) + x$E3 / sigmau)))

    return(list(est = c(sigma1, sigma2, sigmau, sapply(results, `[[`, "alpha")), lik = loglik))
  }

  EM <- function(init, lambda, a, ...) {
    theta <- init
    for (iter in 1:20) {
      step <- M_step(theta, lambda, a, ...)
      if (mean(abs(theta - step$est)) < 1e-3) break
      theta <- step$est
    }
    return(c(theta, lambda, step$lik))
  }

  msecal <- function(wg2, u, mle_mat, truealpha, e2) {
    sapply(1:nrow(mle_mat), function(i) {
      alpha_hat <- mle_mat[i, 4:(3 + ncol(u))]
      z_hat <- wg2 %*% u %*% alpha_hat + e2
      z_true <- wg2 %*% u %*% truealpha + e2
      mean((z_true - z_hat)^2)
    })
  }

  # Choose lambda and a for "small" setting
  a_vals <- c(0, 0.25, 0.75, rep(1, 12))
  lambda_vals <- c(0, 10, 15, 20, seq(30, 200, length.out = 9), 400, 700)

  mle_results <- foreach(i = seq_along(lambda_vals), .combine = rbind) %do% {
    EM(c(sigma1_init, sigma2_init, sigmau_init, truealpha_init), lambda_vals[i], a_vals[i],
       wg1 = wg1, y = y, wg2 = wg2, u = u, truealpha = truealpha_init, e2 = e2,
       m1 = m1, n1 = n1, n2 = n2, p1 = p1, k = k)
  }

  mse <- msecal(wg2, u, mle_results, truealpha_init, e2)
  list(mle = mle_results, mse = mse)
}
