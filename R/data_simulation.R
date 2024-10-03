## Function Details

### 1. **Generating Genotype Matrix**
# The function creates SNP data (`wg1` and `wg2`) for two groups. The genotype is generated under different inheritance models:

# **Additive**: SNPs are modeled additively, assuming a linear combination of alleles.
# **Heterogeneous**: SNPs are binary but reduced when an allele appears more than once.
# **Recessive**: SNPs with a recessive model, where heterozygous variants are assigned a value of 0.
# **Compensatory**: A compensatory model where having two mutated alleles cancels out the genetic effect.

### 2. **Data Normalization**
# After generating the genotype matrices (`wg1` and `wg2`), the function normalizes the data, transforming it to fit within the range [0, 1].

### 3. **Simulating Gene Expression (`y`) and Phenotype (`z`)**
# The gene expression matrix `y` for group 1 and the phenotype vector `z` for group 2 are simulated using the genotype matrix and the specified true effect sizes (`truealpha`).

### 4. **Maximum Likelihood Estimation (MLE) and Expectation-Maximization (EM)**
# The function uses the **M-step** and **EM algorithm** to iteratively estimate model parameters, including SNP effect sizes (`alpha`) and noise variances (`sigma1`, `sigma2`, `sigmau`).

# **M-step**: Optimizes the parameters for gene expression (`y`) and phenotype (`z`) given the current estimates of the SNP effects and noise levels.
# **EM algorithm**: Iteratively updates the parameter estimates using the M-step until convergence.

### 5. **Model Estimation and Cross-validation**
# The function supports k-fold cross-validation to tune regularization parameters (`lambda`) for different models, such as Lasso and Elastic Net. It calculates the mean squared error (MSE) to evaluate model performance.

### 6. **Supported Models**
# The function simulates several models commonly used in genetics research:
# **Lasso**: Performs L1 regularization, useful for feature selection.
# **Elastic Net (EN)**: Combines L1 and L2 regularization for improved prediction accuracy.
# **Multivariate TWAS**: Models multiple gene expression effects on the phenotype.
# **CoMM**: A correlation minimization model to improve phenotype prediction.

### 7. **P-values and Likelihood Ratio Tests**
# For each gene, likelihood ratio tests are conducted to assess the significance of SNP effects. The function returns p-values for the association between gene expression and the phenotype using different methods such as Lasso, Elastic Net, and TWAS.



data_simu <- function(n1,n2,m1,p1,k,sigma1,sigma2,sigmau,truealpha,
                      size=c("small","large","hcsmall","htsmall"),
                      wg_str = c("additive","heterogeneous","recessive","compensatory")){
  # Function to standardize a matrix by column
  standardize <- function(mat) {
    col_means <- colMeans(mat)
    col_sds <- apply(mat, 2, sd)
    standardized_mat <- scale(mat, center = col_means, scale = col_sds)
    return(standardized_mat)
  }
  
  generatewg <- function(n1){
    sequ1 <- rbinom(n1,1,0.2) ## asssume MAF=0.5 or 0.2
    sequ2 <- rbinom(n1,1,0.2)
    allele <- colSums(rbind(sequ1,sequ2))
    return(allele)
  }
  if (wg_str == "additive"){
    wg1 <- foreach(i=1:p1, .combine=cbind, .packages="foreach") %do% generatewg(n1)
    wg2 <- foreach(i=1:p1, .combine=cbind, .packages="foreach") %do% generatewg(n2)
  }else if (wg_str == "heterogeneous"){
    wg1 <- foreach(i=1:p1, .combine=cbind, .packages="foreach") %do% generatewg(n1)
    wg2 <- foreach(i=1:p1, .combine=cbind, .packages="foreach") %do% generatewg(n2)
    wg1[wg1 == 2] <- 1
    wg2[wg2 == 2] <- 1
  }else if (wg_str == "recessive"){
    wg1 <- foreach(i=1:p1, .combine=cbind, .packages="foreach") %do% generatewg(n1)
    wg2 <- foreach(i=1:p1, .combine=cbind, .packages="foreach") %do% generatewg(n2)
    wg1[wg1 == 1] <- 0
    wg2[wg2 == 1] <- 0
  }else if (wg_str == "compensatory"){
    wg1 <- foreach(i=1:p1, .combine=cbind, .packages="foreach") %do% generatewg(n1)
    wg2 <- foreach(i=1:p1, .combine=cbind, .packages="foreach") %do% generatewg(n2)
    wg1[wg1 == 2] <- 0
    wg2[wg2 == 2] <- 0
  }

  normalize_matrix <- function(mat) {
    min_vals <- apply(mat, 2, min) # Find the minimum value of each column
    max_vals <- apply(mat, 2, max) # Find the maximum value of each column
    return((mat - min_vals) / (max_vals - min_vals))
  }
  
  wg1 <- normalize_matrix(wg1)
  wg2 <- normalize_matrix(wg2)
  # Standardize wg1 and wg2
  #wg1 <- standardize(wg1)
  #wg2 <- standardize(wg2)
  if (size == "small"){
    sigman=0.05
  }else if (size=="hcsmall"){
    sigman=0.005
  }else if (size=="htsmall"){
    sigman=0.005
  }else{
    sigman=0.05
  }
  # Initialize the matrix u
  u <- matrix(rnorm(p1 * m1, 0, sqrt(sigman)), p1, m1)
  for (i in 1:m1) {
    block_start <- (i - 1) * 5 + 1
    u[block_start:(block_start + 4), i] <- rnorm(5, 0, sqrt(sigmau))
  }

  e1 <- matrix(rnorm(n1 * m1, 0, sqrt(sigma1)), n1, m1)
  
  e2 <- rnorm(n2,0,sqrt(sigma2))
  
  y <- wg1%*%u +e1
  z <- wg2%*%u%*%truealpha + e2
  
  
  
  M_step <- function(old, lambda, a_para,wg1, y, wg2, u, truealpha, e2, m1, n1, n2,p1, k) {
    # old<- c(sigma1,sigma2,sigmau,truealpha-0.00001)
    sigma1 <- old[1]; sigma2 <- old[2]; sigmau <- old[3]; alpha_g <- old[-c(1:3)]
    a <- a_para
    abc <-  (1 - a) * lambda * k
    ab <- a * lambda 
    wg1_product <- t(wg1) %*% wg1
    wg2_product <- t(wg2) %*% wg2
    sigma_ui_constant <- diag(1, p1, p1)
    
    results <- mclapply(1:m1, function(i) {
      sigma_ui_i <- solve(1 / sigma1 * wg1_product + (alpha_g[i]^2) / sigma2 * wg2_product + 1 / sigmau * sigma_ui_constant)
      mu_ui_i <- sigma_ui_i %*% (1 / sigma1 * t(wg1) %*% y[, i] + alpha_g[i] / sigma2 * t(wg2) %*% (wg2 %*% u[, i] %*% truealpha[i] + e2))
      
      si <- 2 * t(wg2 %*% u[, i] %*% truealpha[i] + e2) %*% wg2 %*% mu_ui_i - a * lambda * sign(alpha_g[i])
      abc <-  (1 - a) * lambda * k[i]
      ab <- a * lambda 
      alphaest <- ifelse(abs(si) <= abc | abs(si) <= ab, 0,(1 - abc / (abs(si))) * 
                           (2 * t(mu_ui_i) %*% t(wg2) %*% wg2 %*% mu_ui_i + 2 * sum(diag(wg2 %*% sigma_ui_i %*% t(wg2))))^(-1) * si)
      newz <- wg2 %*% u[, i] * truealpha[i] + e2
      E1 <- t(y[, i]) %*% y[, i] - 2 * t(y[, i]) %*% wg1 %*% mu_ui_i + t(mu_ui_i) %*% t(wg1) %*% wg1 %*% mu_ui_i + sum(diag(t(wg1) %*% wg1 %*% sigma_ui_i))
      E2 <- t(newz) %*% newz - 2 * alphaest * t(newz) %*% wg2 %*% mu_ui_i + t(as.numeric(alphaest) * wg2 %*% mu_ui_i) %*% (as.numeric(alphaest) * wg2 %*% mu_ui_i) + sum(diag(as.numeric(alphaest)^2 * wg2 %*% sigma_ui_i %*% t(wg2)))
      E3 <- t(mu_ui_i) %*% mu_ui_i + sum(diag(sigma_ui_i))
      list(E1 = E1, E2 = E2, E3 = E3, alphaest = alphaest, mu_ui_i = mu_ui_i)
    }, mc.cores = 30)
    
    sigma1est <- 1 / (m1 * n1) * sum(sapply(results, function(x) x$E1))
    sigma2est <- 1 / (m1 * n2) * sum(sapply(results, function(x) x$E2))
    sigmauest <- 1 / (m1 * p1) * sum(sapply(results, function(x) x$E3))
    log_likelihood_y = -0.5 * n1*m1 * log(2 * pi * sigma1est) - 0.5 * sum(sapply(results, function(x) x$E1)) / sigma1est;
    log_likelihood_z = -0.5 * n2 * log(2 * pi * sigma2est) - 0.5 * sum(sapply(results, function(x) x$E2)) / sigma2est;
    log_likelihood_u = -0.5 * p1*m1 * log(2 * pi * sigmauest) - 0.5 * sum(sapply(results, function(x) x$E3)) / sigmauest;
    log_likelihood = log_likelihood_y + log_likelihood_z + log_likelihood_u
    print(c(sigma1est,sigma2est,sigmauest))
    return(list(est=c(sigma1est, sigma2est, sigmauest, unlist(sapply(results, function(x) x$alphaest)))
                ,lik = log_likelihood))
  }
  
  EM <- function(init, lambda, a_para, wg1, y, wg2, u, truealpha, e2, m1, n1, n2, 
                 p1, k,max_iter = 10, tol = 0.001) {
    # Initialize the parameters
    theta_old <- init-0.0001
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
      cat("Iteration:", iter, "parameter difference:", current_diff, "\n")
      # Check for convergence using the difference between the current and previous differences
      if (abs(current_diff - prev_diff) < tol) break
      
      # Update for the next iteration
      theta_old <- theta_new
      prev_diff <- current_diff
      iter <- iter + 1
      print(prev_diff)
    }
    
    # Return the estimated parameters
    return(list(final_parameters = c(theta_new, lambda,log_likelihood_new)))
  }
  k_fold <- function(kf,a_values,lam_values,y,wg1,e1,wg2,e2,u,m1,truealpha,p1,k){
    # Split data into K-folds
    fold_size_y <- dim(y)[1] / kf
    fold_size_wg2 <- dim(wg2)[1] / kf
    
    # Create a list of indices for each fold
    fold_indices <- lapply(1:kf, function(i) {
      list(y_indices = ((i-1)*fold_size_y + 1):(i*fold_size_y),
           wg2_indices = ((i-1)*fold_size_wg2 + 1):(i*fold_size_wg2))
    })
    
    # Initialize a list to store errors for each lambda value
    all_errors <- list()
    
    for (i in lam_values) {
      all_errors[[paste0("lambda_", i)]] <- mclapply(1:kf, function(fold) {
        # Identify Training and Testing Sets
        training_indices_y <- unlist(lapply(fold_indices[-fold], `[[`, 'y_indices'))
        testing_indices_y <- fold_indices[[fold]]$y_indices
        training_indices_wg2 <- unlist(lapply(fold_indices[-fold], `[[`, 'wg2_indices'))
        testing_indices_wg2 <- fold_indices[[fold]]$wg2_indices
        
        y_train <- y[training_indices_y, ]
        y_test <- y[testing_indices_y, ]
        wg1_train <- wg1[training_indices_y, ]
        e1_train <- e1[training_indices_y,]
        wg2_train <- wg2[training_indices_wg2, ]
        e2_train <- e2[training_indices_wg2]
        wg2_test <- wg2[testing_indices_wg2, ]
        e2_test <- e2[testing_indices_wg2]
        
        sapply(a_values, function(a) {
          initial = c(sigma1, sigma2, sigmau, truealpha)
          result <- EM(initial, lambda = i, a_para = a, wg1 = wg1_train, y = y_train, 
                       wg2 = wg2_train, u = u, truealpha = truealpha, e2 = e2_train, m1 = m1, 
                       n1 = nrow(y_train), n2 = nrow(wg2_train), p1 = p1, k = k)
          alpha_est <- result$final_parameters[-c(1:3,length(result$final_parameters))]
          mean((z_est - z_test)^2)
        })
      }, mc.cores = 30, mc.preschedule = FALSE)
    }
    
    col_means <- lapply(all_errors, function(lambda_errors) {
      sapply(1:length(lambda_errors[[1]]), function(col) {
        mean(sapply(lambda_errors, function(fold_errors) fold_errors[col]))
      })
    })
    # Finding the optimal value of a
    optimal_a <- a_values[sapply(col_means, which.min)]
    return(list(a=optimal_a, colmean = col_means))
  }
  msecal <- function(wg2,mle){
    mse <- rep(0,length(lam_values))
    for (i in 1:length(lam_values)){
      alphaest0 <- mle[i,4:(m1+3)]
      z_est0 <- wg2%*%u%*%alphaest0  + e2
      z <- wg2%*%u%*%truealpha + e2
      mse[i] <- mean((z- z_est0)^2)
    }
    return(mse)
  }
  
  if (size=="small"){
    data_simu <- function(n1,n2,m1,p1,k,sigma1,sigma2,sigmau,truealpha,
                          size=c("small","large","hcsmall","htsmall")){
      # Function to standardize a matrix by column
      standardize <- function(mat) {
        col_means <- colMeans(mat)
        col_sds <- apply(mat, 2, sd)
        standardized_mat <- scale(mat, center = col_means, scale = col_sds)
        return(standardized_mat)
      }
      
      generatewg <- function(n1){
        sequ1 <- rbinom(n1,1,0.2) ## asssume MAF=0.5 or 0.2
        sequ2 <- rbinom(n1,1,0.2)
        allele <- colSums(rbind(sequ1,sequ2))
        return(allele)
      }
      
      wg1 <- foreach(i=1:p1, .combine=cbind, .packages="foreach") %do% generatewg(n1)
      wg2 <- foreach(i=1:p1, .combine=cbind, .packages="foreach") %do% generatewg(n2)
      
      normalize_matrix <- function(mat) {
        min_vals <- apply(mat, 2, min) # Find the minimum value of each column
        max_vals <- apply(mat, 2, max) # Find the maximum value of each column
        return((mat - min_vals) / (max_vals - min_vals))
      }
      
      wg1 <- normalize_matrix(wg1)
      wg2 <- normalize_matrix(wg2)
      # Standardize wg1 and wg2
      #wg1 <- standardize(wg1)
      #wg2 <- standardize(wg2)
      if (size=="hcsmall"){
        sigman=0.005
      }else if (size=="htsmall"){
        sigman=0.005
      }else{
        sigman=0.05
      }
      # Initialize the matrix u
      u <- matrix(rnorm(p1 * m1, 0, sqrt(sigman)), p1, m1)
      for (i in 1:m1) {
        block_start <- (i - 1) * 5 + 1
        u[block_start:(block_start + 4), i] <- rnorm(5, 0, sqrt(sigmau))
      }
      
      
      
      e1 <- matrix(rnorm(n1 * m1, 0, sqrt(sigma1)), n1, m1)
      
      e2 <- rnorm(n2,0,sqrt(sigma2))
      
      y <- wg1%*%u +e1
      z <- wg2%*%u%*%truealpha + e2
      
      
      
      M_step <- function(old, lambda, a_para,wg1, y, wg2, u, truealpha, e2, m1, n1, n2,p1, k) {
        # old<- c(sigma1,sigma2,sigmau,truealpha-0.00001)
        sigma1 <- old[1]; sigma2 <- old[2]; sigmau <- old[3]; alpha_g <- old[-c(1:3)]
        a <- a_para
        abc <-  (1 - a) * lambda * k
        ab <- a * lambda 
        wg1_product <- t(wg1) %*% wg1
        wg2_product <- t(wg2) %*% wg2
        sigma_ui_constant <- diag(1, p1, p1)
        
        results <- mclapply(1:m1, function(i) {
          sigma_ui_i <- solve(1 / sigma1 * wg1_product + (alpha_g[i]^2) / sigma2 * wg2_product + 1 / sigmau * sigma_ui_constant)
          mu_ui_i <- sigma_ui_i %*% (1 / sigma1 * t(wg1) %*% y[, i] + alpha_g[i] / sigma2 * t(wg2) %*% (wg2 %*% u[, i] %*% truealpha[i] + e2))
          
          si <- 2 * t(wg2 %*% u[, i] %*% truealpha[i] + e2) %*% wg2 %*% mu_ui_i - a * lambda * sign(alpha_g[i])
          abc <-  (1 - a) * lambda * k[i]
          ab <- a * lambda 
          alphaest <- ifelse(abs(si) <= abc | abs(si) <= ab, 0,(1 - abc / (abs(si))) * 
                               (2 * t(mu_ui_i) %*% t(wg2) %*% wg2 %*% mu_ui_i + 2 * sum(diag(wg2 %*% sigma_ui_i %*% t(wg2))))^(-1) * si)
          newz <- wg2 %*% u[, i] * truealpha[i] + e2
          E1 <- t(y[, i]) %*% y[, i] - 2 * t(y[, i]) %*% wg1 %*% mu_ui_i + t(mu_ui_i) %*% t(wg1) %*% wg1 %*% mu_ui_i + sum(diag(t(wg1) %*% wg1 %*% sigma_ui_i))
          E2 <- t(newz) %*% newz - 2 * alphaest * t(newz) %*% wg2 %*% mu_ui_i + t(as.numeric(alphaest) * wg2 %*% mu_ui_i) %*% (as.numeric(alphaest) * wg2 %*% mu_ui_i) + sum(diag(as.numeric(alphaest)^2 * wg2 %*% sigma_ui_i %*% t(wg2)))
          E3 <- t(mu_ui_i) %*% mu_ui_i + sum(diag(sigma_ui_i))
          list(E1 = E1, E2 = E2, E3 = E3, alphaest = alphaest, mu_ui_i = mu_ui_i)
        }, mc.cores = 30)
        
        sigma1est <- 1 / (m1 * n1) * sum(sapply(results, function(x) x$E1))
        sigma2est <- 1 / (m1 * n2) * sum(sapply(results, function(x) x$E2))
        sigmauest <- 1 / (m1 * p1) * sum(sapply(results, function(x) x$E3))
        log_likelihood_y = -0.5 * n1*m1 * log(2 * pi * sigma1est) - 0.5 * sum(sapply(results, function(x) x$E1)) / sigma1est;
        log_likelihood_z = -0.5 * n2 * log(2 * pi * sigma2est) - 0.5 * sum(sapply(results, function(x) x$E2)) / sigma2est;
        log_likelihood_u = -0.5 * p1*m1 * log(2 * pi * sigmauest) - 0.5 * sum(sapply(results, function(x) x$E3)) / sigmauest;
        log_likelihood = log_likelihood_y + log_likelihood_z + log_likelihood_u
        print(c(sigma1est,sigma2est,sigmauest))
        return(list(est=c(sigma1est, sigma2est, sigmauest, unlist(sapply(results, function(x) x$alphaest)))
                    ,lik = log_likelihood))
      }
      
      EM <- function(init, lambda, a_para, wg1, y, wg2, u, truealpha, e2, m1, n1, n2, 
                     p1, k,max_iter = 10, tol = 0.001) {
        # Initialize the parameters
        theta_old <- init-0.0001
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
          cat("Iteration:", iter, "parameter difference:", current_diff, "\n")
          # Check for convergence using the difference between the current and previous differences
          if (abs(current_diff - prev_diff) < tol) break
          
          # Update for the next iteration
          theta_old <- theta_new
          prev_diff <- current_diff
          iter <- iter + 1
          print(prev_diff)
        }
        
        # Return the estimated parameters
        return(list(final_parameters = c(theta_new, lambda,log_likelihood_new)))
      }
      k_fold <- function(kf,a_values,lam_values,y,wg1,e1,wg2,e2,u,m1,truealpha,p1,k){
        # Split data into K-folds
        fold_size_y <- dim(y)[1] / kf
        fold_size_wg2 <- dim(wg2)[1] / kf
        
        # Create a list of indices for each fold
        fold_indices <- lapply(1:kf, function(i) {
          list(y_indices = ((i-1)*fold_size_y + 1):(i*fold_size_y),
               wg2_indices = ((i-1)*fold_size_wg2 + 1):(i*fold_size_wg2))
        })
        
        # Initialize a list to store errors for each lambda value
        all_errors <- list()
        
        for (i in lam_values) {
          all_errors[[paste0("lambda_", i)]] <- mclapply(1:kf, function(fold) {
            # Identify Training and Testing Sets
            training_indices_y <- unlist(lapply(fold_indices[-fold], `[[`, 'y_indices'))
            testing_indices_y <- fold_indices[[fold]]$y_indices
            training_indices_wg2 <- unlist(lapply(fold_indices[-fold], `[[`, 'wg2_indices'))
            testing_indices_wg2 <- fold_indices[[fold]]$wg2_indices
            
            y_train <- y[training_indices_y, ]
            y_test <- y[testing_indices_y, ]
            wg1_train <- wg1[training_indices_y, ]
            e1_train <- e1[training_indices_y,]
            wg2_train <- wg2[training_indices_wg2, ]
            e2_train <- e2[training_indices_wg2]
            wg2_test <- wg2[testing_indices_wg2, ]
            e2_test <- e2[testing_indices_wg2]
            
            sapply(a_values, function(a) {
              initial = c(sigma1, sigma2, sigmau, truealpha)
              result <- EM(initial, lambda = i, a_para = a, wg1 = wg1_train, y = y_train, 
                           wg2 = wg2_train, u = u, truealpha = truealpha, e2 = e2_train, m1 = m1, 
                           n1 = nrow(y_train), n2 = nrow(wg2_train), p1 = p1, k = k)
              alpha_est <- result$final_parameters[-c(1:3,length(result$final_parameters))]
              mean((z_est - z_test)^2)
            })
          }, mc.cores = 30, mc.preschedule = FALSE)
        }
        
        col_means <- lapply(all_errors, function(lambda_errors) {
          sapply(1:length(lambda_errors[[1]]), function(col) {
            mean(sapply(lambda_errors, function(fold_errors) fold_errors[col]))
          })
        })
        # Finding the optimal value of a
        optimal_a <- a_values[sapply(col_means, which.min)]
        return(list(a=optimal_a, colmean = col_means))
      }
      msecal <- function(wg2,mle){
        mse <- rep(0,length(lam_values))
        for (i in 1:length(lam_values)){
          alphaest0 <- mle[i,4:(m1+3)]
          z_est0 <- wg2%*%u%*%alphaest0  + e2
          z <- wg2%*%u%*%truealpha + e2
          mse[i] <- mean((z- z_est0)^2)
        }
        return(mse)
      }
      
      if (size=="small"){
        a_mle <- k_fold(kf=5,a_values=seq(0,1,length=5), lam_values=c(seq(0,20,length=7),seq(230,500,length=6),600,750),y,wg1,e1,wg2,e2,u,m1,truealpha,p1,k)$a
        lam_values_l0 = c(0,10,15,20,seq(30,200,length=9),400,700)
        lam_values_l1 = c(0,5,10,14,seq(17,200,length=11))
        lam_values_l2 = c(0,15,30,seq(40,250,length=10),350,550)
      }else if (size=="large"){
        a_mle <- k_fold(kf=5,a_values=seq(0,1,length=5), lam_values=c(seq(0,20,length=7),seq(230,500,length=6),600,750),y,wg1,e1,wg2,e2,u,m1,truealpha,p1,k)$a
        lam_values_l0 = c(0,5,15,seq(20,300,length=10),400,500)
        lam_values_l1 = c(0,3,7,seq(15,200,length=10),300,400)
        lam_values_l2 = c(0,5,10,seq(20,300,length=10),400,500)
      }else if (size=="hcsmall"){
        a_mle <- k_fold(kf=5,a_values=seq(0,1,length=5), lam_values=lam_values_l0,y,wg1,e1,wg2,e2,u,m1,truealpha,p1,k)$a
        lam_values_l0 = c(0,1,2,3,seq(3.5,40,length=10),50)
        lam_values_l1 = c(0,1,2,3,5,seq(6,25,length=9),50)
        lam_values_l2 = c(0,3,6,10,seq(12,50,length=10),100)
      }else if (size=="htsmall"){
        a_mle <- k_fold(kf=5,a_values=seq(0,1,length=5), lam_values=lam_values_l0,y,wg1,e1,wg2,e2,u,m1,truealpha,p1,k)$a
        lam_values_l0 = c(0,3,7,10,seq(15,60,length=10),100)
        lam_values_l1 = c(seq(0,10,length=11),12,15,20,200)
        lam_values_l2 =  c(0,3,6,seq(8,33,length=10),50,100)
      }
      
      
      mle_l0 <- foreach(i = seq_along(lam_values_l0),.combine = rbind) %do% {
        lambda <- lam_values_l0[i]
        a_para <- a_mle[i]
        EM(init=c(sigma1, sigma2, sigmau, truealpha), 
           lambda = lambda, a_para = a_para, 
           wg1, y, wg2, u, truealpha, e2, m1, n1, n2, p1, k)$final_parameters
      }
      
      mle_l1 <- foreach(i=lam_values_l1,.combine = rbind) %do% 
        EM(init=c(sigma1, sigma2, sigmau, truealpha), 
           lambda = i, a_para = 0,  wg1, y, wg2, u, truealpha, e2, m1, n1, n2, p1, k)$final_parameters
      mle_l2 <-  foreach(i=lam_values_l2,.combine = rbind) %do% 
        EM(init=c(sigma1, sigma2, sigmau, truealpha), 
           lambda = i, a_para = 1,  wg1, y, wg2, u, truealpha, e2, m1, n1, n2, p1, k)$final_parameters
      
      ### CoMM
      newu <- matrix(nrow = 5, ncol = m1)
      mle_comm <- matrix(NA,m1,4)
      p_values_lrt <- numeric(length = m1)
      M_step_comm <- function(old, w1, y, w2, z, n1, n2) {
        sigma1 <- old[1]; sigma2 <- old[2]; sigmau <- old[3]; alpha_g <- old[-c(1:3)]
        
        sigma_ui_i <- solve( 1 / sigma1 * t(w1) %*% w1 + (alpha_g^2) / sigma2 * t(w2) %*% w2 + 1 / sigmau * diag(1, 5, 5))
        mu_ui_i <- sigma_ui_i %*% (1 / sigma1 * t(w1) %*% y+ alpha_g / sigma2 * t(w2) %*% z)
        alphaest <- solve(t(mu_ui_i)%*%t(w2)%*%w2%*%mu_ui_i +sum(diag(w2 %*% sigma_ui_i %*% t(w2))))*t(z)%*%w2 %*% mu_ui_i
        E1 <- t(y) %*% y- 2 * t(y) %*% w1 %*% mu_ui_i + t(mu_ui_i) %*% t(w1) %*% w1 %*% mu_ui_i
        E2 <- t(z) %*% z - 2 * alphaest * t(z) %*% w2 %*% mu_ui_i + alphaest^2 * t(mu_ui_i)%*% t(w2)%*%w2 %*%mu_ui_i
        E3 <- t(mu_ui_i) %*% mu_ui_i
        sigma1est <- 1/n1 * (E1 + sum(diag(t(w1) %*% w1 %*% sigma_ui_i)))
        sigma2est <- 1/n2*(E2  + alphaest^2*sum(diag(t(w2) %*% w2 %*% sigma_ui_i)))
        sigmauest <- 1/5*(E3 + sum(diag(sigma_ui_i)))
        
        log_likelihood_y <- -0.5 * n1 * log(2 * pi * sigma1) - 0.5 * E1 / sigma1
        log_likelihood_z <- -0.5 * n2 * log(2 * pi * sigma2) - 0.5 * E2 / sigma2
        log_likelihood_u <- -0.5 * 5 * log(2 * pi * sigmau) - 0.5 * E3/ sigmau
        
        total_log_likelihood <- log_likelihood_y + log_likelihood_z + log_likelihood_u
        return(list(esti = c(sigma1est,sigma2est, sigmauest,alphaest), 
                    log_likelihood = total_log_likelihood))
      }
      EM_comm <- function(init, w1, y, w2, z, n1, n2, max_iter = 20, tol = 0.001) {
        theta_old <- init - 0.0001
        iter <- 1
        prev_diff <- Inf
        
        for (i in 1:max_iter) {
          step_result <- M_step_comm(old = theta_old, w1 = w1, y = y, w2 = w2, z = z, n1 = n1, n2 = n2)
          theta_new <- step_result$esti
          
          current_diff <- mean(abs(theta_new - theta_old))
          if (abs(current_diff - prev_diff) < tol) break
          
          theta_old <- theta_new
          prev_diff <- current_diff
          iter <- iter + 1
          print(prev_diff)
        }
        return(list(final_parameters = theta_new))
      }
      
      loglik <- function(para, w1, y, w2, z, n1, n2){
        sigma1 <- para[1]; sigma2 <- para[2]; sigmau <- para[3]; alpha_g <- para[-c(1:3)]
        
        sigma_ui_i <- solve( 1 / sigma1 * t(w1) %*% w1 + (alpha_g^2) / sigma2 * t(w2) %*% w2 + 1 / sigmau * diag(1, 5, 5))
        mu_ui_i <- sigma_ui_i %*% (1 / sigma1 * t(w1) %*% y+ alpha_g / sigma2 * t(w2) %*% z)
        E1 <- t(y) %*% y- 2 * t(y) %*% w1 %*% mu_ui_i + t(mu_ui_i) %*% t(w1) %*% w1 %*% mu_ui_i
        E2 <- t(z) %*% z - 2 * alpha_g * t(z) %*% w2 %*% mu_ui_i + alpha_g^2 * t(mu_ui_i)%*% t(w2)%*%w2 %*%mu_ui_i
        E3 <- t(mu_ui_i) %*% mu_ui_i
        log_likelihood_y <- -0.5 * n1 * log(2 * pi * sigma1) - 0.5 * E1 / sigma1
        log_likelihood_z <- -0.5 * n2 * log(2 * pi * sigma2) - 0.5 * E2 / sigma2
        log_likelihood_u <- -0.5 * 5 * log(2 * pi * sigmau) - 0.5 * E3/ sigmau
        total_log_likelihood <- log_likelihood_y + log_likelihood_z + log_likelihood_u
        return(total_log_likelihood)
      }
      # Loop over each column of the original matrix
      for (i in 1:m1) {
        # Calculate the row indices for the current subset
        row_indices <- ((i - 1) * 5 + 1):(i * 5)
        # Extract the 5x1 subset and assign it to the corresponding column in newu
        newu[, i] <- u[row_indices, i]
        w1 <- wg1[,((i-1)*5 + 1):(5*i)]
        w2 <- wg2[,((i-1)*5 + 1):(5*i)]
        ynew <- w1 %*% as.matrix(newu[, i]) + e1[,i]
        znew <- truealpha[i] * w2 %*% as.matrix(newu[, i]) + e2
        em1 <- EM_comm(init = c(sigma1, sigma2, sigmau, truealpha[i]),
                       w1 = w1, y = ynew, w2 = w2, z = znew, 
                       n1 = n1, n2 = n2)
        mle_comm[i,] <- em1$final_parameters
        # Alternative hypothesis: Use the estimated alphaest[i]
        log_likelihood_alt <- loglik(c(mle_comm[i,1], mle_comm[i,2], mle_comm[i,3], mle_comm[i,4]), w1 = w1, y = ynew, w2 = w2, z = znew, n1 = n1, n2 = n2)
        
        
        # Null hypothesis: alpha_i is set to 0
        log_likelihood_null <- loglik(c(mle_comm[i,1], mle_comm[i,2], mle_comm[i,3], 0), w1 = w1, y = ynew, w2 = w2, z = znew, n1 = n1, n2 = n2)
        
        # Calculate the likelihood ratio and the p-value
        lr_statistic <- 2 * (log_likelihood_alt - log_likelihood_null)
        p_values_lrt[i] <- pchisq(lr_statistic, df = 1, lower.tail = FALSE)
      }
      
      mle_comm <- c(colMeans(mle_comm[,1:3]),mle_comm[,4]) # length 183
      
      
      
      
      
      
      
      ### Predixcan EN
      library(glmnet)
      library(foreach) # for parallel processing
      
      coefficients_en <- list()
      association_results <- vector("list", length = m1)
      p_values_EN <- rep(NA, m1) # Pre-allocate space for p-values
      
      for (i in 1:m1) {
        # Calculate the row indices for the current subset
        row_indices <- ((i - 1) * 5 + 1):(i * 5)
        
        # Extract the 5x1 subset for the gene
        u_i <- u[row_indices, i]
        w1_i <- wg1[, row_indices] # Genotype data for training
        w2_i <- wg2[, row_indices] # Genotype data for testing/association
        y_i <- w1_i %*% u_i + e1[, i] # Gene expression
        z_i <- truealpha[i] * w2_i %*% u_i + e2 # Phenotype
        
        # Fit the Elastic Net model for the current gene using the training data
        cv_fit <- cv.glmnet(w1_i, y_i, alpha=0, family="gaussian")
        best_lambda <- cv_fit$lambda.min
        fit <- glmnet(w1_i, y_i, alpha=0, lambda=best_lambda)
        
        # Store the coefficients from the Elastic Net model
        coefficients_en[[i]] <- coef(fit, s = "lambda.min")
        
        # Calculate the predicted gene expression levels for the samples in the association study
        predicted_expression <- predict(fit, newx=w2_i)
        
        # Perform the association test between the predicted expression and the phenotype
        association_results[[i]] <- lm(z_i ~ predicted_expression)
        
        # Check if there are predictors before attempting to access their p-value
        if (dim(summary(association_results[[i]])$coefficients)[1] > 1) {
          p_values_EN[i] <- summary(association_results[[i]])$coefficients[2, "Pr(>|t|)"]
        } else {
          # Handle the case when there are no predictors
          p_values_EN[i] <- NA  # Or some other appropriate action
        }
        
      }
      
      
      
      ###Predixcan Lasso
      coefficients_lasso <- list()
      association_results <- vector("list", length = m1)
      p_values_lasso <- rep(NA, m1) # Pre-allocate space for p-values
      
      for (i in 1:m1) {
        # Calculate the row indices for the current subset
        row_indices <- ((i - 1) * 5 + 1):(i * 5)
        
        # Extract the 5x1 subset for the gene
        u_i <- u[row_indices, i]
        w1_i <- wg1[, row_indices] # Genotype data for training
        w2_i <- wg2[, row_indices] # Genotype data for testing/association
        y_i <- w1_i %*% u_i + e1[, i] # Gene expression
        z_i <- truealpha[i] * w2_i %*% u_i + e2 # Phenotype
        
        # Fit the LASSO model for the current gene using the training data
        cv_fit <- cv.glmnet(w1_i, y_i, alpha=0.9, family="gaussian")
        best_lambda <- cv_fit$lambda.min
        fit <- glmnet(w1_i, y_i, alpha=0.9, lambda=best_lambda)
        
        # Store the coefficients from the LASSO model
        coefficients_lasso[[i]] <- coef(fit, s = "lambda.min")
        
        # Calculate the predicted gene expression levels for the samples in the association study
        predicted_expression <- predict(fit, newx=w2_i)
        
        # Perform the association test between the predicted expression and the phenotype
        association_results[[i]] <- lm(z_i ~ predicted_expression)
        
        
        # Check if there are predictors before attempting to access their p-value
        if (dim(summary(association_results[[i]])$coefficients)[1] > 1) {
          p_values_lasso[i] <- summary(association_results[[i]])$coefficients[2, "Pr(>|t|)"]
        } else {
          # Handle the case when there are no predictors
          p_values_lasso[i] <- NA  # Or some other appropriate action
        }
      }
      
      ### multivariate TWAS
      # Assuming m1 is the number of genes, which is 180 in this case
      u_i_est <- matrix(0, ncol = p1, nrow = m1) # Corrected dimensions to 180 x 900
      # Estimating gene expression effects from wg1
      library(MASS)
      for (i in 1:m1) {
        y_i <- y[, i]
        # OLS estimation for each gene's effect
        u_i_est[i, ] <- ginv(t(wg1) %*% wg1)%*% t(wg1) %*% y_i
      }
      predictor <- wg2 %*% t(u_i_est) # Corrected multiplication
      model <- lm(z ~ predictor + 0) 
      # Extracting coefficients (assuming one coefficient per predictor)
      mle_multi <- coef(model)
      p_values_multi <- summary(model)$coefficients[, 4] # 4th column contains p-values
      residuals <- model$residuals
      mse_multi <- mean(residuals^2)
      
      
      result_list <- list(
        mle_l0 = mle_l0,
        mle_l1 = mle_l1,
        mle_l2 = mle_l2,mle_comm=mle_comm, mle_multi = mle_multi,
        p_values_comm =p_values_lrt, p_values_EN=p_values_EN,p_values_lasso=p_values_lasso,
        p_values_multi = p_values_multi,
        mse = rbind(msecal(wg2,mle_l0),msecal(wg2,mle_l1),msecal(wg2,mle_l2)),mse_multi=mse_multi
      )
      return(result_list)
    }
  }else if (size=="large"){
    #a_mle <- c(0.5,15)
    a_mle <-  c(0,0.5,rep(1,13))
    #lam_values_l0 = c(0,5,10,25,50,seq(65,250,length=8),300,400) 
    #lam_values_l1 = c(seq(0,21,length=5),seq(25,90,length=7),seq(105,150,length=3))
    #lam_values_l2 = c(seq(0,90,length=5),seq(110,300,length=9),350) # lar
    lam_values_l0 = c(0,5,15,seq(20,300,length=10),400,500)
    lam_values_l1 = c(0,3,7,seq(15,200,length=10),300,400)
    lam_values_l2 = c(0,5,10,seq(20,300,length=10),400,500)
  }else if (size=="hcsmall"){
    # a_mle <-  c(0, 0.5, rep(1,13))
    # a_mle <- k_fold(kf=5,a_values=seq(0,1,length=5), lam_values=lam_values_l0,y,wg1,e1,wg2,e2,u,m1,truealpha,p1,k)$a
    a_mle <- c(0.00, 0.25, 0.75,rep(1,12)) 
    lam_values_l0 = c(0,1,2,3,seq(3.5,40,length=10),50)
    lam_values_l1 = c(0,1,2,3,5,seq(6,25,length=9),50)
    lam_values_l2 = c(0,3,6,10,seq(12,50,length=10),100)
  }else if (size=="htsmall"){
    #a_mle <-  c(0.00, 0.5, rep(1,13))
    # a_mle <- k_fold(kf=5,a_values=seq(0,1,length=5), lam_values=lam_values_l0,y,wg1,e1,wg2,e2,u,m1,truealpha,p1,k)$a
    a_mle <- c(0.00, 0.50, 0.75,rep(1,12))
    lam_values_l0 = c(0,3,7,10,seq(15,60,length=10),100)
    lam_values_l1 = c(seq(0,10,length=11),12,15,20,200)
    lam_values_l2 =  c(0,3,6,seq(8,33,length=10),50,100)
  }
  
  
  mle_l0 <- foreach(i = seq_along(lam_values_l0),.combine = rbind) %do% {
    lambda <- lam_values_l0[i]
    a_para <- a_mle[i]
    EM(init=c(sigma1, sigma2, sigmau, truealpha), 
       lambda = lambda, a_para = a_para, 
       wg1, y, wg2, u, truealpha, e2, m1, n1, n2, p1, k)$final_parameters
  }
  
  mle_l1 <- foreach(i=lam_values_l1,.combine = rbind) %do% 
    EM(init=c(sigma1, sigma2, sigmau, truealpha), 
       lambda = i, a_para = 0,  wg1, y, wg2, u, truealpha, e2, m1, n1, n2, p1, k)$final_parameters
  mle_l2 <-  foreach(i=lam_values_l2,.combine = rbind) %do% 
    EM(init=c(sigma1, sigma2, sigmau, truealpha), 
       lambda = i, a_para = 1,  wg1, y, wg2, u, truealpha, e2, m1, n1, n2, p1, k)$final_parameters
  
  ### CoMM
  newu <- matrix(nrow = 5, ncol = m1)
  mle_comm <- matrix(NA,m1,4)
  p_values_lrt <- numeric(length = m1)
  M_step_comm <- function(old, w1, y, w2, z, n1, n2) {
    sigma1 <- old[1]; sigma2 <- old[2]; sigmau <- old[3]; alpha_g <- old[-c(1:3)]
    
    sigma_ui_i <- solve( 1 / sigma1 * t(w1) %*% w1 + (alpha_g^2) / sigma2 * t(w2) %*% w2 + 1 / sigmau * diag(1, 5, 5))
    mu_ui_i <- sigma_ui_i %*% (1 / sigma1 * t(w1) %*% y+ alpha_g / sigma2 * t(w2) %*% z)
    alphaest <- solve(t(mu_ui_i)%*%t(w2)%*%w2%*%mu_ui_i +sum(diag(w2 %*% sigma_ui_i %*% t(w2))))*t(z)%*%w2 %*% mu_ui_i
    E1 <- t(y) %*% y- 2 * t(y) %*% w1 %*% mu_ui_i + t(mu_ui_i) %*% t(w1) %*% w1 %*% mu_ui_i
    E2 <- t(z) %*% z - 2 * alphaest * t(z) %*% w2 %*% mu_ui_i + alphaest^2 * t(mu_ui_i)%*% t(w2)%*%w2 %*%mu_ui_i
    E3 <- t(mu_ui_i) %*% mu_ui_i
    sigma1est <- 1/n1 * (E1 + sum(diag(t(w1) %*% w1 %*% sigma_ui_i)))
    sigma2est <- 1/n2*(E2  + alphaest^2*sum(diag(t(w2) %*% w2 %*% sigma_ui_i)))
    sigmauest <- 1/5*(E3 + sum(diag(sigma_ui_i)))
    
    log_likelihood_y <- -0.5 * n1 * log(2 * pi * sigma1) - 0.5 * E1 / sigma1
    log_likelihood_z <- -0.5 * n2 * log(2 * pi * sigma2) - 0.5 * E2 / sigma2
    log_likelihood_u <- -0.5 * 5 * log(2 * pi * sigmau) - 0.5 * E3/ sigmau
    
    total_log_likelihood <- log_likelihood_y + log_likelihood_z + log_likelihood_u
    return(list(esti = c(sigma1est,sigma2est, sigmauest,alphaest), 
                log_likelihood = total_log_likelihood))
  }
  EM_comm <- function(init, w1, y, w2, z, n1, n2, max_iter = 20, tol = 0.001) {
    theta_old <- init - 0.0001
    iter <- 1
    prev_diff <- Inf
    
    for (i in 1:max_iter) {
      step_result <- M_step_comm(old = theta_old, w1 = w1, y = y, w2 = w2, z = z, n1 = n1, n2 = n2)
      theta_new <- step_result$esti
      
      current_diff <- mean(abs(theta_new - theta_old))
      if (abs(current_diff - prev_diff) < tol) break
      
      theta_old <- theta_new
      prev_diff <- current_diff
      iter <- iter + 1
      print(prev_diff)
    }
    return(list(final_parameters = theta_new))
  }
  
  loglik <- function(para, w1, y, w2, z, n1, n2){
    sigma1 <- para[1]; sigma2 <- para[2]; sigmau <- para[3]; alpha_g <- para[-c(1:3)]
    
    sigma_ui_i <- solve( 1 / sigma1 * t(w1) %*% w1 + (alpha_g^2) / sigma2 * t(w2) %*% w2 + 1 / sigmau * diag(1, 5, 5))
    mu_ui_i <- sigma_ui_i %*% (1 / sigma1 * t(w1) %*% y+ alpha_g / sigma2 * t(w2) %*% z)
    E1 <- t(y) %*% y- 2 * t(y) %*% w1 %*% mu_ui_i + t(mu_ui_i) %*% t(w1) %*% w1 %*% mu_ui_i
    E2 <- t(z) %*% z - 2 * alpha_g * t(z) %*% w2 %*% mu_ui_i + alpha_g^2 * t(mu_ui_i)%*% t(w2)%*%w2 %*%mu_ui_i
    E3 <- t(mu_ui_i) %*% mu_ui_i
    log_likelihood_y <- -0.5 * n1 * log(2 * pi * sigma1) - 0.5 * E1 / sigma1
    log_likelihood_z <- -0.5 * n2 * log(2 * pi * sigma2) - 0.5 * E2 / sigma2
    log_likelihood_u <- -0.5 * 5 * log(2 * pi * sigmau) - 0.5 * E3/ sigmau
    total_log_likelihood <- log_likelihood_y + log_likelihood_z + log_likelihood_u
    return(total_log_likelihood)
  }
  # Loop over each column of the original matrix
  for (i in 1:m1) {
    # Calculate the row indices for the current subset
    row_indices <- ((i - 1) * 5 + 1):(i * 5)
    # Extract the 5x1 subset and assign it to the corresponding column in newu
    newu[, i] <- u[row_indices, i]
    w1 <- wg1[,((i-1)*5 + 1):(5*i)]
    w2 <- wg2[,((i-1)*5 + 1):(5*i)]
    ynew <- w1 %*% as.matrix(newu[, i]) + e1[,i]
    znew <- truealpha[i] * w2 %*% as.matrix(newu[, i]) + e2
    em1 <- EM_comm(init = c(sigma1, sigma2, sigmau, truealpha[i]),
                   w1 = w1, y = ynew, w2 = w2, z = znew, 
                   n1 = n1, n2 = n2)
    mle_comm[i,] <- em1$final_parameters
    # Alternative hypothesis: Use the estimated alphaest[i]
    log_likelihood_alt <- loglik(c(mle_comm[i,1], mle_comm[i,2], mle_comm[i,3], mle_comm[i,4]), w1 = w1, y = ynew, w2 = w2, z = znew, n1 = n1, n2 = n2)
    
    
    # Null hypothesis: alpha_i is set to 0
    log_likelihood_null <- loglik(c(mle_comm[i,1], mle_comm[i,2], mle_comm[i,3], 0), w1 = w1, y = ynew, w2 = w2, z = znew, n1 = n1, n2 = n2)
    
    # Calculate the likelihood ratio and the p-value
    lr_statistic <- 2 * (log_likelihood_alt - log_likelihood_null)
    p_values_lrt[i] <- pchisq(lr_statistic, df = 1, lower.tail = FALSE)
  }
  
  mle_comm <- c(colMeans(mle_comm[,1:3]),mle_comm[,4]) # length 183
  
  
  
  
  
  
  
  ### Predixcan EN
  library(glmnet)
  library(foreach) # for parallel processing
  
  coefficients_en <- list()
  association_results <- vector("list", length = m1)
  p_values_EN <- rep(NA, m1) # Pre-allocate space for p-values
  
  for (i in 1:m1) {
    # Calculate the row indices for the current subset
    row_indices <- ((i - 1) * 5 + 1):(i * 5)
    
    # Extract the 5x1 subset for the gene
    u_i <- u[row_indices, i]
    w1_i <- wg1[, row_indices] # Genotype data for training
    w2_i <- wg2[, row_indices] # Genotype data for testing/association
    y_i <- w1_i %*% u_i + e1[, i] # Gene expression
    z_i <- truealpha[i] * w2_i %*% u_i + e2 # Phenotype
    
    # Fit the Elastic Net model for the current gene using the training data
    cv_fit <- cv.glmnet(w1_i, y_i, alpha=0, family="gaussian")
    best_lambda <- cv_fit$lambda.min
    fit <- glmnet(w1_i, y_i, alpha=0, lambda=best_lambda)
    
    # Store the coefficients from the Elastic Net model
    coefficients_en[[i]] <- coef(fit, s = "lambda.min")
    
    # Calculate the predicted gene expression levels for the samples in the association study
    predicted_expression <- predict(fit, newx=w2_i)
    
    # Perform the association test between the predicted expression and the phenotype
    association_results[[i]] <- lm(z_i ~ predicted_expression)
    
    # Check if there are predictors before attempting to access their p-value
    if (dim(summary(association_results[[i]])$coefficients)[1] > 1) {
      p_values_EN[i] <- summary(association_results[[i]])$coefficients[2, "Pr(>|t|)"]
    } else {
      # Handle the case when there are no predictors
      p_values_EN[i] <- NA  # Or some other appropriate action
    }
    
  }
  
  
  
  ###Predixcan Lasso
  coefficients_lasso <- list()
  association_results <- vector("list", length = m1)
  p_values_lasso <- rep(NA, m1) # Pre-allocate space for p-values
  
  for (i in 1:m1) {
    # Calculate the row indices for the current subset
    row_indices <- ((i - 1) * 5 + 1):(i * 5)
    
    # Extract the 5x1 subset for the gene
    u_i <- u[row_indices, i]
    w1_i <- wg1[, row_indices] # Genotype data for training
    w2_i <- wg2[, row_indices] # Genotype data for testing/association
    y_i <- w1_i %*% u_i + e1[, i] # Gene expression
    z_i <- truealpha[i] * w2_i %*% u_i + e2 # Phenotype
    
    # Fit the LASSO model for the current gene using the training data
    cv_fit <- cv.glmnet(w1_i, y_i, alpha=0.9, family="gaussian")
    best_lambda <- cv_fit$lambda.min
    fit <- glmnet(w1_i, y_i, alpha=0.9, lambda=best_lambda)
    
    # Store the coefficients from the LASSO model
    coefficients_lasso[[i]] <- coef(fit, s = "lambda.min")
    
    # Calculate the predicted gene expression levels for the samples in the association study
    predicted_expression <- predict(fit, newx=w2_i)
    
    # Perform the association test between the predicted expression and the phenotype
    association_results[[i]] <- lm(z_i ~ predicted_expression)
    
    
    # Check if there are predictors before attempting to access their p-value
    if (dim(summary(association_results[[i]])$coefficients)[1] > 1) {
      p_values_lasso[i] <- summary(association_results[[i]])$coefficients[2, "Pr(>|t|)"]
    } else {
      # Handle the case when there are no predictors
      p_values_lasso[i] <- NA  # Or some other appropriate action
    }
  }
  
  ### multivariate TWAS
  # Assuming m1 is the number of genes, which is 180 in this case
  u_i_est <- matrix(0, ncol = p1, nrow = m1) # Corrected dimensions to 180 x 900
  # Estimating gene expression effects from wg1
  library(MASS)
  for (i in 1:m1) {
    y_i <- y[, i]
    # OLS estimation for each gene's effect
    u_i_est[i, ] <- ginv(t(wg1) %*% wg1)%*% t(wg1) %*% y_i
  }
  predictor <- wg2 %*% t(u_i_est) # Corrected multiplication
  model <- lm(z ~ predictor + 0) 
  # Extracting coefficients (assuming one coefficient per predictor)
  mle_multi <- coef(model)
  p_values_multi <- summary(model)$coefficients[, 4] # 4th column contains p-values
  residuals <- model$residuals
  mse_multi <- mean(residuals^2)
  
  
  result_list <- list(
    mle_l0 = mle_l0,
    mle_l1 = mle_l1,
    mle_l2 = mle_l2,mle_comm=mle_comm, mle_multi = mle_multi,
    p_values_comm =p_values_lrt, p_values_EN=p_values_EN,p_values_lasso=p_values_lasso,
    p_values_multi = p_values_multi,
    mse = rbind(msecal(wg2,mle_l0),msecal(wg2,mle_l1),msecal(wg2,mle_l2)),mse_multi=mse_multi
  )
  return(result_list)
}
