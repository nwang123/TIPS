\name{data_simu}
\alias{data_simu}

\title{
TIPS
}
\description{
This is the package code for TIPS: a novel pathway-guided transcriptome-wide association studies method to reveal biological processes underlying complex polygenic traits.
}
\usage{
data_simu(n1, n2, m1, p1, k, sigma1, sigma2, sigmau, truealpha, size = c("small", "large", "hcsmall", "htsmall"))
}

\arguments{
  \item{n1}{
sample size for gene expression data matrix
}
  \item{n2}{
sample size for phenotype data
}
  \item{m1}{
number of genes
}
  \item{p1}{
number of snps
}
  \item{k}{
group structure of pathway
}
  \item{sigma1}{
Variance for e1
}
  \item{sigma2}{
Variance for e2
}
  \item{sigmau}{
variance for u
}
  \item{truealpha}{
the effects of gene expression on the phenotype
}
  \item{size}{
different scenarios 
}
}

\author{
Neng Wang
}

\examples{
## The function is currently defined as
function (n1, n2, m1, p1, k, sigma1, sigma2, sigmau, truealpha, 
    size = c("small", "large", "hcsmall", "htsmall")) 
{
    standardize <- function(mat) {
        col_means <- colMeans(mat)
        col_sds <- apply(mat, 2, sd)
        standardized_mat <- scale(mat, center = col_means, scale = col_sds)
        return(standardized_mat)
    }
    generatewg <- function(n1) {
        sequ1 <- rbinom(n1, 1, 0.2)
        sequ2 <- rbinom(n1, 1, 0.2)
        allele <- colSums(rbind(sequ1, sequ2))
        return(allele)
    }
    wg1 <- foreach(i = 1:p1, .combine = cbind, .packages = "foreach") \%do\% 
        generatewg(n1)
    wg2 <- foreach(i = 1:p1, .combine = cbind, .packages = "foreach") \%do\% 
        generatewg(n2)
    normalize_matrix <- function(mat) {
        min_vals <- apply(mat, 2, min)
        max_vals <- apply(mat, 2, max)
        return((mat - min_vals)/(max_vals - min_vals))
    }
    wg1 <- normalize_matrix(wg1)
    wg2 <- normalize_matrix(wg2)
    if (size == "hcsmall") {
        sigman = 0.005
    }
    else if (size == "htsmall") {
        sigman = 0.005
    }
    else {
        sigman = 0.05
    }
    u <- matrix(rnorm(p1 * m1, 0, sqrt(sigman)), p1, m1)
    for (i in 1:m1) {
        block_start <- (i - 1) * 5 + 1
        u[block_start:(block_start + 4), i] <- rnorm(5, 0, sqrt(sigmau))
    }
    e1 <- matrix(rnorm(n1 * m1, 0, sqrt(sigma1)), n1, m1)
    e2 <- rnorm(n2, 0, sqrt(sigma2))
    y <- wg1 \%*\% u + e1
    z <- wg2 \%*\% u \%*\% truealpha + e2
    M_step <- function(old, lambda, a_para, wg1, y, wg2, u, truealpha, 
        e2, m1, n1, n2, p1, k) {
        sigma1 <- old[1]
        sigma2 <- old[2]
        sigmau <- old[3]
        alpha_g <- old[-c(1:3)]
        a <- a_para
        abc <- (1 - a) * lambda * k
        ab <- a * lambda
        wg1_product <- t(wg1) \%*\% wg1
        wg2_product <- t(wg2) \%*\% wg2
        sigma_ui_constant <- diag(1, p1, p1)
        results <- mclapply(1:m1, function(i) {
            sigma_ui_i <- solve(1/sigma1 * wg1_product + (alpha_g[i]^2)/sigma2 * 
                wg2_product + 1/sigmau * sigma_ui_constant)
            mu_ui_i <- sigma_ui_i \%*\% (1/sigma1 * t(wg1) \%*\% 
                y[, i] + alpha_g[i]/sigma2 * t(wg2) \%*\% (wg2 \%*\% 
                u[, i] \%*\% truealpha[i] + e2))
            si <- 2 * t(wg2 \%*\% u[, i] \%*\% truealpha[i] + e2) \%*\% 
                wg2 \%*\% mu_ui_i - a * lambda * sign(alpha_g[i])
            abc <- (1 - a) * lambda * k[i]
            ab <- a * lambda
            alphaest <- ifelse(abs(si) <= abc | abs(si) <= ab, 
                0, (1 - abc/(abs(si))) * (2 * t(mu_ui_i) \%*\% 
                  t(wg2) \%*\% wg2 \%*\% mu_ui_i + 2 * sum(diag(wg2 \%*\% 
                  sigma_ui_i \%*\% t(wg2))))^(-1) * si)
            newz <- wg2 \%*\% u[, i] * truealpha[i] + e2
            E1 <- t(y[, i]) \%*\% y[, i] - 2 * t(y[, i]) \%*\% wg1 \%*\% 
                mu_ui_i + t(mu_ui_i) \%*\% t(wg1) \%*\% wg1 \%*\% mu_ui_i + 
                sum(diag(t(wg1) \%*\% wg1 \%*\% sigma_ui_i))
            E2 <- t(newz) \%*\% newz - 2 * alphaest * t(newz) \%*\% 
                wg2 \%*\% mu_ui_i + t(as.numeric(alphaest) * wg2 \%*\% 
                mu_ui_i) \%*\% (as.numeric(alphaest) * wg2 \%*\% 
                mu_ui_i) + sum(diag(as.numeric(alphaest)^2 * 
                wg2 \%*\% sigma_ui_i \%*\% t(wg2)))
            E3 <- t(mu_ui_i) \%*\% mu_ui_i + sum(diag(sigma_ui_i))
            list(E1 = E1, E2 = E2, E3 = E3, alphaest = alphaest, 
                mu_ui_i = mu_ui_i)
        }, mc.cores = 30)
        sigma1est <- 1/(m1 * n1) * sum(sapply(results, function(x) x$E1))
        sigma2est <- 1/(m1 * n2) * sum(sapply(results, function(x) x$E2))
        sigmauest <- 1/(m1 * p1) * sum(sapply(results, function(x) x$E3))
        log_likelihood_y = -0.5 * n1 * m1 * log(2 * pi * sigma1est) - 
            0.5 * sum(sapply(results, function(x) x$E1))/sigma1est
        log_likelihood_z = -0.5 * n2 * log(2 * pi * sigma2est) - 
            0.5 * sum(sapply(results, function(x) x$E2))/sigma2est
        log_likelihood_u = -0.5 * p1 * m1 * log(2 * pi * sigmauest) - 
            0.5 * sum(sapply(results, function(x) x$E3))/sigmauest
        log_likelihood = log_likelihood_y + log_likelihood_z + 
            log_likelihood_u
        print(c(sigma1est, sigma2est, sigmauest))
        return(list(est = c(sigma1est, sigma2est, sigmauest, 
            unlist(sapply(results, function(x) x$alphaest))), 
            lik = log_likelihood))
    }
    EM <- function(init, lambda, a_para, wg1, y, wg2, u, truealpha, 
        e2, m1, n1, n2, p1, k, max_iter = 10, tol = 0.001) {
        theta_old <- init - 1e-04
        iter <- 1
        prev_diff <- Inf
        for (i in 1:max_iter) {
            step_results <- M_step(old = theta_old, lambda = lambda, 
                a_para = a_para, wg1 = wg1, y = y, wg2 = wg2, 
                u = u, truealpha = truealpha, e2 = e2, m1 = m1, 
                n1 = n1, n2 = n2, p1 = p1, k = k)
            theta_new <- step_results$est
            log_likelihood_new <- step_results$lik
            current_diff <- mean(abs(theta_new - theta_old))
            cat("Iteration:", iter, "parameter difference:", 
                current_diff, "\n")
            if (abs(current_diff - prev_diff) < tol) 
                break
            theta_old <- theta_new
            prev_diff <- current_diff
            iter <- iter + 1
            print(prev_diff)
        }
        return(list(final_parameters = c(theta_new, lambda, log_likelihood_new)))
    }
    k_fold <- function(kf, a_values, lam_values, y, wg1, e1, 
        wg2, e2, u, m1, truealpha, p1, k) {
        fold_size_y <- dim(y)[1]/kf
        fold_size_wg2 <- dim(wg2)[1]/kf
        fold_indices <- lapply(1:kf, function(i) {
            list(y_indices = ((i - 1) * fold_size_y + 1):(i * 
                fold_size_y), wg2_indices = ((i - 1) * fold_size_wg2 + 
                1):(i * fold_size_wg2))
        })
        all_errors <- list()
        for (i in lam_values) {
            all_errors[[paste0("lambda_", i)]] <- mclapply(1:kf, 
                function(fold) {
                  training_indices_y <- unlist(lapply(fold_indices[-fold], 
                    `[[`, "y_indices"))
                  testing_indices_y <- fold_indices[[fold]]$y_indices
                  training_indices_wg2 <- unlist(lapply(fold_indices[-fold], 
                    `[[`, "wg2_indices"))
                  testing_indices_wg2 <- fold_indices[[fold]]$wg2_indices
                  y_train <- y[training_indices_y, ]
                  y_test <- y[testing_indices_y, ]
                  wg1_train <- wg1[training_indices_y, ]
                  e1_train <- e1[training_indices_y, ]
                  wg2_train <- wg2[training_indices_wg2, ]
                  e2_train <- e2[training_indices_wg2]
                  wg2_test <- wg2[testing_indices_wg2, ]
                  e2_test <- e2[testing_indices_wg2]
                  sapply(a_values, function(a) {
                    initial = c(sigma1, sigma2, sigmau, truealpha)
                    result <- EM(initial, lambda = i, a_para = a, 
                      wg1 = wg1_train, y = y_train, wg2 = wg2_train, 
                      u = u, truealpha = truealpha, e2 = e2_train, 
                      m1 = m1, n1 = nrow(y_train), n2 = nrow(wg2_train), 
                      p1 = p1, k = k)
                    alpha_est <- result$final_parameters[-c(1:3, 
                      length(result$final_parameters))]
                    mean((z_est - z_test)^2)
                  })
                }, mc.cores = 30, mc.preschedule = FALSE)
        }
        col_means <- lapply(all_errors, function(lambda_errors) {
            sapply(1:length(lambda_errors[[1]]), function(col) {
                mean(sapply(lambda_errors, function(fold_errors) fold_errors[col]))
            })
        })
        optimal_a <- a_values[sapply(col_means, which.min)]
        return(list(a = optimal_a, colmean = col_means))
    }
    msecal <- function(wg2, mle) {
        mse <- rep(0, length(lam_values))
        for (i in 1:length(lam_values)) {
            alphaest0 <- mle[i, 4:(m1 + 3)]
            z_est0 <- wg2 \%*\% u \%*\% alphaest0 + e2
            z <- wg2 \%*\% u \%*\% truealpha + e2
            mse[i] <- mean((z - z_est0)^2)
        }
        return(mse)
    }
    if (size == "small") {
        data_simu <- function(n1, n2, m1, p1, k, sigma1, sigma2, 
            sigmau, truealpha, size = c("small", "large", "hcsmall", 
                "htsmall")) {
            standardize <- function(mat) {
                col_means <- colMeans(mat)
                col_sds <- apply(mat, 2, sd)
                standardized_mat <- scale(mat, center = col_means, 
                  scale = col_sds)
                return(standardized_mat)
            }
            generatewg <- function(n1) {
                sequ1 <- rbinom(n1, 1, 0.2)
                sequ2 <- rbinom(n1, 1, 0.2)
                allele <- colSums(rbind(sequ1, sequ2))
                return(allele)
            }
            wg1 <- foreach(i = 1:p1, .combine = cbind, .packages = "foreach") \%do\% 
                generatewg(n1)
            wg2 <- foreach(i = 1:p1, .combine = cbind, .packages = "foreach") \%do\% 
                generatewg(n2)
            normalize_matrix <- function(mat) {
                min_vals <- apply(mat, 2, min)
                max_vals <- apply(mat, 2, max)
                return((mat - min_vals)/(max_vals - min_vals))
            }
            wg1 <- normalize_matrix(wg1)
            wg2 <- normalize_matrix(wg2)
            if (size == "hcsmall") {
                sigman = 0.005
            }
            else if (size == "htsmall") {
                sigman = 0.005
            }
            else {
                sigman = 0.05
            }
            u <- matrix(rnorm(p1 * m1, 0, sqrt(sigman)), p1, 
                m1)
            for (i in 1:m1) {
                block_start <- (i - 1) * 5 + 1
                u[block_start:(block_start + 4), i] <- rnorm(5, 
                  0, sqrt(sigmau))
            }
            e1 <- matrix(rnorm(n1 * m1, 0, sqrt(sigma1)), n1, 
                m1)
            e2 <- rnorm(n2, 0, sqrt(sigma2))
            y <- wg1 \%*\% u + e1
            z <- wg2 \%*\% u \%*\% truealpha + e2
            M_step <- function(old, lambda, a_para, wg1, y, wg2, 
                u, truealpha, e2, m1, n1, n2, p1, k) {
                sigma1 <- old[1]
                sigma2 <- old[2]
                sigmau <- old[3]
                alpha_g <- old[-c(1:3)]
                a <- a_para
                abc <- (1 - a) * lambda * k
                ab <- a * lambda
                wg1_product <- t(wg1) \%*\% wg1
                wg2_product <- t(wg2) \%*\% wg2
                sigma_ui_constant <- diag(1, p1, p1)
                results <- mclapply(1:m1, function(i) {
                  sigma_ui_i <- solve(1/sigma1 * wg1_product + 
                    (alpha_g[i]^2)/sigma2 * wg2_product + 1/sigmau * 
                    sigma_ui_constant)
                  mu_ui_i <- sigma_ui_i \%*\% (1/sigma1 * t(wg1) \%*\% 
                    y[, i] + alpha_g[i]/sigma2 * t(wg2) \%*\% (wg2 \%*\% 
                    u[, i] \%*\% truealpha[i] + e2))
                  si <- 2 * t(wg2 \%*\% u[, i] \%*\% truealpha[i] + 
                    e2) \%*\% wg2 \%*\% mu_ui_i - a * lambda * sign(alpha_g[i])
                  abc <- (1 - a) * lambda * k[i]
                  ab <- a * lambda
                  alphaest <- ifelse(abs(si) <= abc | abs(si) <= 
                    ab, 0, (1 - abc/(abs(si))) * (2 * t(mu_ui_i) \%*\% 
                    t(wg2) \%*\% wg2 \%*\% mu_ui_i + 2 * sum(diag(wg2 \%*\% 
                    sigma_ui_i \%*\% t(wg2))))^(-1) * si)
                  newz <- wg2 \%*\% u[, i] * truealpha[i] + e2
                  E1 <- t(y[, i]) \%*\% y[, i] - 2 * t(y[, i]) \%*\% 
                    wg1 \%*\% mu_ui_i + t(mu_ui_i) \%*\% t(wg1) \%*\% 
                    wg1 \%*\% mu_ui_i + sum(diag(t(wg1) \%*\% wg1 \%*\% 
                    sigma_ui_i))
                  E2 <- t(newz) \%*\% newz - 2 * alphaest * t(newz) \%*\% 
                    wg2 \%*\% mu_ui_i + t(as.numeric(alphaest) * 
                    wg2 \%*\% mu_ui_i) \%*\% (as.numeric(alphaest) * 
                    wg2 \%*\% mu_ui_i) + sum(diag(as.numeric(alphaest)^2 * 
                    wg2 \%*\% sigma_ui_i \%*\% t(wg2)))
                  E3 <- t(mu_ui_i) \%*\% mu_ui_i + sum(diag(sigma_ui_i))
                  list(E1 = E1, E2 = E2, E3 = E3, alphaest = alphaest, 
                    mu_ui_i = mu_ui_i)
                }, mc.cores = 30)
                sigma1est <- 1/(m1 * n1) * sum(sapply(results, 
                  function(x) x$E1))
                sigma2est <- 1/(m1 * n2) * sum(sapply(results, 
                  function(x) x$E2))
                sigmauest <- 1/(m1 * p1) * sum(sapply(results, 
                  function(x) x$E3))
                log_likelihood_y = -0.5 * n1 * m1 * log(2 * pi * 
                  sigma1est) - 0.5 * sum(sapply(results, function(x) x$E1))/sigma1est
                log_likelihood_z = -0.5 * n2 * log(2 * pi * sigma2est) - 
                  0.5 * sum(sapply(results, function(x) x$E2))/sigma2est
                log_likelihood_u = -0.5 * p1 * m1 * log(2 * pi * 
                  sigmauest) - 0.5 * sum(sapply(results, function(x) x$E3))/sigmauest
                log_likelihood = log_likelihood_y + log_likelihood_z + 
                  log_likelihood_u
                print(c(sigma1est, sigma2est, sigmauest))
                return(list(est = c(sigma1est, sigma2est, sigmauest, 
                  unlist(sapply(results, function(x) x$alphaest))), 
                  lik = log_likelihood))
            }
            EM <- function(init, lambda, a_para, wg1, y, wg2, 
                u, truealpha, e2, m1, n1, n2, p1, k, max_iter = 10, 
                tol = 0.001) {
                theta_old <- init - 1e-04
                iter <- 1
                prev_diff <- Inf
                for (i in 1:max_iter) {
                  step_results <- M_step(old = theta_old, lambda = lambda, 
                    a_para = a_para, wg1 = wg1, y = y, wg2 = wg2, 
                    u = u, truealpha = truealpha, e2 = e2, m1 = m1, 
                    n1 = n1, n2 = n2, p1 = p1, k = k)
                  theta_new <- step_results$est
                  log_likelihood_new <- step_results$lik
                  current_diff <- mean(abs(theta_new - theta_old))
                  cat("Iteration:", iter, "parameter difference:", 
                    current_diff, "\n")
                  if (abs(current_diff - prev_diff) < tol) 
                    break
                  theta_old <- theta_new
                  prev_diff <- current_diff
                  iter <- iter + 1
                  print(prev_diff)
                }
                return(list(final_parameters = c(theta_new, lambda, 
                  log_likelihood_new)))
            }
            k_fold <- function(kf, a_values, lam_values, y, wg1, 
                e1, wg2, e2, u, m1, truealpha, p1, k) {
                fold_size_y <- dim(y)[1]/kf
                fold_size_wg2 <- dim(wg2)[1]/kf
                fold_indices <- lapply(1:kf, function(i) {
                  list(y_indices = ((i - 1) * fold_size_y + 1):(i * 
                    fold_size_y), wg2_indices = ((i - 1) * fold_size_wg2 + 
                    1):(i * fold_size_wg2))
                })
                all_errors <- list()
                for (i in lam_values) {
                  all_errors[[paste0("lambda_", i)]] <- mclapply(1:kf, 
                    function(fold) {
                      training_indices_y <- unlist(lapply(fold_indices[-fold], 
                        `[[`, "y_indices"))
                      testing_indices_y <- fold_indices[[fold]]$y_indices
                      training_indices_wg2 <- unlist(lapply(fold_indices[-fold], 
                        `[[`, "wg2_indices"))
                      testing_indices_wg2 <- fold_indices[[fold]]$wg2_indices
                      y_train <- y[training_indices_y, ]
                      y_test <- y[testing_indices_y, ]
                      wg1_train <- wg1[training_indices_y, ]
                      e1_train <- e1[training_indices_y, ]
                      wg2_train <- wg2[training_indices_wg2, 
                        ]
                      e2_train <- e2[training_indices_wg2]
                      wg2_test <- wg2[testing_indices_wg2, ]
                      e2_test <- e2[testing_indices_wg2]
                      sapply(a_values, function(a) {
                        initial = c(sigma1, sigma2, sigmau, truealpha)
                        result <- EM(initial, lambda = i, a_para = a, 
                          wg1 = wg1_train, y = y_train, wg2 = wg2_train, 
                          u = u, truealpha = truealpha, e2 = e2_train, 
                          m1 = m1, n1 = nrow(y_train), n2 = nrow(wg2_train), 
                          p1 = p1, k = k)
                        alpha_est <- result$final_parameters[-c(1:3, 
                          length(result$final_parameters))]
                        mean((z_est - z_test)^2)
                      })
                    }, mc.cores = 30, mc.preschedule = FALSE)
                }
                col_means <- lapply(all_errors, function(lambda_errors) {
                  sapply(1:length(lambda_errors[[1]]), function(col) {
                    mean(sapply(lambda_errors, function(fold_errors) fold_errors[col]))
                  })
                })
                optimal_a <- a_values[sapply(col_means, which.min)]
                return(list(a = optimal_a, colmean = col_means))
            }
            msecal <- function(wg2, mle) {
                mse <- rep(0, length(lam_values))
                for (i in 1:length(lam_values)) {
                  alphaest0 <- mle[i, 4:(m1 + 3)]
                  z_est0 <- wg2 \%*\% u \%*\% alphaest0 + e2
                  z <- wg2 \%*\% u \%*\% truealpha + e2
                  mse[i] <- mean((z - z_est0)^2)
                }
                return(mse)
            }
            if (size == "small") {
                a_mle <- k_fold(kf = 5, a_values = seq(0, 1, 
                  length = 5), lam_values = c(seq(0, 20, length = 7), 
                  seq(230, 500, length = 6), 600, 750), y, wg1, 
                  e1, wg2, e2, u, m1, truealpha, p1, k)$a
                lam_values_l0 = c(0, 10, 15, 20, seq(30, 200, 
                  length = 9), 400, 700)
                lam_values_l1 = c(0, 5, 10, 14, seq(17, 200, 
                  length = 11))
                lam_values_l2 = c(0, 15, 30, seq(40, 250, length = 10), 
                  350, 550)
            }
            else if (size == "large") {
                a_mle <- k_fold(kf = 5, a_values = seq(0, 1, 
                  length = 5), lam_values = c(seq(0, 20, length = 7), 
                  seq(230, 500, length = 6), 600, 750), y, wg1, 
                  e1, wg2, e2, u, m1, truealpha, p1, k)$a
                lam_values_l0 = c(0, 5, 15, seq(20, 300, length = 10), 
                  400, 500)
                lam_values_l1 = c(0, 3, 7, seq(15, 200, length = 10), 
                  300, 400)
                lam_values_l2 = c(0, 5, 10, seq(20, 300, length = 10), 
                  400, 500)
            }
            else if (size == "hcsmall") {
                a_mle <- k_fold(kf = 5, a_values = seq(0, 1, 
                  length = 5), lam_values = lam_values_l0, y, 
                  wg1, e1, wg2, e2, u, m1, truealpha, p1, k)$a
                lam_values_l0 = c(0, 1, 2, 3, seq(3.5, 40, length = 10), 
                  50)
                lam_values_l1 = c(0, 1, 2, 3, 5, seq(6, 25, length = 9), 
                  50)
                lam_values_l2 = c(0, 3, 6, 10, seq(12, 50, length = 10), 
                  100)
            }
            else if (size == "htsmall") {
                a_mle <- k_fold(kf = 5, a_values = seq(0, 1, 
                  length = 5), lam_values = lam_values_l0, y, 
                  wg1, e1, wg2, e2, u, m1, truealpha, p1, k)$a
                lam_values_l0 = c(0, 3, 7, 10, seq(15, 60, length = 10), 
                  100)
                lam_values_l1 = c(seq(0, 10, length = 11), 12, 
                  15, 20, 200)
                lam_values_l2 = c(0, 3, 6, seq(8, 33, length = 10), 
                  50, 100)
            }
            mle_l0 <- foreach(i = seq_along(lam_values_l0), .combine = rbind) \%do\% 
                {
                  lambda <- lam_values_l0[i]
                  a_para <- a_mle[i]
                  EM(init = c(sigma1, sigma2, sigmau, truealpha), 
                    lambda = lambda, a_para = a_para, wg1, y, 
                    wg2, u, truealpha, e2, m1, n1, n2, p1, k)$final_parameters
                }
            mle_l1 <- foreach(i = lam_values_l1, .combine = rbind) \%do\% 
                EM(init = c(sigma1, sigma2, sigmau, truealpha), 
                  lambda = i, a_para = 0, wg1, y, wg2, u, truealpha, 
                  e2, m1, n1, n2, p1, k)$final_parameters
            mle_l2 <- foreach(i = lam_values_l2, .combine = rbind) \%do\% 
                EM(init = c(sigma1, sigma2, sigmau, truealpha), 
                  lambda = i, a_para = 1, wg1, y, wg2, u, truealpha, 
                  e2, m1, n1, n2, p1, k)$final_parameters
            newu <- matrix(nrow = 5, ncol = m1)
            mle_comm <- matrix(NA, m1, 4)
            p_values_lrt <- numeric(length = m1)
            M_step_comm <- function(old, w1, y, w2, z, n1, n2) {
                sigma1 <- old[1]
                sigma2 <- old[2]
                sigmau <- old[3]
                alpha_g <- old[-c(1:3)]
                sigma_ui_i <- solve(1/sigma1 * t(w1) \%*\% w1 + 
                  (alpha_g^2)/sigma2 * t(w2) \%*\% w2 + 1/sigmau * 
                  diag(1, 5, 5))
                mu_ui_i <- sigma_ui_i \%*\% (1/sigma1 * t(w1) \%*\% 
                  y + alpha_g/sigma2 * t(w2) \%*\% z)
                alphaest <- solve(t(mu_ui_i) \%*\% t(w2) \%*\% w2 \%*\% 
                  mu_ui_i + sum(diag(w2 \%*\% sigma_ui_i \%*\% t(w2)))) * 
                  t(z) \%*\% w2 \%*\% mu_ui_i
                E1 <- t(y) \%*\% y - 2 * t(y) \%*\% w1 \%*\% mu_ui_i + 
                  t(mu_ui_i) \%*\% t(w1) \%*\% w1 \%*\% mu_ui_i
                E2 <- t(z) \%*\% z - 2 * alphaest * t(z) \%*\% w2 \%*\% 
                  mu_ui_i + alphaest^2 * t(mu_ui_i) \%*\% t(w2) \%*\% 
                  w2 \%*\% mu_ui_i
                E3 <- t(mu_ui_i) \%*\% mu_ui_i
                sigma1est <- 1/n1 * (E1 + sum(diag(t(w1) \%*\% 
                  w1 \%*\% sigma_ui_i)))
                sigma2est <- 1/n2 * (E2 + alphaest^2 * sum(diag(t(w2) \%*\% 
                  w2 \%*\% sigma_ui_i)))
                sigmauest <- 1/5 * (E3 + sum(diag(sigma_ui_i)))
                log_likelihood_y <- -0.5 * n1 * log(2 * pi * 
                  sigma1) - 0.5 * E1/sigma1
                log_likelihood_z <- -0.5 * n2 * log(2 * pi * 
                  sigma2) - 0.5 * E2/sigma2
                log_likelihood_u <- -0.5 * 5 * log(2 * pi * sigmau) - 
                  0.5 * E3/sigmau
                total_log_likelihood <- log_likelihood_y + log_likelihood_z + 
                  log_likelihood_u
                return(list(esti = c(sigma1est, sigma2est, sigmauest, 
                  alphaest), log_likelihood = total_log_likelihood))
            }
            EM_comm <- function(init, w1, y, w2, z, n1, n2, max_iter = 20, 
                tol = 0.001) {
                theta_old <- init - 1e-04
                iter <- 1
                prev_diff <- Inf
                for (i in 1:max_iter) {
                  step_result <- M_step_comm(old = theta_old, 
                    w1 = w1, y = y, w2 = w2, z = z, n1 = n1, 
                    n2 = n2)
                  theta_new <- step_result$esti
                  current_diff <- mean(abs(theta_new - theta_old))
                  if (abs(current_diff - prev_diff) < tol) 
                    break
                  theta_old <- theta_new
                  prev_diff <- current_diff
                  iter <- iter + 1
                  print(prev_diff)
                }
                return(list(final_parameters = theta_new))
            }
            loglik <- function(para, w1, y, w2, z, n1, n2) {
                sigma1 <- para[1]
                sigma2 <- para[2]
                sigmau <- para[3]
                alpha_g <- para[-c(1:3)]
                sigma_ui_i <- solve(1/sigma1 * t(w1) \%*\% w1 + 
                  (alpha_g^2)/sigma2 * t(w2) \%*\% w2 + 1/sigmau * 
                  diag(1, 5, 5))
                mu_ui_i <- sigma_ui_i \%*\% (1/sigma1 * t(w1) \%*\% 
                  y + alpha_g/sigma2 * t(w2) \%*\% z)
                E1 <- t(y) \%*\% y - 2 * t(y) \%*\% w1 \%*\% mu_ui_i + 
                  t(mu_ui_i) \%*\% t(w1) \%*\% w1 \%*\% mu_ui_i
                E2 <- t(z) \%*\% z - 2 * alpha_g * t(z) \%*\% w2 \%*\% 
                  mu_ui_i + alpha_g^2 * t(mu_ui_i) \%*\% t(w2) \%*\% 
                  w2 \%*\% mu_ui_i
                E3 <- t(mu_ui_i) \%*\% mu_ui_i
                log_likelihood_y <- -0.5 * n1 * log(2 * pi * 
                  sigma1) - 0.5 * E1/sigma1
                log_likelihood_z <- -0.5 * n2 * log(2 * pi * 
                  sigma2) - 0.5 * E2/sigma2
                log_likelihood_u <- -0.5 * 5 * log(2 * pi * sigmau) - 
                  0.5 * E3/sigmau
                total_log_likelihood <- log_likelihood_y + log_likelihood_z + 
                  log_likelihood_u
                return(total_log_likelihood)
            }
            for (i in 1:m1) {
                row_indices <- ((i - 1) * 5 + 1):(i * 5)
                newu[, i] <- u[row_indices, i]
                w1 <- wg1[, ((i - 1) * 5 + 1):(5 * i)]
                w2 <- wg2[, ((i - 1) * 5 + 1):(5 * i)]
                ynew <- w1 \%*\% as.matrix(newu[, i]) + e1[, i]
                znew <- truealpha[i] * w2 \%*\% as.matrix(newu[, 
                  i]) + e2
                em1 <- EM_comm(init = c(sigma1, sigma2, sigmau, 
                  truealpha[i]), w1 = w1, y = ynew, w2 = w2, 
                  z = znew, n1 = n1, n2 = n2)
                mle_comm[i, ] <- em1$final_parameters
                log_likelihood_alt <- loglik(c(mle_comm[i, 1], 
                  mle_comm[i, 2], mle_comm[i, 3], mle_comm[i, 
                    4]), w1 = w1, y = ynew, w2 = w2, z = znew, 
                  n1 = n1, n2 = n2)
                log_likelihood_null <- loglik(c(mle_comm[i, 1], 
                  mle_comm[i, 2], mle_comm[i, 3], 0), w1 = w1, 
                  y = ynew, w2 = w2, z = znew, n1 = n1, n2 = n2)
                lr_statistic <- 2 * (log_likelihood_alt - log_likelihood_null)
                p_values_lrt[i] <- pchisq(lr_statistic, df = 1, 
                  lower.tail = FALSE)
            }
            mle_comm <- c(colMeans(mle_comm[, 1:3]), mle_comm[, 
                4])
            library(glmnet)
            library(foreach)
            coefficients_en <- list()
            association_results <- vector("list", length = m1)
            p_values_EN <- rep(NA, m1)
            for (i in 1:m1) {
                row_indices <- ((i - 1) * 5 + 1):(i * 5)
                u_i <- u[row_indices, i]
                w1_i <- wg1[, row_indices]
                w2_i <- wg2[, row_indices]
                y_i <- w1_i \%*\% u_i + e1[, i]
                z_i <- truealpha[i] * w2_i \%*\% u_i + e2
                cv_fit <- cv.glmnet(w1_i, y_i, alpha = 0, family = "gaussian")
                best_lambda <- cv_fit$lambda.min
                fit <- glmnet(w1_i, y_i, alpha = 0, lambda = best_lambda)
                coefficients_en[[i]] <- coef(fit, s = "lambda.min")
                predicted_expression <- predict(fit, newx = w2_i)
                association_results[[i]] <- lm(z_i ~ predicted_expression)
                if (dim(summary(association_results[[i]])$coefficients)[1] > 
                  1) {
                  p_values_EN[i] <- summary(association_results[[i]])$coefficients[2, 
                    "Pr(>|t|)"]
                }
                else {
                  p_values_EN[i] <- NA
                }
            }
            coefficients_lasso <- list()
            association_results <- vector("list", length = m1)
            p_values_lasso <- rep(NA, m1)
            for (i in 1:m1) {
                row_indices <- ((i - 1) * 5 + 1):(i * 5)
                u_i <- u[row_indices, i]
                w1_i <- wg1[, row_indices]
                w2_i <- wg2[, row_indices]
                y_i <- w1_i \%*\% u_i + e1[, i]
                z_i <- truealpha[i] * w2_i \%*\% u_i + e2
                cv_fit <- cv.glmnet(w1_i, y_i, alpha = 0.9, family = "gaussian")
                best_lambda <- cv_fit$lambda.min
                fit <- glmnet(w1_i, y_i, alpha = 0.9, lambda = best_lambda)
                coefficients_lasso[[i]] <- coef(fit, s = "lambda.min")
                predicted_expression <- predict(fit, newx = w2_i)
                association_results[[i]] <- lm(z_i ~ predicted_expression)
                if (dim(summary(association_results[[i]])$coefficients)[1] > 
                  1) {
                  p_values_lasso[i] <- summary(association_results[[i]])$coefficients[2, 
                    "Pr(>|t|)"]
                }
                else {
                  p_values_lasso[i] <- NA
                }
            }
            u_i_est <- matrix(0, ncol = p1, nrow = m1)
            library(MASS)
            for (i in 1:m1) {
                y_i <- y[, i]
                u_i_est[i, ] <- ginv(t(wg1) \%*\% wg1) \%*\% t(wg1) \%*\% 
                  y_i
            }
            predictor <- wg2 \%*\% t(u_i_est)
            model <- lm(z ~ predictor + 0)
            mle_multi <- coef(model)
            p_values_multi <- summary(model)$coefficients[, 4]
            residuals <- model$residuals
            mse_multi <- mean(residuals^2)
            result_list <- list(mle_l0 = mle_l0, mle_l1 = mle_l1, 
                mle_l2 = mle_l2, mle_comm = mle_comm, mle_multi = mle_multi, 
                p_values_comm = p_values_lrt, p_values_EN = p_values_EN, 
                p_values_lasso = p_values_lasso, p_values_multi = p_values_multi, 
                mse = rbind(msecal(wg2, mle_l0), msecal(wg2, 
                  mle_l1), msecal(wg2, mle_l2)), mse_multi = mse_multi)
            return(result_list)
        }
    }
    else if (size == "large") {
        a_mle <- c(0, 0.5, rep(1, 13))
        lam_values_l0 = c(0, 5, 15, seq(20, 300, length = 10), 
            400, 500)
        lam_values_l1 = c(0, 3, 7, seq(15, 200, length = 10), 
            300, 400)
        lam_values_l2 = c(0, 5, 10, seq(20, 300, length = 10), 
            400, 500)
    }
    else if (size == "hcsmall") {
        a_mle <- c(0, 0.25, 0.75, rep(1, 12))
        lam_values_l0 = c(0, 1, 2, 3, seq(3.5, 40, length = 10), 
            50)
        lam_values_l1 = c(0, 1, 2, 3, 5, seq(6, 25, length = 9), 
            50)
        lam_values_l2 = c(0, 3, 6, 10, seq(12, 50, length = 10), 
            100)
    }
    else if (size == "htsmall") {
        a_mle <- c(0, 0.5, 0.75, rep(1, 12))
        lam_values_l0 = c(0, 3, 7, 10, seq(15, 60, length = 10), 
            100)
        lam_values_l1 = c(seq(0, 10, length = 11), 12, 15, 20, 
            200)
        lam_values_l2 = c(0, 3, 6, seq(8, 33, length = 10), 50, 
            100)
    }
    mle_l0 <- foreach(i = seq_along(lam_values_l0), .combine = rbind) \%do\% 
        {
            lambda <- lam_values_l0[i]
            a_para <- a_mle[i]
            EM(init = c(sigma1, sigma2, sigmau, truealpha), lambda = lambda, 
                a_para = a_para, wg1, y, wg2, u, truealpha, e2, 
                m1, n1, n2, p1, k)$final_parameters
        }
    mle_l1 <- foreach(i = lam_values_l1, .combine = rbind) \%do\% 
        EM(init = c(sigma1, sigma2, sigmau, truealpha), lambda = i, 
            a_para = 0, wg1, y, wg2, u, truealpha, e2, m1, n1, 
            n2, p1, k)$final_parameters
    mle_l2 <- foreach(i = lam_values_l2, .combine = rbind) \%do\% 
        EM(init = c(sigma1, sigma2, sigmau, truealpha), lambda = i, 
            a_para = 1, wg1, y, wg2, u, truealpha, e2, m1, n1, 
            n2, p1, k)$final_parameters
    newu <- matrix(nrow = 5, ncol = m1)
    mle_comm <- matrix(NA, m1, 4)
    p_values_lrt <- numeric(length = m1)
    M_step_comm <- function(old, w1, y, w2, z, n1, n2) {
        sigma1 <- old[1]
        sigma2 <- old[2]
        sigmau <- old[3]
        alpha_g <- old[-c(1:3)]
        sigma_ui_i <- solve(1/sigma1 * t(w1) \%*\% w1 + (alpha_g^2)/sigma2 * 
            t(w2) \%*\% w2 + 1/sigmau * diag(1, 5, 5))
        mu_ui_i <- sigma_ui_i \%*\% (1/sigma1 * t(w1) \%*\% y + alpha_g/sigma2 * 
            t(w2) \%*\% z)
        alphaest <- solve(t(mu_ui_i) \%*\% t(w2) \%*\% w2 \%*\% mu_ui_i + 
            sum(diag(w2 \%*\% sigma_ui_i \%*\% t(w2)))) * t(z) \%*\% 
            w2 \%*\% mu_ui_i
        E1 <- t(y) \%*\% y - 2 * t(y) \%*\% w1 \%*\% mu_ui_i + t(mu_ui_i) \%*\% 
            t(w1) \%*\% w1 \%*\% mu_ui_i
        E2 <- t(z) \%*\% z - 2 * alphaest * t(z) \%*\% w2 \%*\% mu_ui_i + 
            alphaest^2 * t(mu_ui_i) \%*\% t(w2) \%*\% w2 \%*\% mu_ui_i
        E3 <- t(mu_ui_i) \%*\% mu_ui_i
        sigma1est <- 1/n1 * (E1 + sum(diag(t(w1) \%*\% w1 \%*\% sigma_ui_i)))
        sigma2est <- 1/n2 * (E2 + alphaest^2 * sum(diag(t(w2) \%*\% 
            w2 \%*\% sigma_ui_i)))
        sigmauest <- 1/5 * (E3 + sum(diag(sigma_ui_i)))
        log_likelihood_y <- -0.5 * n1 * log(2 * pi * sigma1) - 
            0.5 * E1/sigma1
        log_likelihood_z <- -0.5 * n2 * log(2 * pi * sigma2) - 
            0.5 * E2/sigma2
        log_likelihood_u <- -0.5 * 5 * log(2 * pi * sigmau) - 
            0.5 * E3/sigmau
        total_log_likelihood <- log_likelihood_y + log_likelihood_z + 
            log_likelihood_u
        return(list(esti = c(sigma1est, sigma2est, sigmauest, 
            alphaest), log_likelihood = total_log_likelihood))
    }
    EM_comm <- function(init, w1, y, w2, z, n1, n2, max_iter = 20, 
        tol = 0.001) {
        theta_old <- init - 1e-04
        iter <- 1
        prev_diff <- Inf
        for (i in 1:max_iter) {
            step_result <- M_step_comm(old = theta_old, w1 = w1, 
                y = y, w2 = w2, z = z, n1 = n1, n2 = n2)
            theta_new <- step_result$esti
            current_diff <- mean(abs(theta_new - theta_old))
            if (abs(current_diff - prev_diff) < tol) 
                break
            theta_old <- theta_new
            prev_diff <- current_diff
            iter <- iter + 1
            print(prev_diff)
        }
        return(list(final_parameters = theta_new))
    }
    loglik <- function(para, w1, y, w2, z, n1, n2) {
        sigma1 <- para[1]
        sigma2 <- para[2]
        sigmau <- para[3]
        alpha_g <- para[-c(1:3)]
        sigma_ui_i <- solve(1/sigma1 * t(w1) \%*\% w1 + (alpha_g^2)/sigma2 * 
            t(w2) \%*\% w2 + 1/sigmau * diag(1, 5, 5))
        mu_ui_i <- sigma_ui_i \%*\% (1/sigma1 * t(w1) \%*\% y + alpha_g/sigma2 * 
            t(w2) \%*\% z)
        E1 <- t(y) \%*\% y - 2 * t(y) \%*\% w1 \%*\% mu_ui_i + t(mu_ui_i) \%*\% 
            t(w1) \%*\% w1 \%*\% mu_ui_i
        E2 <- t(z) \%*\% z - 2 * alpha_g * t(z) \%*\% w2 \%*\% mu_ui_i + 
            alpha_g^2 * t(mu_ui_i) \%*\% t(w2) \%*\% w2 \%*\% mu_ui_i
        E3 <- t(mu_ui_i) \%*\% mu_ui_i
        log_likelihood_y <- -0.5 * n1 * log(2 * pi * sigma1) - 
            0.5 * E1/sigma1
        log_likelihood_z <- -0.5 * n2 * log(2 * pi * sigma2) - 
            0.5 * E2/sigma2
        log_likelihood_u <- -0.5 * 5 * log(2 * pi * sigmau) - 
            0.5 * E3/sigmau
        total_log_likelihood <- log_likelihood_y + log_likelihood_z + 
            log_likelihood_u
        return(total_log_likelihood)
    }
    for (i in 1:m1) {
        row_indices <- ((i - 1) * 5 + 1):(i * 5)
        newu[, i] <- u[row_indices, i]
        w1 <- wg1[, ((i - 1) * 5 + 1):(5 * i)]
        w2 <- wg2[, ((i - 1) * 5 + 1):(5 * i)]
        ynew <- w1 \%*\% as.matrix(newu[, i]) + e1[, i]
        znew <- truealpha[i] * w2 \%*\% as.matrix(newu[, i]) + 
            e2
        em1 <- EM_comm(init = c(sigma1, sigma2, sigmau, truealpha[i]), 
            w1 = w1, y = ynew, w2 = w2, z = znew, n1 = n1, n2 = n2)
        mle_comm[i, ] <- em1$final_parameters
        log_likelihood_alt <- loglik(c(mle_comm[i, 1], mle_comm[i, 
            2], mle_comm[i, 3], mle_comm[i, 4]), w1 = w1, y = ynew, 
            w2 = w2, z = znew, n1 = n1, n2 = n2)
        log_likelihood_null <- loglik(c(mle_comm[i, 1], mle_comm[i, 
            2], mle_comm[i, 3], 0), w1 = w1, y = ynew, w2 = w2, 
            z = znew, n1 = n1, n2 = n2)
        lr_statistic <- 2 * (log_likelihood_alt - log_likelihood_null)
        p_values_lrt[i] <- pchisq(lr_statistic, df = 1, lower.tail = FALSE)
    }
    mle_comm <- c(colMeans(mle_comm[, 1:3]), mle_comm[, 4])
    library(glmnet)
    library(foreach)
    coefficients_en <- list()
    association_results <- vector("list", length = m1)
    p_values_EN <- rep(NA, m1)
    for (i in 1:m1) {
        row_indices <- ((i - 1) * 5 + 1):(i * 5)
        u_i <- u[row_indices, i]
        w1_i <- wg1[, row_indices]
        w2_i <- wg2[, row_indices]
        y_i <- w1_i \%*\% u_i + e1[, i]
        z_i <- truealpha[i] * w2_i \%*\% u_i + e2
        cv_fit <- cv.glmnet(w1_i, y_i, alpha = 0, family = "gaussian")
        best_lambda <- cv_fit$lambda.min
        fit <- glmnet(w1_i, y_i, alpha = 0, lambda = best_lambda)
        coefficients_en[[i]] <- coef(fit, s = "lambda.min")
        predicted_expression <- predict(fit, newx = w2_i)
        association_results[[i]] <- lm(z_i ~ predicted_expression)
        if (dim(summary(association_results[[i]])$coefficients)[1] > 
            1) {
            p_values_EN[i] <- summary(association_results[[i]])$coefficients[2, 
                "Pr(>|t|)"]
        }
        else {
            p_values_EN[i] <- NA
        }
    }
    coefficients_lasso <- list()
    association_results <- vector("list", length = m1)
    p_values_lasso <- rep(NA, m1)
    for (i in 1:m1) {
        row_indices <- ((i - 1) * 5 + 1):(i * 5)
        u_i <- u[row_indices, i]
        w1_i <- wg1[, row_indices]
        w2_i <- wg2[, row_indices]
        y_i <- w1_i \%*\% u_i + e1[, i]
        z_i <- truealpha[i] * w2_i \%*\% u_i + e2
        cv_fit <- cv.glmnet(w1_i, y_i, alpha = 0.9, family = "gaussian")
        best_lambda <- cv_fit$lambda.min
        fit <- glmnet(w1_i, y_i, alpha = 0.9, lambda = best_lambda)
        coefficients_lasso[[i]] <- coef(fit, s = "lambda.min")
        predicted_expression <- predict(fit, newx = w2_i)
        association_results[[i]] <- lm(z_i ~ predicted_expression)
        if (dim(summary(association_results[[i]])$coefficients)[1] > 
            1) {
            p_values_lasso[i] <- summary(association_results[[i]])$coefficients[2, 
                "Pr(>|t|)"]
        }
        else {
            p_values_lasso[i] <- NA
        }
    }
    u_i_est <- matrix(0, ncol = p1, nrow = m1)
    library(MASS)
    for (i in 1:m1) {
        y_i <- y[, i]
        u_i_est[i, ] <- ginv(t(wg1) \%*\% wg1) \%*\% t(wg1) \%*\% y_i
    }
    predictor <- wg2 \%*\% t(u_i_est)
    model <- lm(z ~ predictor + 0)
    mle_multi <- coef(model)
    p_values_multi <- summary(model)$coefficients[, 4]
    residuals <- model$residuals
    mse_multi <- mean(residuals^2)
    result_list <- list(mle_l0 = mle_l0, mle_l1 = mle_l1, mle_l2 = mle_l2, 
        mle_comm = mle_comm, mle_multi = mle_multi, p_values_comm = p_values_lrt, 
        p_values_EN = p_values_EN, p_values_lasso = p_values_lasso, 
        p_values_multi = p_values_multi, mse = rbind(msecal(wg2, 
            mle_l0), msecal(wg2, mle_l1), msecal(wg2, mle_l2)), 
        mse_multi = mse_multi)
    return(result_list)
  }
}

