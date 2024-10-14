# Overview

# This repository contains an R function, `simulate_data()`, that generates synthetic genetic and gene expression data. It is designed for researchers and data scientists working in areas such as statistical genetics, genomics, and machine learning. The function simulates genotype matrices (`wg1`, `wg2`), gene expression data (`y`), and phenotype data (`z`) for downstream analysis.

## Features

# - **Genotype Simulation**: Generates genotype data under different inheritance models:
#     - Additive
#     - Heterogeneous
#     - Recessive
#     - Compensatory
# - **Variance Control**: Parameters to control the variances for both genotype data and noise, allowing the user to simulate different levels of genetic and environmental variation.
# - **Customizable Sizes**: Simulate datasets of various sizes and genetic architectures, including small and large datasets.
# - **Error Terms**: Customizable error terms to reflect biological variability in gene expression and phenotype data.

# Function Description

### `simulate_data()`

# This function simulates genotype and gene expression data, allowing the user to specify sample sizes, number of genes, SNPs, inheritance model, and other parameters.

#### Parameters:

# - `n1`, `n2`: Sample sizes for two groups.
# - `m1`: Number of genes to simulate.
# - `p1`: Number of SNPs to simulate.
# - `k`: SNP loadings, defining how many SNPs affect each gene.
# - `sigma1`, `sigma2`, `sigmau`: Variances for the error terms in gene expression and genotype matrices.
# - `truealpha`: True effect sizes used for simulating phenotypic data.
# - `size`: Size of the dataset (`small`, `large`, `hcsmall`, `htsmall`).
# - `wg_str`: Type of genetic inheritance model (`additive`, `heterogeneous`, `recessive`, `compensatory`).

#### Returns:

# The function returns a list containing the following elements:

# - `wg1`, `wg2`: Simulated genotype matrices.
# - `y`: Simulated gene expression data.
# - `z`: Simulated phenotype data.
# - Other relevant parameters (`n1`, `n2`, `m1`, `p1`, etc.).


simulate_data <- function(n1, n2, m1, p1, sigma1, sigma2, sigmau, truealpha,
                          size = c("small", "large", "hcsmall", "htsmall"),
                          wg_str = c("additive", "heterogeneous", "recessive", "compensatory")) {
  
  library(foreach)
  
  # Function to generate genotype data
  generatewg <- function(n) {
    sequ1 <- rbinom(n, 1, 0.2)  # Assume MAF=0.2
    sequ2 <- rbinom(n, 1, 0.2)
    allele <- sequ1 + sequ2
    return(allele)
  }
  
  # Generate genotype matrices wg1 and wg2 based on inheritance model
  if (wg_str == "additive") {
    wg1 <- foreach(i = 1:p1, .combine = cbind) %do% generatewg(n1)
    wg2 <- foreach(i = 1:p1, .combine = cbind) %do% generatewg(n2)
  } else if (wg_str == "heterogeneous") {
    wg1 <- foreach(i = 1:p1, .combine = cbind) %do% generatewg(n1)
    wg2 <- foreach(i = 1:p1, .combine = cbind) %do% generatewg(n2)
    wg1[wg1 == 2] <- 1
    wg2[wg2 == 2] <- 1
  } else if (wg_str == "recessive") {
    wg1 <- foreach(i = 1:p1, .combine = cbind) %do% generatewg(n1)
    wg2 <- foreach(i = 1:p1, .combine = cbind) %do% generatewg(n2)
    wg1[wg1 == 1] <- 0
    wg2[wg2 == 1] <- 0
  } else if (wg_str == "compensatory") {
    wg1 <- foreach(i = 1:p1, .combine = cbind) %do% generatewg(n1)
    wg2 <- foreach(i = 1:p1, .combine = cbind) %do% generatewg(n2)
    wg1[wg1 == 2] <- 0
    wg2[wg2 == 2] <- 0
  }
  
  # Normalize genotype matrices
  normalize_matrix <- function(mat) {
    min_vals <- apply(mat, 2, min)
    max_vals <- apply(mat, 2, max)
    return((mat - min_vals) / (max_vals - min_vals))
  }
  
  wg1 <- normalize_matrix(wg1)
  wg2 <- normalize_matrix(wg2)
  
  # Determine variance based on size
  if (size == "small") {
    sigman <- 0.05
  } else if (size == "hcsmall" || size == "htsmall") {
    sigman <- 0.005
  } else {
    sigman <- 0.05
  }
  
  # Initialize the matrix u_init
  u_init <- matrix(rnorm(p1 * m1, 0, sqrt(sigman)), p1, m1)
  for (i in 1:m1) {
    block_start <- (i - 1) * 5 + 1
    u_init[block_start:(block_start + 4), i] <- rnorm(5, 0, sqrt(sigmau))
  }
  
  # Generate error terms
  e1 <- matrix(rnorm(n1 * m1, 0, sqrt(sigma1)), n1, m1)
  e2 <- rnorm(n2, 0, sqrt(sigma2))
  
  # Generate response variables
  y <- wg1 %*% u_init + e1  # Gene expressions
  z <- wg2 %*% u_init %*% truealpha + e2  # Phenotype
  
  # Return the simulated data as a list
  simulated_data <- list(
    wg1 = wg1,
    wg2 = wg2,
    u_init = u_init,
    y = y,
    z = z
  )
  
  return(simulated_data)
}
