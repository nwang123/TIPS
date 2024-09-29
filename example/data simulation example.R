# Gene Expression Simulation Example

This repository contains an example of simulating gene expression data with various biological parameters using R. The code simulates data for two groups with differing sample sizes and includes variance parameters for noise and effect size simulation.
## Simulation Overview

This simulation generates gene expression or SNP data for two groups (`n1 = 300`, `n2 = 600`), and it simulates 180 genes with varying numbers of SNPs (5 SNPs per gene). The true effect sizes for these SNPs are controlled by the `truealpha` parameter.

### Parameters:
- **n1**: Number of samples in the first group (e.g., 300 ).
- **n2**: Number of samples in the second group (e.g., 600 ).
- **m1**: Number of genes being simulated (180).
- **p1**: Number of variables (5 SNPs per gene, hence 900 SNPs in total).
- **k**: SNP loadings for each gene. Some groups have higher impact SNPs (`sqrt(25)`), while others have smaller effects.
- **sigma1**, **sigma2**, **sigmau**: Variance components for noise in the data.
- **truealpha**: A vector representing the true effect sizes of the SNPs.

n1 = 300
n2 = 600
m1 = 180
p1 = 5 * m1
k = c(rep(sqrt(25), 25 * 5), rep(1, 55))
sigma1 = 0.2
sigma2 = 0.2
sigmau = 0.02
truealpha <- c(rep(0, 25), rep(0.2, 25), rep(0.02, 25), c(rep(0, 10), rep(0.2, 10), rep(0, 5)),
               c(rep(0, 10), rep(0.02, 5), rep(0, 10)), 
               c(rep(0.2, 10), rep(0.02, 10),  rep(0.2, 5), rep(0.02, 0), rep(0, 30)))

simu1_1 <- data_simu(n1, n2, m1, p1, k, sigma1, sigma2, sigmau, truealpha, size = "small", wg_str = "additive")
