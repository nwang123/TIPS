# This simulation generates gene expression or SNP data for reference panel (`n1 = 300`) and GWAS data (`n2 = 600`), and it simulates 180 genes with varying numbers of SNPs (5 SNPs per gene). The true effect sizes for these SNPs are controlled by the `truealpha` parameter.
### Parameters:
# - **n1**: Number of samples in the first group (e.g., 300 ).
# - **n2**: Number of samples in the second group (e.g., 600 ).
# - **m1**: Number of genes being simulated (180).
# - **p1**: Number of variables (5 SNPs per gene, hence 900 SNPs in total).
# - **k**: SNP loadings for each gene. Some groups have higher impact SNPs (`sqrt(25)`), while others have smaller effects.
# - **sigma1**, **sigma2**, **sigmau**: Variance components for noise in the data.
# - **truealpha**: A vector representing the true effect sizes of the SNPs.


# Define parameters
n1 <- 100  # Sample size for gene expression data
n2 <- 200  # Sample size for phenotype data
m1 <- 50   # Number of genes
p1 <- 5 * m1  # Number of SNPs
k <- c(rep(sqrt(25), 25 * 5), rep(1, 55))  # pathway grouping information
sigma1 <- 0.2  # Variance for gene expression
sigma2 <- 0.2  # Variance for phenotype
sigmau <- 0.02  # Variance for SNPs' effect on gene expression
truealpha <- rep(1, m1)  # True alpha values
size <- "small"  # Size parameter for variance
wg_str <- "additive"  # Inheritance model (additive)

# Step 1: Simulate the data
simulated_data <- simulate_data(n1 = n1, n2 = n2, m1 = m1, p1 = p1, 
                                sigma1 = sigma1, sigma2 = sigma2, 
                                sigmau = sigmau, truealpha = truealpha, 
                                size = size, wg_str = wg_str)

# Step 2: Run the TIPS algorithm
# Using the simulated data and the penalization parameter k
result <- tips(simulated_data, k = k, 
               sigma1_init = sigma1, sigma2_init = sigma2, 
               sigmau_init = sigmau, truealpha_init = truealpha)

# View the final estimated parameters from the TIPS algorithm
print(result)
