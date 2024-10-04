\name{simulate_data}
\alias{simulate_data}
\title{Simulate Gene Expression and SNP Data}
\description{
  The \code{simulate_data} function simulates gene expression and SNP data for two groups with different inheritance models. The function supports multiple inheritance models such as additive, heterogeneous, recessive, and compensatory, and provides gene expression and phenotype data.
}
\usage{
simulate_data(n1, n2, m1, p1, k, sigma1, sigma2, sigmau, truealpha, 
              size = c("small", "large", "hcsmall", "htsmall"), 
              wg_str = c("additive", "heterogeneous", "recessive", "compensatory"))
}
\arguments{
  \item{n1}{Number of samples in the first group.}
  \item{n2}{Number of samples in the second group.}
  \item{m1}{Number of genes being simulated.}
  \item{p1}{Number of SNPs being simulated.}
  \item{k}{SNP loadings for each gene.}
  \item{sigma1}{Variance component for noise in gene expression data.}
  \item{sigma2}{Variance component for noise in phenotype data.}
  \item{sigmau}{Variance component for SNP effects.}
  \item{truealpha}{True effect sizes of the SNPs on the phenotype.}
  \item{size}{Size category for determining noise variance. Options are \code{"small"}, \code{"large"}, \code{"hcsmall"}, and \code{"htsmall"}.}
  \item{wg_str}{String specifying the inheritance model. Options are \code{"additive"}, \code{"heterogeneous"}, \code{"recessive"}, and \code{"compensatory"}.}
}
\value{
  A list containing the following elements:
  \item{wg1}{Genotype matrix for the first group.}
  \item{wg2}{Genotype matrix for the second group.}
  \item{u}{Effect size matrix for the SNPs.}
  \item{e1}{Noise matrix for the gene expression data.}
  \item{e2}{Noise vector for the phenotype data.}
  \item{y}{Simulated gene expression data for the first group.}
  \item{z}{Simulated phenotype data for the second group.}
}
\examples{
n1 <- 300
n2 <- 600
m1 <- 180
p1 <- 900
k <- rep(1, m1)
sigma1 <- 0.1
sigma2 <- 0.2
sigmau <- 0.05
truealpha <- runif(m1, 0.5, 1.5)
sim_data <- simulate_data(n1, n2, m1, p1, k, sigma1, sigma2, sigmau, truealpha,size = "small",wg_str ="additive")
}
\seealso{
  \code{\link{TIPS}} for further analysis on the simulated data.
}
\author{
Neng
}
