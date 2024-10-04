\name{TIPS}
\alias{TIPS}
\title{Implementation of the TIPS Algorithm for Gene Expression Analysis}
\description{
The \code{TIPS} function applies the TIPS (Two-stage Integrative Polygenic Score) algorithm to estimate genetic effects using simulated or real gene expression and SNP data. This function handles different inheritance models and performs multi-lambda EM algorithm computations.
}
\usage{
TIPS(simulated_data)
}
\arguments{
  \item{simulated_data}{A list containing simulated genotype, gene expression, and phenotype data. The data must be generated using the \code{\link{simulate_data}} function.}
}
\details{
This function first extracts the genotype (wg1, wg2), gene expression (y), and phenotype (z) data from the input list. It applies the EM (Expectation-Maximization) algorithm to estimate effect sizes for each gene's SNP block, and evaluates multiple lambda values to regularize the model. The method supports different inheritance models (additive, heterogeneous, recessive, compensatory) and computes MSE values to assess the performance of the model.
}
\value{
  A list containing:
  \item{mle_l0}{Estimated parameters for the l0 penalty model (multi-lambda EM).}
  \item{mle_l1}{Estimated parameters for the l1 penalty model (multi-lambda EM).}
  \item{mle_l2}{Estimated parameters for the l2 penalty model (multi-lambda EM).}
  \item{mse_l0}{MSE for the l0 penalty model.}
  \item{mse_l1}{MSE for the l1 penalty model.}
  \item{mse_l2}{MSE for the l2 penalty model.}
}
\examples{
simulated_data <- simulate_data(300, 600, 180, 900, rep(1, 180), 0.1, 0.2, 0.05, runif(180, 0.5, 1.5))
result <- TIPS(simulated_data)
}
\seealso{
  \code{\link{simulate_data}} for generating the data.
}
\author{
Neng
}
