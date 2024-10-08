% This function implements the M-step of the Expectation-Maximization (EM) algorithm for estimating parameters
% in a genetic data model. It estimates the variance components (`sigma1`, `sigma2`, `sigmau`) and the effect sizes
% (`alpha`) for each gene by solving a series of linear equations.
%
% Inputs:
% - old: Initial estimates for the parameters (sigma1, sigma2, sigmau, alpha).
% - a_para: A regularization parameter controlling the amount of penalization.
% - wg1, wg2: Genotypic matrices for two groups (n1 for group 1, n2 for group 2).
% - z: Phenotypic data for group 2.
% - y: Gene expression data for group 1.
% - lambda: Regularization parameter for Lasso or Ridge-type shrinkage.
% - k: A vector representing the SNP loadings for each gene.
% - m1: The number of genes.
% - n1, n2: The number of samples in groups 1 and 2, respectively.
% - p1: The number of SNPs.
%
% Outputs:
% - A vector containing the estimated variance components (`sigma1`, `sigma2`, `sigmau`), the estimated effect sizes (`alphaest`), log-likelihood, and other auxiliary variables (e.g., E2).
%
% The function processes each gene's SNP block, computes the necessary matrix operations, and updates the parameter estimates using the M-step of the EM algorithm. The results are used to estimate residuals and log-likelihood for convergence checking in subsequent EM iterations.
%
% The function leverages parallel processing (via `parfor`) for efficiency, allowing for the concurrent computation of gene-specific parameters.

function [output] = mstep_updated(old, a_para, wg1, wg2, z, y, lambda, k, m1, n1, n2, p1)
% Initialize parameters
sigma2 = old(2);
sigmau = old(3);
alpha_g = old(4:end);
sigma1 = old(1);


% Initialize matrices and vectors
results = repmat(struct('E1', [], 'E2', [], 'E3', [], 'alphaest', []), m1, 1);
index_gene_geno = readtable('index_gene_geno.csv');
index_gene_array = table2array(index_gene_geno);
start_indices = [1, (cumsum(index_gene_array(1:end-1)) + 1)'];
end_indices = cumsum(index_gene_array);
% Process each gene's SNP block
parfor i = 1:m1
     disp(['Index: ' num2str(i) ' in m1 using mstep_updated function.']);
    % Define block indices for the current gene
    block_start = start_indices(i);
    block_end = end_indices(i);

    % Here you can use block_start and block_end to define your block_indices
    block_indices = block_start:block_end;
    % Extract blocks for wg1 and wg2
    wg1_block = wg1(:, block_indices);
    wg2_block = wg2(:, block_indices);

    % Compute products
    wg1_product = wg1_block' * wg1_block;
    wg2_product = wg2_block' * wg2_block;

    col = size(wg1_block)

    % Form matrix A
    A = (1 / sigma1) * wg1_product + (alpha_g(i)^2 / sigma2) * wg2_product + (1 / sigmau) * eye(col(2));

    % Regularization to ensure numerical stability
    A = A + 1e-4 * eye(col(2));

    % Compute the inverse of A
    sigma_ui_i = inv(A);

    % Compute mu_ui_i for the current gene
    mu_ui_i = sigma_ui_i * ((1 / sigma1) * wg1_block' * y(:, i) + (alpha_g(i) / sigma2) * wg2_block' * z);

    % Calculate si for significance testing
    si = 2 * z' * wg2_block * mu_ui_i - a_para * lambda * sign(alpha_g(i));
    abc_i = (1 - a_para) * lambda * k(i);

    % Determine alphaest based on si
    if abs(si) <= abc_i || abs(si) <= (a_para * lambda)
        alphaest = 0;
    else
        alphaest = (1 - abc_i / abs(si)) * (2 * mu_ui_i' * wg2_product * mu_ui_i + 2 * trace(wg2_block * sigma_ui_i * wg2_block'))^(-1) * si;
    end

    % Calculate E1, E2, E3
    E1 = y(:, i)' * y(:, i) - 2 * y(:, i)' * wg1_block * mu_ui_i + mu_ui_i' * wg1_product * mu_ui_i + trace(wg1_product * sigma_ui_i);
    E2 = z' * z - 2 * alphaest * z' * wg2_block * mu_ui_i + alphaest^2 * (mu_ui_i' * wg2_product * mu_ui_i) + trace((alphaest^2) * wg2_product * sigma_ui_i);
    E3 = mu_ui_i' * mu_ui_i + trace(sigma_ui_i);

    % Store results
    results(i) = struct('E1', E1, 'E2', E2, 'E3', E3, 'alphaest', alphaest);
end

% Compute final estimates
sigma1est = 1 / (m1 * n1) * sum([results.E1]);
sigma2est = 1 / (m1 * n2) * sum([results.E2]);
sigmauest = 1 / ( p1) * sum([results.E3]);


log_likelihood_y = -0.5 * n1 * m1 * log(2 * pi * sigma1est) - 0.5 * sum([results.E1]) / sigma1est;
log_likelihood_z = -0.5 * n2 * log(2 * pi * sigma2est) - 0.5 * sum([results.E2]) / sigma2est;
log_likelihood_u = -0.5 * p1 * m1 * log(2 * pi * sigmauest) - 0.5 * sum([results.E3]) / sigmauest;
log_likelihood = log_likelihood_y + log_likelihood_z + log_likelihood_u;
all_E2_values = arrayfun(@(x) x.E2, results);
E2_null = z' * z;
 output = [sigma1est, sigma2est, sigmauest, [results.alphaest],all_E2_values' ,E2_null,log_likelihood];
end
