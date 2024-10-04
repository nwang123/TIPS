%  MATLAB Code: Gene Expression Data Simulation and Estimation

%  This MATLAB script performs the simulation of genetic data, specifically working with SNP (Single Nucleotide Polymorphism) data and gene expression. The code applies an EM (Expectation-Maximization) algorithm to estimate effect sizes (`alpha`) and variances (`sigma1`, `sigma2`, `sigmau`) using both genotypic and phenotypic data.

%  Key Variables:
% **n1, n2**: Number of samples in the first group (e.g., 205) and second group (e.g., 16,470).
% **m1**: Number of genes (e.g., 440).
% **p1**: Number of SNPs (e.g., 6,769).
% **alpha_initial**: Initial effect size estimates for each gene.
% **k**: SNP loadings, indicating the relative influence of different SNPs on gene expression.
% **wg1_mat, wg2_mat**: Genotypic matrices for two groups, `wg1_mat` for group 1 and `wg2_mat` for group 2.
% **Y_mat**: Gene expression matrix for group 1.
% **newz_mat**: Standardized phenotype data for group 2.

% The script begins by initializing key variables, such as the sample sizes `n1`, `n2`, the number of genes `m1`, the number of SNPs `p1`, and the initial effect sizes `alpha_initial`. SNP loadings `k` are predefined for different genes, representing the number of SNPs that contribute to each gene's expression.

% n1 = 205; n2 = 16470  ; m1 = 440; p1  = 6769; 
% alpha_initial = normrnd(0, 0.01, [m1, 1]);
% old = [0.01,0.01,0.01, alpha_initial'];
% The following is the Pathway grouping information.
% k = [repmat(7, 1, 7), repmat(8, 1, 8), repmat(6, 1, 12), repmat(5, 1, 5), ...
         repmat(6, 1, 6), repmat(5, 1, 5), repmat(9, 1, 9), repmat(6, 1, 6), ...
         repmat(8, 1, 8), repmat(5, 1, 5), repmat(6, 1, 6), repmat(10, 1, 10), ...
         repmat(12, 1, 12), repmat(7, 1, 7), repmat(6, 1, 6), repmat(7, 1, 7), ...
         repmat(12, 1, 12), repmat(5, 1, 5),repmat(17, 1, 17),repmat(5, 1, 10),...
         repmat(14, 1, 14), repmat(19, 1, 19), repmat(8, 1, 8), repmat(5, 1, 5),...
         repmat(7, 1, 7), repmat(18, 1, 18), repmat(5, 1, 5), repmat(6, 1, 6),...
         repmat(5, 1, 5), repmat(7, 1, 14), repmat(20, 1, 20), repmat(12, 1, 12),...
         repmat(6, 1, 6), repmat(18, 1, 18), repmat(1, 1, 120)];
% Y_mat = readtable("y_gene_brain_chu1.csv");
% Y_mat = table2array(Y_mat);
% wg1_mat = readtable("w1_brain_chu1.csv");
% wg1_mat = table2array(wg1_mat);
% Perform the matrix operations
% u = pinv(wg1_mat' * wg1_mat) * wg1_mat' * Y_mat;
% wg2_mat = readtable("w2_brain_chu1.csv");
% wg2_mat = table2array(wg2_mat);

% X = wg2_mat * u;
% newz_mat = load('z_brain.mat');
% newz_mat  = newz_mat.w1_brain;
% newz_mat = (newz_mat - mean(newz_mat)) ./ std(newz_mat);
% Compute the OLS estimate for alpha
% alpha1 = pinv(X' * X) * X' * newz_mat; 

% Estimating residuals
% e1 = Y_mat - wg1_mat * u;
% e2 = newz_mat - X * alpha1;

% Estimating sigmas
% sigma1 = std(e1);
% sigma2 = std(e2);
% sigmau = std(u);
% old = [mean(sigma1), mean(sigma2), mean(sigmau), alpha1'];

% a = 0.5;
% lam = 75; % Set lambda to 0
         
% results = EM_updated(old, a, wg1_mat, wg2_mat, newz_mat, Y_mat, lam, k, m1, n1, n2, p1);
