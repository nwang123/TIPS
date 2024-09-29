% This function performs the Maximization (M-step) of the Expectation-Maximization (EM) algorithm.
% It estimates the variance parameters (sigma1, sigma2, sigmau) and the effect size alpha_g for genetic association analysis. The function works by solving a system of linear equations to find the estimates of the hidden variable `u` and uses these estimates to update the parameter estimates.
%
% Inputs:
% - old: Vector of previous parameter estimates (sigma1, sigma2, sigmau, alpha_g).
% - w1, w2: Genotype matrices for group 1 and group 2, respectively.
% - y: Gene expression data for group 1.
% - z: Phenotypic data for group 2.
% - n1, n2: Number of samples in groups 1 and 2, respectively.
%
% Outputs:
% - esti: A vector containing the updated estimates for sigma1, sigma2, sigmau, alpha_g, and 
%   the sum of squared residuals for both null (E2null) and alternative models (E2).
%
% The function calculates the inverse of the weighted covariance matrix (`A`) for each iteration, uses it to compute the mean (`mu_ui_i`), and updates the parameter estimates based on the current estimates of `mu_ui_i`. A small regularization term (`delta`) is added to ensure numerical stability.
%
% Finally, the function returns the updated estimates for sigma1, sigma2, sigmau, alpha_g, and the log-likelihood contributions for the alternative and null models.
function [esti] = M_step_comm(old, w1, y, w2, z, n1, n2)
    sigma1 = old(1); 
    sigma2 = old(2); 
    sigmau = old(3); 
    alpha_g = old(4);

    A = (1 / sigma1) * (w1' * w1) + (alpha_g.^2) / sigma2 * (w2' * w2) + 1 / sigmau * eye(size(w1, 2));
    delta = 1e-4; % Small regularization term
    A = A + delta * eye(size(A));
    I = eye(size(A)); % Create an identity matrix of the same size as A
    sigma_ui_i = A \ I; % Compute the inverse
    mu_ui_i = sigma_ui_i * (1/sigma1 * w1' * y + alpha_g/sigma2 * w2' * z);
    alphaest = inv(mu_ui_i'*w2'*w2*mu_ui_i + sum(diag(w2'* w2*sigma_ui_i)))*(z'*w2 * mu_ui_i);
    E1 = y' * y - 2 * y' * w1 * mu_ui_i + mu_ui_i' * w1' * w1 * mu_ui_i + sum(diag(w1'*w1*sigma_ui_i));
    E2 = z' * z - 2 * alphaest * z' * w2 * mu_ui_i + alphaest^2 * mu_ui_i' * w2' * w2 * mu_ui_i+alphaest^2 *sum(diag(w1'*w1*sigma_ui_i));
    E3 = mu_ui_i' * mu_ui_i+ sum(diag(sigma_ui_i));
    sigma1est = 1/n1 * E1 ;
    sigma2est = 1/n2 * E2 ;
    sigmauest = 1/5 * E3 ;
    E2null = z' * z
    
    esti = [sigma1est, sigma2est, sigmauest, alphaest,E2,E2null];
end
