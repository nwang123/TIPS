% This function implements the Expectation-Maximization (EM) algorithm for estimating parameters in a two-group genetic data model. The function iteratively updates the parameters (sigma1, sigma2, sigmau, and alpha_g) using results from the Maximization step (`M_step_comm`), until convergence criteria are met (based on tolerance `tol` or maximum number of iterations `max_iter`).
%
% Inputs:
% - init: Initial parameter estimates [sigma1, sigma2, sigmau, alpha_g].
% - w1, w2: Genotype matrices for group 1 and group 2, respectively.
% - y: Gene expression data for group 1.
% - z: Phenotypic data for group 2.
% - n1, n2: Number of samples in groups 1 and 2, respectively.
% - max_iter: Maximum number of iterations (default: 10).
% - tol: Convergence tolerance (default: 0.001).
%
% Outputs:
% - final_parameters: Final estimates of sigma1, sigma2, sigmau, alpha_g, and log-likelihood values for both the alternative (E2) and null (E2null) models.
%
% The function iteratively calls the M-step (`M_step_comm`) to update parameter estimates. It checks for convergence by comparing the difference between current and previous parameter estimates. 
% When the difference is less than the specified tolerance, the iteration stops and the final estimates are returned.
function final_parameters = EM_comm(init, w1, y, w2, z, n1, n2, max_iter, tol)
    theta_old = init - 0.0001;
    prev_diff = Inf;
    
    for i = 1:max_iter
        step_results = M_step_comm(theta_old, w1, y, w2, z, n1, n2);
        theta_new = step_results(1:4);
        like = step_results(5:6);

        current_diff = mean(abs(theta_new - theta_old));
        if abs(current_diff - prev_diff) < tol
            break;
        end
        theta_old = theta_new;
        prev_diff = current_diff;
        fprintf('%f\n', prev_diff);
    end
    final_parameters = [theta_new, like] ;
end
