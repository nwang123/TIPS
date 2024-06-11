function [final_parameters] = EM_updated(old,  a_para, wg1, wg2, z, y, lambda, k, m1, n1, n2, p1, max_iter, tol) % [final_parameters, history]
    if nargin < 13
        max_iter = 10; % default value if not provided
    end
   
    if nargin < 14
        tol = 0.001; % default value if not provided
    end

    % Initialize the parameters
    theta_old = old - 0.0001;
    prev_diff = Inf; % Initialize the previous difference to an arbitrarily large value
    %iter = 1;

    % Create a list to store parameter estimates for each iteration
    %theta_history = cell(1, max_iter);
%    %theta_history{iter} = theta_old;
    
    % Loop over the maximum number of iterations
    current_diff = NaN(1,max_iter);
    prev_diff = Inf(1,max_iter);
    theta_history = cell(1, max_iter);
    theta_old = cell(1,max_iter);
    theta_old{1} = old(1:(3+m1));
    for iter = 1:max_iter
        disp(['start index = ' num2str(iter) ' in iter using EM_updated function.']);
        step_results = mstep_updated(theta_old{iter}, a_para, wg1, wg2, z, y, lambda, k, m1, n1, n2, p1);
        theta_history{iter} = step_results(1:(3+m1))
        like = step_results((4+m1):end);
    % Compute the current difference
        current_diff(iter) = mean(abs(theta_history{iter} - theta_old{iter}));
    % Check for convergence using the difference between the current and previous differences
        if abs(current_diff(iter) - prev_diff(iter)) < tol
            disp(['Break the loop when index = ' num2str(iter) ' due to difference < tol using EM_updated function.']);
            final_parameters = [theta_history{iter},lambda,like];
            disp(['Finished using EM_updated function.']);
            break;
        else
            theta_old{(iter+1)} = theta_history{iter};
            prev_diff((iter+1)) = current_diff(iter);
            disp(['Finished index : ' num2str(iter) ' in iter using EM_updated function.']);
        end
    end


     % Return the estimated parameters
    %final_parameters = [theta_history{res_iter},lambda]; %, lambda];
    %history = theta_history(1:max_iter); % return only the non-empty cells
    disp(['Finished EM_updated function.']);
end




%    for iter = 1:max_iter
%        theta_new = mstep(theta_old, a_para, wg1, wg2, z, y, lambda, k, m1, n1, n2, p1);
%        theta_history{iter} = theta_new;

%        % Compute the current difference
%        current_diff = mean(abs(theta_new - theta_old));

%        % Check for convergence using the difference between the current and previous differences
%        if abs(current_diff - prev_diff) < tol
%            break;
%        end

%        % Update for the next iteration
%        theta_old = theta_new;
%        prev_diff = current_diff;
%        iter = iter + 1;
%        disp(prev_diff);
%    end

%    % Return the estimated parameters
%    final_parameters = [theta_new, lambda];
%    history = theta_history(1:iter); % return only the non-empty cells
%end
