% Cramer-Rao Lower Bound for RMS error

function filter_bound = CRLB (target_trajectories, ownship, init_covs)

global transition_matrix covariance_matrix sigma_theta n_x 

num_runs = size(init_covs,3);
all_steps = size(target_trajectories,2);

filter_bound = zeros(1, all_steps);

K_11 = transition_matrix'/covariance_matrix*transition_matrix;
K_12 = -transition_matrix'/covariance_matrix;
K_22 = inv(covariance_matrix);

for run = 1:num_runs
    
    filter_FIM = zeros(n_x, n_x, all_steps);
    pred_FIM = zeros(n_x, n_x, all_steps);
    relative_states_run = target_trajectories(:,:,run) - ownship;
    H = [relative_states_run(2,:)./(relative_states_run(1,:).^2 + relative_states_run(2,:).^2); 
        -relative_states_run(1,:)./(relative_states_run(1,:).^2 + relative_states_run(2,:).^2);
        zeros(2, all_steps)];
    
    for k = 1:all_steps
        
        if k == 1
            pred_FIM(:,:,k) = K_22;
            filter_FIM(:,:,k) = inv(init_covs(:,:,:,run));
        else
            pred_FIM(:,:,k) = K_22 - K_12'/(K_11 + filter_FIM(:,:,k-1))*K_12;
            filter_FIM(:,:,k) = pred_FIM(:,:,k) + H(:,k)/(sigma_theta^2)*H(:,k)';
        end     
        
        inverse = inv(filter_FIM(:,:,k));
        filter_bound(k) = filter_bound(k) + sqrt(inverse(1,1) + inverse(2,2))/num_runs;
 
    end  
    
end
