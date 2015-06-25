% Gaussian Sum Cubature Kalman Filter
% -- measure of nonlinearity -- Li
% -- method of merging components -- no interactions between CKFs

function [target_est, target_cov, nees] = GSCKF (ownship, measurement, target, target_speed_bar, const, th)

global n_x N_F cubature_points covariance_matrix transition_matrix sigma_theta 

% Outputs - estimated target state and covariance
target_est = zeros(n_x, length(measurement));
target_cov = zeros(n_x, n_x, length(measurement));
num_components = ones(1, length(measurement))*N_F;
nonlinearities = zeros(N_F, length(measurement)-1);

% Storing values
states_NF_pred = zeros(n_x, N_F, length(measurement));        % Predicted States
cov_NF_pred = zeros(n_x, n_x, N_F, length(measurement));      % Predicted Covariance
states_NF = zeros(n_x, N_F, length(measurement));             % Updated States
cov_NF = zeros(n_x, n_x, N_F, length(measurement));           % Updated Covariance
weights_NF = zeros(1, N_F, length(measurement));              % Updated Weights

n = 3;      % Number of gaussians in the mixture

% Initialisation
for m = 1:N_F;
    [states_NF(:,m,1), cov_NF(:,:,m,1)] = rp_filter_initialization (ownship(:,1), measurement(1), m, target_speed_bar);
    weights_NF(:,m,1) = 1/N_F;
end

[target_est(:,1), target_cov(:,:,1)] = gaussian_mixture (weights_NF(:,:,1), states_NF(:,:,1), cov_NF(:,:,:,1));
% num_components(1) = N_F;

% For each time step
for k = 2:length(measurement)
    
    input_vector = ownship(:,k) - transition_matrix*ownship(:,k-1);
    
    for m = 1:N_F

        % Prediction
        states_NF_pred(:,m,k) = transition_matrix*states_NF(:,m,k-1) - input_vector;
        cov_NF_pred(:,:,m,k) = transition_matrix*cov_NF(:,:,m,k-1)*transition_matrix' + covariance_matrix;
        
        % Update
        points_pred = chol(cov_NF_pred(:,:,m,k), 'lower')*cubature_points + repmat(states_NF_pred(:,m,k),1,2*n_x);
        
        points_meas_pred = atan2(points_pred(1,:), points_pred(2,:));
        
        nonlinearities(m,k-1) = nonlinearity_degree (points_pred, states_NF_pred(:,m,k), points_meas_pred, bearings_average(points_meas_pred, ones(1,2*n_x)/(2*n_x)), cov_NF_pred(:,:,m,k), ones(1,2*n_x)/(2*n_x));

        if (nonlinearities(m,k-1) > th)
            
            updated_state = zeros(n_x,n);
            updated_cov = zeros(n_x,n_x,n);
            updated_weight = zeros(1,n);
            num_components(k) = num_components(k) + n - 1;
          
            vector = dir_nonlinearity (points_meas_pred,atan2(states_NF_pred(1,m,k), states_NF_pred(2,m,k)), points_pred,states_NF_pred(:,m,k));
            [predicted_state1, state_prediction_cov1, weights1] = gaussian_split3 (states_NF_pred(:,m,k), cov_NF_pred(:,:,m,k),n,const,vector);
            
            for a = 1:n

                points_pred = chol(state_prediction_cov1(:,:,a), 'lower')*cubature_points + repmat(predicted_state1(:,a), 1, 2*n_x);
                points_meas_pred = atan2(points_pred(1,:), points_pred(2,:));
                pred_measurement = bearings_average(points_meas_pred, ones(1, 2*n_x)/(2*n_x));

                difference = wraparound(points_meas_pred - pred_measurement*ones(1,2*n_x));
                innov_covariance = difference*difference'/(2*n_x) + sigma_theta^2;
                cross_covariance = (points_pred - predicted_state1(:,a)*ones(1,2*n_x))*difference'/(2*n_x);

                kalman_gain = cross_covariance/innov_covariance;
                innovation = wraparound(measurement(k) - pred_measurement);
                updated_state(:,a) = predicted_state1(:,a) + kalman_gain*innovation;
                updated_cov(:,:,a) = state_prediction_cov1(:,:,a) - kalman_gain*innov_covariance*kalman_gain';

                likelihood = exp(-0.5*(innovation^2)/innov_covariance)/sqrt(2*pi*innov_covariance);
                updated_weight(a) = likelihood*weights_NF(:,m,k-1)*weights1(a);
                
            end
       
            weights_NF(:,m,k) = sum(updated_weight);
            [states_NF(:,m,k), cov_NF(:,:,m,k)] = gaussian_mixture(updated_weight, updated_state, updated_cov);
                
        else
            
            pred_measurement = bearings_average(points_meas_pred, ones(1, 2*n_x)/(2*n_x));
        
            difference = wraparound(points_meas_pred - pred_measurement*ones(1,2*n_x));
            innov_covariance = difference*difference'/(2*n_x) + sigma_theta^2;
            cross_covariance = (points_pred - repmat(states_NF_pred(:,m,k),1,2*n_x))*difference'/(2*n_x);
            kalman_gain = cross_covariance/innov_covariance;
        
            innovation = wraparound(measurement(k) - pred_measurement);
            states_NF(:,m,k) = states_NF_pred(:,m,k) + kalman_gain*innovation;
            cov_NF(:,:,m,k) = cov_NF_pred(:,:,m,k) - kalman_gain*innov_covariance*kalman_gain';
        
            likelihood = exp(-0.5*(innovation^2)/innov_covariance)/sqrt(2*pi*innov_covariance);
            weights_NF(:,m,k) = likelihood*weights_NF(:,m,k-1);

        end
        
    end

    weights_NF(:,:,k) = weights_NF(:,:,k)/sum(weights_NF(:,:,k));
    
    % Combining outputs
    [target_est(:,k), target_cov(:,:,k)] = gaussian_mixture(weights_NF(:,:,k), states_NF(:,:,k), cov_NF(:,:,:,k));
    
    % Make sure there is no component with negligible weights
    [states_NF(:,:,k), cov_NF(:,:,:,k), weights_NF(:,:,k)] = check_components(states_NF(:,:,k), cov_NF(:,:,:,k), weights_NF(:,:,k), nonlinearities(:,k-1));

end

target_est = target_est + ownship;

% NEES
nees = zeros(1, length(measurement));
for k = 1:length(measurement)
    nees(k) = (target(:,k) - target_est(:,k))'/target_cov(:,:,k)*(target(:,k) - target_est(:,k));
end
