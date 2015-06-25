% Cubature Kalman Filter for the bearings-only tracking problem

function [target_est, target_cov, nees] = CKF (ownship, measurement, target, target_range_bar, target_speed_bar)

global n_x cubature_points covariance_matrix transition_matrix sigma_theta

% Outputs - estimated target state and covariance
target_est = zeros(n_x, length(measurement));
target_cov = zeros(n_x, n_x, length(measurement));

% Storing values
states_pred = zeros(n_x, length(measurement));          % Predicted States
cov_pred = zeros(n_x, n_x, length(measurement));        % Predicted Covariance

% Initialisation
[target_est(:,1), target_cov(:,:,1)] = filter_initialization (ownship(:,1), measurement(1), target_speed_bar, target_range_bar);

% For each time step
for k = 2:length(measurement)
    
    input_vector = ownship(:,k) - transition_matrix*ownship(:,k-1);
    
    % Prediction
    states_pred(:,k) = transition_matrix*target_est(:,k-1) - input_vector;
    cov_pred(:,:,k) = transition_matrix*target_cov(:,:,k-1)*transition_matrix' + covariance_matrix;
        
    % Update
    points_pred = chol(cov_pred(:,:,k), 'lower')*cubature_points + repmat(states_pred(:,k), 1, 2*n_x);
    points_meas_pred = atan2(points_pred(1,:), points_pred(2,:));
    pred_measurement = bearings_average(points_meas_pred, ones(1, 2*n_x)/(2*n_x));
        
    difference = wraparound(points_meas_pred - pred_measurement*ones(1,2*n_x));
    innov_covariance = difference*difference'/(2*n_x) + sigma_theta^2;
    cross_covariance = (points_pred - repmat(states_pred(:,k),1,2*n_x))*difference'/(2*n_x);
    kalman_gain = cross_covariance/innov_covariance;
        
    innovation = wraparound(measurement(k) - pred_measurement);
    target_est(:,k) = states_pred(:,k) + kalman_gain*innovation;
    target_cov(:,:,k) = cov_pred(:,:,k) - kalman_gain*innov_covariance*kalman_gain';
        
end

target_est = target_est + ownship;

% NEES
nees = zeros(1, length(measurement));
for k = 1:length(measurement)
    nees(k) = (target(:,k) - target_est(:,k))'/target_cov(:,:,k)*(target(:,k) - target_est(:,k));
end
