% Bearings-only tracking using the Gaussian-sum cubature Kalman filter with
% improved robustness as described in
% P.H.Leong, S.Arulampalam, T.A.Lamahewa, and T.D.Abhayapala, "Gaussian-sum
% cubature Kalman filter with improved robustness for bearings-only
% tracking", IEEE Signal Processing Letters, vol.21, no.5, pp.513-517,
% May 2014.

clc
close all
clear all

% data generation for bearings-only tracking - target with constant
% velocity
% 'm' - moderately nonlinear scenario
% 'h' - highly nonlinear scenario

scenario = 'm';   

M = 500;                            % Number of independent runs

global N_F 
N_F = 5;                            % Number of independent CKF
eta = 0.7;                          % Displacement parameter of the Gaussian components
th = 0.1;                           % Threshold for splitting a component

% Constants
global n_x knot_to_kmps; 
n_x = 4;                            % Dimensionality of the state
knot_to_kmps = 5.14e-4;             % Conversion from knots to km/s

% Model
global transition_matrix covariance_matrix covariance_matrix_sqrt;

global sigma_theta sigma_r sigma_s sigma_c; 

if (scenario == 'm') % moderately nonlinear scenario
    
    % Simulation parameters
    num_steps = 32;                     % Number of time steps
    T = 60;                             % Sampling interval (seconds)
    transition_matrix = [1, 0, T, 0; 0, 1, 0, T; 0, 0, 1, 0; 0, 0, 0, 1];
    
    % Variables
    sigma_r = 2;                        % standard deviation of range estimate (km)
    sigma_s = 2*knot_to_kmps;           % standard deviation of target speed (km/s)
    sigma_c = pi/sqrt(12);              % standard deviation of target course (radians)
    sigma_theta = 1.5*pi/180;           % standard deviation of bearing/angle measurement (radians)
    q_tilde = 9e-12;
    
    ownship = zeros(n_x, num_steps+1);
    init_target = zeros(n_x, 1);
    % Ownship trajectory
    ownship_speed = 5*knot_to_kmps;     % in km/s
    ownship_init_course = 140*pi/180;   
    ownship_final_course = 20*pi/180;
    ownship(3:4, 1) = [ownship_speed*sin(ownship_init_course); ownship_speed*cos(ownship_init_course)];  

    % True target trajectory
    init_target(1:2, 1) = [4.92; 0.85];
    target_speed_true = 4*knot_to_kmps;
    target_course_true = -140*pi/180;
    init_target(3:4, 1) = [target_speed_true*sin(target_course_true); target_speed_true*cos(target_course_true)];

    prev_ownship_course = ownship_init_course;
    for a = 1:num_steps
        ownship(:,a+1) = transition_matrix*ownship(:,a);
        if ((a>=13) && (a<=16))
            current_ownship_course = prev_ownship_course + (ownship_final_course-ownship_init_course)/4;
            ownship(3,a+1) = ownship_speed*sin(current_ownship_course);
            ownship(4,a+1) = ownship_speed*cos(current_ownship_course);
            prev_ownship_course = current_ownship_course;
        end  
    end

elseif (scenario == 'h') % highly nonlinear scenario
    
    % Simulation parameters
    num_steps = 30;                     % Number of time steps
    T = 60;                             % Sampling interval (seconds)
    transition_matrix = [1, 0, T, 0; 0, 1, 0, T; 0, 0, 1, 0; 0, 0, 0, 1];
    
    % Variables
    sigma_r = 4;                        % standard deviation of range estimate (km)
    sigma_s = 2*knot_to_kmps;           % standard deviation of target speed (km/s)
    sigma_c = pi/sqrt(12);              % standard deviation of target course (radians)
    sigma_theta = 2*pi/180;             % standard deviation of bearing/angle measurement (radians)
    q_tilde = 9e-12;
    
    ownship = zeros(n_x, num_steps+1);
    init_target = zeros(n_x, 1);
    % Ownship trajectory
    ownship_speed = 5*knot_to_kmps;     % in km/s
    ownship_init_course = -80*pi/180;   
    ownship_final_course = 146*pi/180;
    ownship(3:4, 1) = [ownship_speed*sin(ownship_init_course); ownship_speed*cos(ownship_init_course)];  

    % True target trajectory
    init_target(1:2, 1) = [7.0909; 7.0526];
    target_speed_true = 15*knot_to_kmps;
    target_course_true = -135.4*pi/180;
    init_target(3:4, 1) = [target_speed_true*sin(target_course_true); target_speed_true*cos(target_course_true)];

    for a = 1:num_steps
        ownship(:,a+1) = transition_matrix*ownship(:,a);
        if (a == 15)
            ownship(3,a+1) = ownship_speed*sin(ownship_final_course);
            ownship(4,a+1) = ownship_speed*cos(ownship_final_course);
        end  
    end
 
else
    error('Error. Enter m or h for scenario.')
    break
end

covariance_matrix = [(T^3)/3, 0, (T^2)/2, 0; 0, (T^3)/3, 0, (T^2)/2; (T^2)/2, 0, T, 0; 0, (T^2)/2, 0, T]*q_tilde;
covariance_matrix_sqrt = sqrtm(covariance_matrix);

target_trajectories = zeros(n_x, num_steps+1, M);
measurements = zeros(M, num_steps+1);

target_range_bar = zeros(1, M);
target_speed_bar = zeros(1, M);

for run = 1:M;

    target_trajectories(:,1,run) = init_target;
    relative_state1 = init_target - ownship(:,1);
    measurements(run, 1) = atan2(relative_state1(1), relative_state1(2)) + sigma_theta*randn;
    
    % Initialization of target range and speed
    target_range_bar(:,run) = sqrt(sum(relative_state1(1:2).^2)) + sigma_r*randn; % in km   
    while (target_range_bar(:,run) <= 0.6)
        target_range_bar(:,run) = sqrt(sum(relative_state1(1:2).^2)) + sigma_r*randn;
    end
    target_speed_bar(:,run) = target_speed_true + sigma_s*randn;    % in km/s
    while (target_speed_bar(:,run) <= 0.6*knot_to_kmps)
        target_speed_bar(:,run) = target_speed_true + sigma_s*randn; 
    end
            
    for k = 1:num_steps
        target_trajectories(:,k+1,run) = transition_matrix*target_trajectories(:,k,run) + covariance_matrix_sqrt*randn(n_x,1);
        relative_state = target_trajectories(:,k+1,run) - ownship(:,k+1);
        measurements(run, k+1) = atan2(relative_state(1), relative_state(2)) + sigma_theta*randn;            
    end
    
end

figure;
plot(ownship(1,:), ownship(2,:)); hold on
plot(target_trajectories(1,:,1), target_trajectories(2,:,1), 'r-.');
xlabel('x (km)')
ylabel('y (km)')
legend('Observer', 'Target')
title('Target and Observer Trajectories')
text(ownship(1,1)+0.1, ownship(2,1), 'start');
text(target_trajectories(1,1,1), target_trajectories(2,1,1)-0.25, 'start');
axis equal

%%%%%%%%%%%%%
global common_ratio r_min
relative_state1 = target_trajectories(:,1,1) - ownship(:,1);
init_range = sqrt(sum(relative_state1(1:2).^2));
r_max = init_range + 3*sigma_r;
r_min = init_range - 3*sigma_r;
if(r_min<0)
    r_min = 0.1*init_range;
end
common_ratio = (r_max/r_min)^(1/N_F);

global cubature_points;
cubature_points = sqrt(n_x)*[eye(n_x), -eye(n_x)];

% Output of the GSCKF and its improved version
target_states_CKF_all = zeros(n_x, num_steps+1, M);
state_cov_CKF_all = zeros(n_x, n_x, num_steps+1, M);
target_states_GSCKF_all = zeros(n_x, num_steps+1, M);
state_cov_GSCKF_all = zeros(n_x, n_x, num_steps+1, M);

nees_CKF = zeros(1, num_steps+1, M);
nees_GSCKF = zeros(1, num_steps+1, M);

CKF_time = 0;
GSCKF_time = 0;

% Filtering
for run = 1:M
    
    run
   
    % Cubature Kalman Filter
    CKF_start = cputime;
    [target_states_CKF_all(:,:,run), state_cov_CKF_all(:,:,:,run), nees_CKF(:,:,run)] = CKF (ownship, measurements(run,:), target_trajectories(:,:,run), target_range_bar(:,run), target_speed_bar(:,run));
    CKF_time = CKF_time + cputime - CKF_start;
    
    % Gaussian-Sum Cubature Kalman Filter (Improved)
    GSCKF_start = cputime;
    [target_states_GSCKF_all(:,:,run), state_cov_GSCKF_all(:,:,:,run), nees_GSCKF(:,:,run)] = GSCKF (ownship, measurements(run,:), target_trajectories(:,:,run), target_speed_bar(:,run), eta, th);
    GSCKF_time = GSCKF_time + cputime - GSCKF_start;
       
end

% Cramer-Rao Lower Bound (CRLB)
CRLB_filter = CRLB (target_trajectories, ownship, state_cov_CKF_all(:,:,1,:));

% RMS error
RMS_error_CKF = sqrt(mean(sum((target_states_CKF_all(1:2,:,:) - target_trajectories(1:2,:,1:M)).^2,1),3));
RMS_error_GSCKF = sqrt(mean(sum((target_states_GSCKF_all(1:2,:,:) - target_trajectories(1:2,:,1:M)).^2,1),3));

% Plotting RMS Error together with its CRLB
figure;
plot((0:num_steps)*T/60, RMS_error_CKF, 'g+-'); hold on
plot((0:num_steps)*T/60, RMS_error_GSCKF, 'bs-'); hold on
plot((0:num_steps)*T/60, CRLB_filter, 'r-', 'Linewidth', 1.5); 
xlabel('Time (Minute)')
ylabel('RMS Position Error (km)')
legend('CKF', 'GSCKF', 'CRLB')
axis([16 30 0 4.5])
