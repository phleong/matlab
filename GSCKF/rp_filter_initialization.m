% Initializaing the filtering algorithm with a number of independent
% filters, each with a different range estimate

function [init_estimate, init_state_cov] = rp_filter_initialization (init_ownship, init_bearing_bar, n, target_speed_bar)

global sigma_theta sigma_s sigma_c r_min common_ratio

% Prior knowledge of target trajectory
init_range_bar = r_min*(common_ratio^(n-1) + common_ratio^n)/2;
target_course_bar = init_bearing_bar + pi;               % in radians

sigma_r_n = r_min*(common_ratio^n - common_ratio^(n-1))/sqrt(12);

% Initialization of state vector and its covariance
init_estimate = [init_range_bar*sin(init_bearing_bar); 
                 init_range_bar*cos(init_bearing_bar); 
                 target_speed_bar*sin(target_course_bar) - init_ownship(3); 
                 target_speed_bar*cos(target_course_bar) - init_ownship(4)];

P_xx = (init_range_bar*sigma_theta*cos(init_bearing_bar))^2 + (sigma_r_n*sin(init_bearing_bar))^2;
P_yy = (init_range_bar*sigma_theta*sin(init_bearing_bar))^2 + (sigma_r_n*cos(init_bearing_bar))^2;
P_xy = (sigma_r_n^2-(init_range_bar*sigma_theta)^2)*sin(init_bearing_bar)*cos(init_bearing_bar);
P_xx_dot = (target_speed_bar*sigma_c*cos(target_course_bar))^2 + (sigma_s*sin(target_course_bar))^2;
P_yy_dot = (target_speed_bar*sigma_c*sin(target_course_bar))^2 + (sigma_s*cos(target_course_bar))^2;
P_xy_dot = (sigma_s^2 - (target_speed_bar*sigma_c)^2)*sin(target_course_bar)*cos(target_course_bar);

init_state_cov = [P_xx, P_xy, 0, 0; P_xy, P_yy, 0, 0; 0, 0, P_xx_dot, P_xy_dot; 0, 0, P_xy_dot, P_yy_dot];
