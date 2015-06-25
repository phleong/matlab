% Initialization of the filter for bearings-only tracking

function [init_estimate, init_state_cov] = filter_initialization (init_ownship, init_bearing_bar, target_speed_bar, target_range_bar)

global sigma_r sigma_theta sigma_s sigma_c

% Prior knowledge of target trajectory
target_course_bar = init_bearing_bar + pi;               % in radians

% Initialization of state vector and its covariance
init_estimate = [target_range_bar*sin(init_bearing_bar); 
                 target_range_bar*cos(init_bearing_bar); 
                 target_speed_bar*sin(target_course_bar) - init_ownship(3); 
                 target_speed_bar*cos(target_course_bar) - init_ownship(4)];

P_xx = (target_range_bar*sigma_theta*cos(init_bearing_bar))^2 + (sigma_r*sin(init_bearing_bar))^2;
P_yy = (target_range_bar*sigma_theta*sin(init_bearing_bar))^2 + (sigma_r*cos(init_bearing_bar))^2;
P_xy = (sigma_r^2-(target_range_bar*sigma_theta)^2)*sin(init_bearing_bar)*cos(init_bearing_bar);
P_xx_dot = (target_speed_bar*sigma_c*cos(target_course_bar))^2 + (sigma_s*sin(target_course_bar))^2;
P_yy_dot = (target_speed_bar*sigma_c*sin(target_course_bar))^2 + (sigma_s*cos(target_course_bar))^2;
P_xy_dot = (sigma_s^2 - (target_speed_bar*sigma_c)^2)*sin(target_course_bar)*cos(target_course_bar);

init_state_cov = [P_xx, P_xy, 0, 0; P_xy, P_yy, 0, 0; 0, 0, P_xx_dot, P_xy_dot; 0, 0, P_xy_dot, P_yy_dot];
