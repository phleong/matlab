% Calculates the degree of nonlinearity of the CKF using 'Measure of nonlinearity for
% stochastic systems' by X.R.Li

function nonlinearity = nonlinearity_degree (states, state_center, transformed_states, transformed_center, cov_state, weights)

diff = wraparound(transformed_states - transformed_center);
C_g = weights.*diff*diff';
C_gx = weights.*diff*(states - state_center*ones(1,length(transformed_states)))';

nonlinearity = sqrt(1 - ((C_gx/cov_state*C_gx')/C_g));
