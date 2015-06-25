% The Gaussian splitting method in `a versatile Gaussian splitting approach
% to nonlinear state estimation and its application to noise-robust ASR'

function [states, covariances, weights] = gaussian_split3 (state, covariance, num, const, eigenvector)

K = (num-1)/2;      % Total number of components = 2K+1

dim = length(state);

max_V = eigenvector;
eta = const*sqrt(3);
norm_max_V = max_V/sqrt(max_V'/covariance*max_V);

c = norm_max_V'/covariance*norm_max_V;
exp_term = exp(-0.5*(1:K).^2*eta^2*c);
first_term = 6*((sum((1:K).^2.*exp_term)^2)/(sum((1:K).^4.*exp_term)));
second_term = 2*sum(exp_term);
gamma = first_term - second_term;

states = zeros(dim, 2*K+1);
weights = zeros(1, 2*K+1);

states(:,1) = state;
weights(:,1) = gamma;

for k = 2:K+1
    states(:,k) = state + eta*(k-1)*norm_max_V;
    weights(:,k) = exp(-0.5*(k-1)^2*eta^2*c);
end
for k = K+2:2*K+1
    states(:,k) = state - eta*(k-K-1)*norm_max_V;
    weights(:,k) = weights(:,k-K);
end

weights = weights/sum(weights);
beta = 2*eta^2*sum(weights(2:K+1).*(1:K).^2);

new_covariance = covariance - beta*(norm_max_V*norm_max_V');
covariances = repmat(new_covariance, [1, 1, 2*K+1]);
