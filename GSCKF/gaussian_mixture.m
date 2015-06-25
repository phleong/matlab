% Gaussian mixture of several Gaussian components

function [state, covariance] = gaussian_mixture (weights, state_m, cov_m)

global n_x

if(sum(weights) == 0)
    state = state_m(:,1);
    covariance = cov_m(:,:,1);
else
    weights = weights/sum(weights);  % Make sure the weights are normalised
    state = sum(repmat(weights, n_x, 1).*state_m, 2);
    cov = zeros(n_x, n_x);
    for g = 1:length(weights);
        cov = cov + weights(g)*(cov_m(:,:,g) + (state_m(:,g) - state)*(state_m(:,g) - state)');
    end
    covariance = cov;
end



