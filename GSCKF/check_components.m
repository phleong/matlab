% When one of the components have negligible weight, it is discarded and
% one of the components will be split into 2 components

function [new_states, new_covs, new_weights] = check_components (old_states, old_covs, old_weights, nonlinearities)

new_states = zeros(size(old_states));
new_covs = zeros(size(old_covs));
new_weights = zeros(size(old_weights));

% The threshold for discarding a component
threshold = 1e-4;

[sorted_weights, ind] = sort(old_weights, 'descend');

if sorted_weights(end) < threshold
    
    [~, ind_nonlinearity]= sort(nonlinearities(ind(1:end-1)), 'descend');
    states = old_states(:,ind(1:end-1));
    covs = old_covs(:,:,ind(1:end-1));
    weights = old_weights(:,ind(1:end-1));
    [new_states(:,1:2), new_covs(:,:,1:2), new_weights(:,1:2)] = gaussian_split (states(:,ind_nonlinearity(1)), covs(:,:,ind_nonlinearity(1)), weights(:,ind_nonlinearity(1)), 2);
    new_states(:,3:end) = states(:,ind_nonlinearity(2:end));
    new_covs(:,:,3:end) = covs(:,:,ind_nonlinearity(2:end));
    new_weights(:,3:end) = weights(:,ind_nonlinearity(2:end));
    
    new_weights = new_weights/sum(new_weights);
    
else
    
    new_states = old_states;
    new_covs = old_covs;
    new_weights = old_weights;
    
end