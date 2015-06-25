% Splies a Gaussian component into 2 or 3 components

function [means, covs, weights] = gaussian_split (mean, covariance, ori_weight, n)

[V, D] = eig(covariance);

[~,Dx] = find(D == max(max(D)));

max_D = D(Dx, Dx);
max_V = V(:,Dx);

v = 0.5;

means = zeros([length(mean), n]);
covs = zeros([size(covariance), n]);
weights = zeros(1, n);

if(length(max_D)>1)
    max_V = max_V(:,1);
    max_D = max_D(1);
end

means(:,1) = mean + v*sqrt(max_D)*max_V;
means(:,2) = mean - v*sqrt(max_D)*max_V;

if (n == 2)
    
    weights(1) = 0.5*ori_weight;
    weights(2) = 0.5*ori_weight;
    covs(:,:,1) = covariance - v^2*max_D*(max_V*max_V');
    covs(:,:,2) = covs(:,:,1);
    
elseif (n == 3)
    
    weights(1) = 1/6*ori_weight;
    weights(2) = 1/6*ori_weight;
    weights(3) = 4/6*ori_weight;
    means(:,3) = mean;
    covs(:,:,1) = covariance - (1/3)*v^2*max_D*(max_V*max_V');
    covs(:,:,2) = covs(:,:,1);
    covs(:,:,3) = covs(:,:,1);
    
else
    
    error('Invalid number of Gaussian mixtures')
    
end