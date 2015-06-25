% Determines the direction of nonlinearity

function vector = dir_nonlinearity (points, center, points_x, center_x)

dim = length(points)/2;

matrix = zeros(dim, dim);
for n = 1:dim;
    dist = 0.5*(wraparound(points(n) + points(n+dim) - 2*center))^2;
    phi = (points_x(:,n) - center_x)/norm(points_x(:,n) - center_x);
    matrix = matrix + dist*(phi*phi');
end

[V, D] = eig(matrix);
[~,Dx] = find(D == max(max(D)));

vector = V(:,Dx);
