% Wraparound theta so they lie between -pi and pi

function wraparound_theta = wraparound (theta)

new_theta = theta;
[xx1] = find(new_theta > pi);
new_theta(xx1) = new_theta(xx1) - 2*pi;
[xx2] = find(new_theta < -pi);
new_theta(xx2) = new_theta(xx2) + 2*pi;

wraparound_theta = new_theta;
