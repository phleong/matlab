% This function calculates the average of a list of bearing values, with their corresponding weights

function average = bearings_average(bearings, weights)

% To make sure weights are normalized
weights = weights/sum(weights);

average = bearings(1);
for k = 2:length(weights)
    diff = wraparound(bearings(k) - bearings(1));
    average = average + diff*weights(k);
end

average = wraparound(average);