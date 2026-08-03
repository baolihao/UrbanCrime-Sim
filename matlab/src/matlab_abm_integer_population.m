function population = matlab_abm_integer_population(expected)
%MATLAB_ABM_INTEGER_POPULATION Randomly round a density with a fixed total.
% The total equals ceil(sum(expected(:))). Fractional agents are sampled
% without replacement with weights proportional to their fractional parts.

validateattributes(expected, {'numeric'}, {'real', 'finite', 'nonnegative'});
population = floor(expected);
remainder = ceil(sum(expected, 'all')) - sum(population, 'all');
if remainder == 0
    return
end
fractions = expected(:) - floor(expected(:));
eligible = find(fractions > 0);
if remainder > numel(eligible)
    error('urbancrime:PopulationInvariant', ...
        'Fractional population cannot supply the requested integer total.');
end
weights = fractions(eligible);
chosen = zeros(remainder, 1);
for index = 1:remainder
    threshold = rand() * sum(weights);
    local_index = find(cumsum(weights) >= threshold, 1, 'first');
    chosen(index) = eligible(local_index);
    weights(local_index) = 0;
end
population(chosen) = population(chosen) + 1;
end
