% Authors:
% The code has been implemented by Alessandro Muscoloni (2017).

% Released under MIT License
% Copyright (c) 2017 A. Muscoloni, C. V. Cannistraci

% NB: the code contains parallel computation.

% USAGE EXAMPLE

% load example network
load('network.mat', 'x')

% generate 1000 randomized networks using the CM null-model
m = 1000;
x_rand = cell(m,1);
parfor i = 1:m
    x_rand{i} = randomize_network(x, 'CM');
end
save('null_model_CM.mat', 'x_rand')

% perform the rich-clubness test
[pvalue, deg_peak, peak, peak_rand, richclub, richclub_rand_mean, richclub_normdiff, deg_uni] = richclub_test(x, x_rand);
save('richclub_test.mat', 'pvalue', 'deg_peak', 'peak', 'peak_rand', 'richclub', 'richclub_rand_mean', 'richclub_normdiff', 'deg_uni')