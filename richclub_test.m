function [pvalue, deg_peak, peak, peak_rand, richclub, richclub_rand_mean, richclub_normdiff, deg_uni] = richclub_test(x, x_rand)

% Reference:
% Muscoloni, A. and Cannistraci, C.V. (2017)
% "Rich-clubness test: how to determine whether a complex network
% has or doesn't have a rich-club?". arXiv:1704.03526

% Authors:
% The code has been implemented by Alessandro Muscoloni (2017).

% Released under MIT License
% Copyright (c) 2017 A. Muscoloni, C. V. Cannistraci

% NB: the code contains parallel computation. 

% It performs the statistical test for rich-clubness giving in output a pvalue
% that indicates whether the network contains a significant rich-club,
% and a degree-cut that allows to extract the rich-club subnetwork.

%%% INPUT %%%
% x - adjacency matrix of the network, which must be symmetric, zero-diagonal, unweighted.
%
% x_rand - cell array whose elements are the adjacency matrices of the randomized networks.
%          NB: it is recommended to input at least 1000 randomized networks.

%%% OUTPUT %%%
% pvalue - pvalue computed considering the peak of the normalized richclub coefficient
%          for the input network with respect to the empirical distribution
%          of the peaks in the randomized networks.
% deg_peak - degree corresponding to the peak;
%          if the pvalue is significant, the richclub subnetwork is
%          composed of the nodes with degree greater than deg_peak.
% peak - maximum value of the normalized richclub coefficient in the input network.
% peak_rand - vector containing the peak for each randomized network.
% richclub - vector containing for each degree value in the input network
%          the richclub coefficient (the last degree value is not considered,
%          since there are not nodes with degree greater than the maximum);
%          if the last k-subnetwork is composed of one node, the richclub
%          coefficient will be NaN.
% richclub_rand_mean - vector containing for each degree value in the input network
%          the mean richclub coefficient of the randomized networks.
% richclub_normdiff - vector containing for each degree value in the input network
%          the normalized richclub coefficient (the last degree value is not considered).
% deg_uni - vector of unique degree values in increasing order.

% check input
narginchk(2,2)
validateattributes(x, {'numeric'}, {'square','binary'});
if ~issymmetric(x)
    error('The input matrix must be symmetric.')
end
if any(x(speye(size(x))==1))
    error('The input matrix must be zero-diagonal.')
end
validateattributes(x_rand, {'cell'}, {'vector'});

% initialization
deg = full(sum(x,1));
deg_uni = unique(deg)';
d = length(deg_uni);
m = length(x_rand);
richclub = zeros(d-1,1);
richclub_rand_mean = zeros(d-1,1);
richclub_normdiff = zeros(d-1,1);
richclub_rand_normdiff = zeros(d-1,m);

% for each degree value
for j = 1:d-1
    
    k = deg_uni(j);
    
    % for each randomized network compute the richclub coefficient
    % as the density of the subnetwork of nodes with degree greater than k
    richclub_rand_k = zeros(m,1);
    parfor i = 1:m
        x_ik = x_rand{i}(deg>k,deg>k);
        richclub_rand_k(i) = mean(x_ik(triu(true(size(x_ik)),1)));
    end
    
    % normalization factor: mean richclub coefficient of the randomized networks
    richclub_rand_mean(j) = mean(richclub_rand_k);
    
    % for each randomized network compute the normalized richclub coefficient
    richclub_rand_normdiff(j,:) = richclub_rand_k - richclub_rand_mean(j);
    
    % compute the richclub and normalized richclub coefficients for the input network
    x_k = x(deg>k,deg>k);
    richclub(j) = mean(x_k(triu(true(size(x_k)),1)));
    richclub_normdiff(j) = richclub(j) - richclub_rand_mean(j);
    
end

% for each randomized network compute the peak of the normalized richclub coefficient
peak_rand = max(richclub_rand_normdiff, [], 1);

% compute the peak of the normalized richclub coefficient
% and the corresponding degree for the input network
[peak, deg_peak] = max(richclub_normdiff);
deg_peak = deg_uni(deg_peak);

% compute the pvalue for the peak
pvalue = mean(peak_rand >= peak);
