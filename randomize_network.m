function [x, eff] = randomize_network(x, null_model, iters, max_rej)

% References:
% 1) Cannistraci-Muscoloni null-model
%   Muscoloni, A. and Cannistraci, C.V. (2017)
%   "Rich-clubness test: how to determine whether a complex network
%    has or doesn't have a rich-club?". arXiv:1704.03526
% 2) Maslov-Sneppen null-model
%   Maslov, S. and Sneppen, K. (2002)
%   "Specificity and Stability in Topology of Protein Networks". Science 296:910

% Authors:
% The code has been implemented by Alessandro Muscoloni (2017).

% Released under MIT License
% Copyright (c) 2017 A. Muscoloni, C. V. Cannistraci

% It performs a randomization of the network in input preserving the node degrees.
% The randomization consists of an iterative random reshuffle of link pairs.
% At each iteration two links are randomly sampled (uniformly or nonuniformly,
% according to the Maslov-Sneppen or Cannistraci-Muscoloni null-model respectively)
% and one endpoint of the first is randomly exchanged with one endpoint of the second.
% If the link that would be created already exists, the attempt is rejected.

%%% INPUT %%%
% x - adjacency matrix of the network, which must be symmetric, zero-diagonal, unweighted.
%
% null_model - null-model used for the randomization of the network, the alternatives are:
%      'CM' -> (Cannistraci-Muscoloni) links sampled nonuniformly;
%              one link according to probabilities proportional to the degree-product
%              and one link with inverse probabilities.
%      'MS' -> (Maslov-Sneppen) links sampled uniformly.
%
% iters - number of iterations of the randomization procedure;
%         [optional] if not given or empty, it is set by default to 10*edges.
%
% max_rej - number of consecutive rejections that stops the procedure;
%           [optional] if not given, it is set by default to Inf,
%           to indicate that the procedure will not stop due to rejections.

%%% OUTPUT %%%
% x - adjacency matrix of the randomized network.
%
% eff - number of effective iterations, in case the procedure is stopped
%       due to a maximum number of consecutive rejections.

% Possible usage:
% randomize_network(x, null_model)                  [default iters and max_rej]
% randomize_network(x, null_model, iters)           [default max_rej]
% randomize_network(x, null_model, [], max_rej)     [default iters]
% randomize_network(x, null_model, iters, max_rej)

% check input
narginchk(2,4)
validateattributes(x, {'numeric'}, {'square','binary'});
if ~issymmetric(x)
    error('The input matrix must be symmetric.')
end
if any(x(speye(size(x))==1))
    error('The input matrix must be zero-diagonal.')
end
if strcmp(null_model,'CM')
    nonuniform = 1;
elseif strcmp(null_model,'MS')
    nonuniform = 0;
else
    error('Possible null-models: ''CM'',''MS''.');
end
if ~exist('iters', 'var')
    iters = [];
elseif ~isempty(iters)
    validateattributes(iters, {'numeric'}, {'scalar','integer','positive','finite'});
end
if ~exist('max_rej', 'var')
    max_rej = Inf;
elseif ~isinf(max_rej)
    validateattributes(max_rej, {'numeric'}, {'scalar','integer','positive','finite'});
end

% initialization
n = size(x,1);
[i,j] = find(triu(x,1));
E = length(i);
eff = 0;
rej = 0;
if isempty(iters)
    iters = 10*E;
end

if nonuniform
    % compute weights
    deg = full(sum(x,1));
    degprod = repmat(deg,n,1) .* repmat(deg',1,n);
    degprod(speye(size(degprod))==1) = 0;
    degprod_rev = abs(degprod - min(degprod(triu(true(size(degprod)),1))) - max(degprod(triu(true(size(degprod)),1))));
    degprod_rev(speye(size(degprod_rev))==1) = 0;
    w1 = degprod(sub2ind(size(degprod),i,j));
    w2 = degprod_rev(sub2ind(size(degprod_rev),i,j));
end

% randomization
while eff < iters
    
    % stop if reached maximum number of consecutive rejections
    if rej==max_rej
        break;
    end
    
    % random sample two different links
    if nonuniform
        % nonuniform probabilities
        e1 = randsample(E,1,1,w1);
        e2 = randsample(E,1,1,w2);
        while e2==e1
            e2 = randsample(E,1,1,w2);
        end
    else
        % uniform probabilities
        e1 = randi(E);
        e2 = randi(E);
        while e2==e1
            e2 = randi(E);
        end
    end
    a=i(e1); b=j(e1);
    c=i(e2); d=j(e2);
    
    % all four nodes must be different
    if a==c || a==d || b==c || b==d
        rej = rej+1;
        continue;
    end
    
    % choose randomly between two alternatives
    % 1) exchange b and d -> a-d, b-c
    % 2) exchange b and c -> a-c, b-d
    if rand > 0.5
        if ~(x(a,d) || x(b,c))
            x(a,b)=0; x(b,a)=0; x(c,d)=0; x(d,c)=0;
            x(a,d)=1; x(d,a)=1; x(b,c)=1; x(c,b)=1;
            j(e1)=d; j(e2)=b;
        else
            rej = rej+1;
            continue;
        end
    else
        if ~(x(a,c) || x(b,d))
            x(a,b)=0; x(b,a)=0; x(c,d)=0; x(d,c)=0;
            x(a,c)=1; x(c,a)=1; x(b,d)=1; x(d,b)=1;
            j(e1)=c; i(e2)=b;
        else
            rej = rej+1;
            continue;
        end
    end
    
    if nonuniform
        % update weights
        w1(e1) = degprod(i(e1),j(e1));
        w2(e1) = degprod_rev(i(e1),j(e1));
        w1(e2) = degprod(i(e2),j(e2));
        w2(e2) = degprod_rev(i(e2),j(e2));
    end
    
    eff = eff+1;
    rej = 0;
end
