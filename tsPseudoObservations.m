function [ U ] = tsPseudoObservations( X )
%PSEUDOOBSERVATIONS Uniforms input sample to pseudo-observations.
%   Based on empirical CDF function described in [1]. We use n+1 for
%   division to keep empirical CDF lower than 1.
%
%   References:
%       [1] Berg, D. Bakken, H. (2006) Copula Goodness-of-fit Tests: A
%       Comparative Study

[n, d] = size(X);
U = zeros(n, d);

for i=1:d
    U(:,i) = tsRankmax(X(:,i)) / (n + 1);
end

end