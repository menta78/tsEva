function [ C ] = tsEmpirical( U )

% tsEmpirical Empirical copula
%[ C ] = tsEmpirical( U )
%           retuns a variable C containing the empirical copula
%           corresponding with the U variable including the uniform
%           variates

% originally obtained from https://github.com/mscavnicky/copula-matlab

[n,d] = size(U);

C = zeros(n, 1);
for i=1:n
    S = ones(n, 1);

    for j=1:d
        S = S .* (U(:, j) <= U(i, j));
    end
    C(i) = sum(S) / n;
end

end



