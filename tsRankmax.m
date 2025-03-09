function [ R ] = tsRankmax( X )
%tsRankmax Returns vector of one-based ranks for each element

% [ R ] = tsRankmax( X )
%   For the groups of same element, the maximum rank is returned. This can
%   be viewed as number of elements smaller or equal than given number.


% input:
%  X                           - a variable of type double containing
%                                samples in data space



% output:
%  R                          - A variable of type double containing
%                                 ranks of each element in X
%
%

% M.H.Bahmanpour, 2025

%REFERENCES

% [1] Bahmanpour, M.H., Mentaschi, L., Tilloy, A., Vousdoukas, M.,
%     Federico, I., Coppini, G., and Feyen, L., 2025,
%     Transformed-Stationary EVA 2.0: A Generalized Framework for
%     Non-stationary Joint Extreme Analysis (submitted to Hydrology and
%     Earth System Sciences; Feb 2025)
% [2] Mentaschi, L., Vousdoukas, M. I., Voukouvalas, E., Sartini, L.,
%     Feyen, L., Besio, G., & Alfieri, L. (2016). The
%     transformed-stationary approach: a generic and simplified methodology
%     for non-stationary extreme value analysis. Hydrology and Earth System
%     Sciences, 20(9), 3527–3547. https://doi.org/10.5194/hess-20-3527-2016
% [3] Genest, C., Rémillard, B., Beaudoin, D., Goodness-of-fit tests
%      for copulas: A review and a power study (Open Access),(2009) Insurance:
%      Mathematics and Economics, 44 (2), pp. 199-213, doi:
%      10.1016/j.insmatheco.2007.10.005
% [4] Hofert, M., Kojadinovic, I., Mächler, M. & Yan, J. Elements of
%      Copula Modeling with R (Springer, New York, 2018).
% [5] Berg, D. Bakken, H. (2006) Copula Goodness-of-fit Tests: A
%       Comparative Study
%%%%%%%%%%%%%%%%%%%%%%

% setting the default parameters

% Number of elements
n = size(X, 1);
% Preallocate ranks vector
R = zeros(n, 1);
% Sort the array and retrieve indices
[S, I] =  sort(X);
% Rank of the previous element
r = n;
% Value of the previous element
prev = S(n);

for i=n:-1:1
    x = S(i);
    if x == prev
        R(I(i)) = r;
    else
        prev = x;
        r = i;
        R(I(i)) = r;
    end
end

end