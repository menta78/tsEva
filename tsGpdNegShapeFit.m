function [paramEsts, paramCIs] = tsGpdNegShapeFit(data, alphaCI)
%TSGPDNEGSHAPEFIT Fit GPD with negative shape constraint (k < 0)
%   [paramEsts, paramCIs] = tsGpdNegShapeFit(data, alphaCI)
%   mimics gpfit but constrains the shape parameter to be negative.

    if nargin < 2
        alphaCI = 0.05;  % Default 95% CI
    end

    % Define the distribution log-likelihood
    gpdPdf = @(x, k, sigma) max(gppdf(x, k, sigma), realmin);

    % Initial guesses
    startK = -0.1;
    startSigma = std(data) / sqrt(2);  % Crude scale estimate

    % Fit using MLE with constraint k < 0
    % [paramEsts, ~, exitflag, output, ~, hessian] 
    [paramEsts, paramCIs] = mle(data, ...
        'pdf', gpdPdf, ...
        'start', [startK, startSigma], ...
        'lowerbound', [-Inf, 0], ...
        'upperbound', [0, Inf], ...
        'alpha', alphaCI, ...
        'options', statset('MaxIter', 1e4, 'Display', 'off'));
end
