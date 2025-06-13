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
    paramEsts = mle(data, ...
        'pdf', gpdPdf, ...
        'start', [startK, startSigma], ...
        'lowerbound', [-Inf, 0], ...
        'upperbound', [0, Inf], ...
        'options', statset('MaxIter', 1e4, 'Display', 'off'));
    % covMatrix could be computed as hessian^-1 using the fact that dt ~ t-E(t) 
    covMatrix = mlecov(paramEsts, data, 'pdf', gpdPdf);

    % Compute standard errors
    se = sqrt(diag(covMatrix));

    % Compute confidence intervals
    z = norminv(1 - alphaCI / 2);
    k = paramEsts(1);
    sigma = paramEsts(2);
    kCI = [k - z * se(1), k + z * se(1)];
    sigmaCI = [sigma - z * se(2), sigma + z * se(2)];

    % Output structure matches gpfit
    paramCIs = [kCI; sigmaCI]';

end
