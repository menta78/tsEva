function u = tsCopulaRnd(family, copulaPar, N, uProb)

% generates the random vector from a copula distribution.
% if multivariate t or gaussian or bivariate archimedean, it uses copularnd
% else (this would be the case of multivariate archimedean) it implements a
% simple C-Vine copula
%
% parameters:
%   family: copula family
%   copulaPar: rho for gaussian, alpha for archimedean
%   N: sample size
%   uProb: pseudo-observation data (that is, transformed in uniform probability space) 
%   from which the copula was built. Necessary in case the
%   c-vine copula must be constructed

if strcmpi(family, 't')
    error("copula family unsupported: t");
end

cndCopulaRnd0 = strcmpi(family, 'gaussian');
cndCopulaRnd1 = isscalar(copulaPar) || size(copulaPar, 1) == 2;
cndCopulaRnd = cndCopulaRnd0 | cndCopulaRnd1;

if cndCopulaRnd
    if ~isscalar(copulaPar), copulaPar = copulaPar(1,2); end
    u = copularnd(family, copulaPar, N);
else
    if ~strcmpi(family, 'gumbel')
        error("multivariate archimedean only supported for gumbel");
    end
    % Resampling via a C-vine copula:
    % conceptually equivalent to sequentially sampling one variable at a time,
    % estimating the conditional distribution at each step (using pair-copulas)
    % and drawing the next variable conditioned on the previous ones.

    alpha = copulaPar;
    order = tsGumbelCVine.cvineOrder(alpha);          % root-first ordering
    theta = tsGumbelCVine.fit(uProb, alpha, order);   % estimate θ on all trees (pseudo-obs)
    u = tsGumbelCVine.simulate(N, order, theta);      % sample (the sim I gave you)
end
end
