function [C, se] = tsCopulaCdfFromSamples(u, Usample)
% tsCopulaCdfFromSamples  Empirical copula CDF from sample points
% u       : (q x d) query points in [0,1]^d
% Usample : (M x d) sample points from the fitted copula
% C       : (q x 1) empirical CDF estimates
% se      : (q x 1) standard errors sqrt(C(1-C)/M) [optional]

    if nargin < 2, error('Need u (qxd) and Usample (Mxd).'); end
    [q,d] = size(u);
    if size(Usample,2) ~= d
        error('Dimension mismatch: size(u,2)=%d, size(Usample,2)=%d', d, size(Usample,2));
    end

    epsv = 1e-12;
    u       = min(max(u,       epsv), 1-epsv);
    Usample = min(max(Usample, epsv), 1-epsv);

    M    = size(Usample,1);
    U3   = reshape(Usample, [M, d, 1]);  
    u3   = permute(u, [3, 2, 1]);        

    mask = all(U3 <= u3, 2);             
    C    = squeeze(mean(mask, 1));    

    if nargout > 1
        se = sqrt(C .* max(1 - C, 0) / M);
    end
end