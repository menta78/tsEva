function [gofStatistics] = tsCopulaUncertainty(copulaAnalysis)
%tsCopulaUncertainty estimation of copula goodness-of-fit
% [gofStatistics] = tsCopulaUncertainty(copulaAnalysis)
%                     returns a variable of type structure containing various parameters
%                     related to the goodness-of-fit of the copula



% A battery of goodnes-of-fit parameters are provided by tsCopulaUncertainty function.
% These include different assessment of the correlation parameters,
% Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC), Log-likelihood
% and a gof statistic calculated based on Cramer-Von mises statistic (CvM).

% input:
%  copulaAnalysis                           - a variable of type structure provided as the output of tsCopulaCompoundGPD or
%                                              tsCopulaCompoundGPDMontecarlo functions



% output:
%  gofStatistics:                           - A variable of type structure containing:
%                                               xxx                         --
%
%

% M.H.Bahmanpour, 2024

%   References:
%       [1] Mentaschi, L., Vousdoukas, M., Voukouvalas, E., Sartini, L., Feyen, L., Besio, G., and Alfieri, L.:
%           The transformed-stationary approach: a generic and simplified methodology for non-stationary extreme value analysis,
%           Hydrol. Earth Syst. Sci., 20, 3527–3547, https://doi.org/10.5194/hess-20-3527-2016, 2016
%       [2] Genest, C., Rémillard, B., Beaudoin, D., Goodness-of-fit tests for copulas:
%           A review and a power study (Open Access),(2009) Insurance:
%           Mathematics and Economics, 44 (2), pp. 199-213, doi: 10.1016/j.insmatheco.2007.10.005

%%%%%%%%%%%%%%%%%%%%%%

% Some parts of the code are reworked from https://github.com/mscavnicky/copula-matlab


% setting the default parameters

uResampled=copulaAnalysis.resampleProb; %from Monte-Carlo simulations
uSample=copulaAnalysis.jointExtremeMonovariateProbNS; %based on original samples


copulaFamily=copulaAnalysis.copulaParam.family;
copulaParam=copulaAnalysis.copulaParam;


if iscell(uResampled)
    %time-varying (non-stationary) copula
    [s1Sample,s2Sample]=cellfun(@(x) size(x),uSample);
    [s1Monte,~]=cellfun(@(x) size(x),uResampled);

    if strcmpi(copulaFamily, 'Gaussian')
        rho=copulaAnalysis.copulaParam.rho;
        yCDFSample=cellfun(@(x,y) copulacdf('Gaussian',x,y),uSample,rho,'UniformOutput',0);
        yCDFMonte=cellfun(@(x,y) copulacdf('Gaussian',x,y),uResampled,rho,'UniformOutput',0);

        copulaParam.numParams = s2Sample.*(s2Sample-1) / 2; %number of parameters of a multivariate Gaussian copula
    elseif strcmpi(copulaFamily, 't')
        rho=copulaAnalysis.copulaParam.rho;
        nu=copulaAnalysis.copulaParam.nu;
        yCDFSample=cellfun(@(x,y,z) copulacdf('t',x,y,z),uSample,rho,nu,'UniformOutput',0);
        yCDFMonte=cellfun(@(x,y,z) copulacdf('t',x,y,z),uResampled,rho,nu,'UniformOutput',0);
        copulaParam.numParams = 1 + s2Sample.*(s2Sample-1) / 2;
    elseif strcmpi(copulaFamily, 'Gumbel') || strcmpi(copulaFamily, 'Clayton') || strcmpi(copulaFamily, 'Frank')
        % in MATLAB, d-dimensional Gumbel, Clayton, and Frank copulas have
        % only one parameter currently; d-dimensional Archimedean copula
        % could be developed later as part of future works
        alpha=copulaAnalysis.copulaParam.rho;
        copulaFamilyc=repmat({copulaFamily},1,size(uSample,2));
        yCDFSample=cellfun(@(x,y,z) copulacdf(z,x,y),uSample,alpha,copulaFamilyc,'UniformOutput',0);
        yCDFMonte=cellfun(@(x,y,z) copulacdf(z,x,y),uResampled,alpha,copulaFamilyc,'UniformOutput',0);
        copulaParam.numParams = ones(1,size(uSample,2));
    end
    % else
    % Compute the log-likelihood

    llSample=cellfun(@(x) sum(log(x)),yCDFSample);
    llMonte=cellfun(@(x) sum(log(x)),yCDFMonte);


    k = copulaParam.numParams;
    % Compute the AIC
    aicSample = -2*llSample + (2*s1Sample.*k)./(s1Sample-k-1);
    aicMonte = -2*llMonte + (2*s1Monte.*k)./(s1Monte-k-1);
    % Compute the BIC
    bicSample = -2*llSample + k.*log(s1Sample);
    bicMonte = -2*llMonte + k.*log(s1Monte);
    eSample = tsRosenblattTransform(copulaParam, uSample );
    eMonte = tsRosenblattTransform(copulaParam, uResampled);

    % Compute the SnC statistics - SnC is a measure of the departure

    % see genest, et al., 2009, equation (9)
    sncSample=cellfun(@(x) sum((tsEmpirical(x) - prod(x, 2)) .^ 2),eSample,'UniformOutput',0);
    sncMonte=cellfun(@(x) sum((tsEmpirical(x) - prod(x, 2)) .^ 2),eMonte,'UniformOutput',0);

    % compute a range of correlation parameters; both for the samples and
    % for the resampled (or monte-carlo) values

    corrKendallSample=cellfun(@(x) corr(x,'type','Kendall'),uSample,'UniformOutput',0);
    corrSpearmanSample=cellfun(@(x) corr(x,'type','Spearman'),uSample,'UniformOutput',0);
    corrPearsonSample=cellfun(@(x) corr(x,'type','Pearson'),uSample,'UniformOutput',0);

    corrKendallMonte=cellfun(@(x) corr(x,'type','Kendall'),uResampled,'UniformOutput',0);
    corrSpearmanMonte=cellfun(@(x) corr(x,'type','Spearman'),uResampled,'UniformOutput',0);
    corrPearsonMonte=cellfun(@(x) corr(x,'type','Pearson'),uResampled,'UniformOutput',0);

    gofStatistics.sncSample=sncSample;
    gofStatistics.aicSample=aicSample;
    gofStatistics.bicSample=bicSample;
    gofStatistics.llSample=llSample;
    gofStatistics.sncMonte=sncMonte;
    gofStatistics.aicMonte=aicMonte;
    gofStatistics.bicMonte=bicMonte;
    gofStatistics.llMonte=llMonte;
    gofStatistics.corrKendallSample=corrKendallSample;
    gofStatistics.corrSpearmanSample=corrSpearmanSample;
    gofStatistics.corrPearsonSample=corrPearsonSample;
    gofStatistics.corrKendallMonte=corrKendallMonte;
    gofStatistics.corrSpearmanMonte=corrSpearmanMonte;
    gofStatistics.corrPearsonMonte=corrPearsonMonte;

else
    %stationary copula

    [s1Sample,s2Sample]=size(uSample);
    [S1Monte,~]=size(uResampled);
    if strcmpi(copulaFamily,'Gaussian')

        rho=copulaAnalysis.copulaParam.rho;
        yCDFSample=copulacdf('Gaussian',uSample,rho);
        yCDFMonte=copulacdf('Gaussian',uResampled,rho);
        %number of parameters of a multivariate Gaussian copula
        copulaParam.numParams = s2Sample.*(s2Sample-1) / 2;

    elseif strcmpi(copulaFamily, 't')

        rho=copulaAnalysis.copulaParam.rho;
        nu=copulaAnalysis.copulaParam.nu;
        yCDFSample=copulacdf('t',uSample,rho,nu);
        yCDFMonte=copulacdf('t',uResampled,rho,nu);
        copulaParam.numParams = 1 + s2Sample*(s2Sample-1) / 2;

    elseif strcmpi(copulaFamily, 'Gumbel') || strcmpi(copulaFamily, 'Clayton') || strcmpi(copulaFamily, 'Frank')

        alpha=copulaAnalysis.copulaParam.rho;
        yCDFSample=copulacdf(copulaFamily,uSample,alpha);
        yCDFMonte=copulacdf(copulaFamily,uResampled,alpha);
        copulaParam.numParams = 1;

    end

    % Compute the log-likelihood


    llSample = sum(log(yCDFSample));
    llMonte= sum(log(yCDFMonte));
    k = copulaParam.numParams;

    % Compute the AIC
    aicSample = -2*llSample + (2*s1Sample*k)/(s1Sample-k-1);
    aicMonte = -2*llMonte + (2*S1Monte*k)/(S1Monte-k-1);
    % Compute the BIC
    bicSample = -2*llSample + k*log(s1Sample);
    bicMonte = -2*llMonte + k*log(S1Monte);

    %compute SnC (a variant of CvM test Cramer-Von Mises)
    eSample = tsRosenblattTransform(copulaParam, uSample);
    eMonte = tsRosenblattTransform(copulaParam, uResampled);
    sncSample = sum((tsEmpirical(eSample) - prod(eSample, 2)) .^ 2);
    sncMonte = sum((tsEmpirical(eMonte) - prod(eMonte, 2)) .^ 2);

    corrKendallSample=corr(uSample,'type','Kendall');
    corrSpearmanSample=corr(uSample,'type','Spearman');
    corrPearsonSample=corr(uSample,'type','Pearson');
    corrKendallMonte=corr(uResampled,'type','Kendall');
    corrSpearmanMonte=corr(uResampled,'type','Spearman');
    corrPearsonMonte=corr(uResampled,'type','Pearson');

    gofStatistics.sncSample=sncSample;
    gofStatistics.sncMonte=sncMonte;
    gofStatistics.aicSample=aicSample;
    gofStatistics.aicMonte=aicMonte;
    gofStatistics.bicSample=bicSample;
    gofStatistics.bicMonte=bicMonte;
    gofStatistics.llSample=llSample;
    gofStatistics.llMonte=llMonte;
    gofStatistics.corrKendallSample=corrKendallSample;
    gofStatistics.corrKendallMonte=corrKendallMonte;
    gofStatistics.corrSpearmanSample=corrSpearmanSample;
    gofStatistics.corrSpearmanMonte=corrSpearmanMonte;
    gofStatistics.corrPearsonSample=corrPearsonSample;
    gofStatistics.corrPearsonMonte=corrPearsonMonte;
end

end


function [ Tc ] = tsRosenblattTransform( copulaparams, U )
%tsRosenblattTransform Rosenblatt's transformation
%
%[ Tc ] = tsRosenblattTransform( copulaparams, U )
%           returns a variable Tc corresponding to the Rosenblatt
%           transformation of probability values U
%
%
% this function performs Rosenblatt's probability integral under the null hypothesis that
% the data are generated using the given copula (as dictated by
% copulaparams)
%
%   References:
%       [1] Breymann, Dependence Structures for Multivariate High-Frequency
%       Data in Finance, 2003

% partly reworked from https://github.com/mscavnicky/copula-matlab



if iscell(U)

    [s1,s2]=cellfun(@(x) size(x),U,'UniformOutput',0);

    Tc=cellfun(@(x,x1) zeros(x,x1),s1,s2,'UniformOutput',0);

    uFirstColumn=cellfun(@(x) x(:,1),U,'UniformOutput',0);
    TcSecondColumn=cellfun(@(x,x1) x(:,2),Tc,'UniformOutput',0);
    Tc=cellfun(@(x,x1) [x,x1],uFirstColumn,TcSecondColumn,'UniformOutput',0);
    copulaparamsC={};

    rho=copulaparams.rho;
    if strcmpi(copulaparams.family,'t')
        nu=copulaparams.nu;
    end
    numPar=copulaparams.numParams;

    for ij=1:size(U,2)

        copulaparamx=copulaparams;
        copulaparamx.rho=rho{ij};
        copulaparamx.numParams=numPar(ij);
        if strcmpi(copulaparams.family,'t')
            copulaparamx.nu=nu{ij};
        end
        copulaparamsC=[copulaparamsC,{copulaparamx}];

    end

    d=unique([s2{:}]);

    for ii=2:d % for d-dimensional copula

        ithVectorCell=repmat({ii},1,size(U,2));
        conditCopulaCell=cellfun(@(x,y,z) tsCopulaCnd(x,y,z),copulaparamsC,U,ithVectorCell,'UniformOutput',0);
        Tc=cellfun(@(x,x1) [x(:,1),x1(:,1)],U,conditCopulaCell,'UniformOutput',0);

    end
else

    [s1,s2]=size(U);

    Tc = zeros(s1, s2);

    Tc(:,1) = U(:,1);

    for i=2:s2
        Tc(:,i) = tsCopulaCnd( copulaparams, U, i );
    end

end
end


function [ Y ] = tsCopulaCnd( copulaparams, U, m )
%tsCopulaCnd Conditional distribution function of different copula families.
%[ Y ] = tsCopulaCnd( copulaparams, U, m )
%           retuns a variable Y containing conditional distribution function
%   This function computes conditional CDF of d-dimensional copula, where m-th variable
%   is conditined upon the first m-1 variables.

% partly reworked from https://github.com/mscavnicky/copula-matlab

family = copulaparams.family;

% switch family
if strcmpi(family,'independent')
    Y = U(:,m);
elseif strcmpi(family,'gaussian')
    Y = tsGaussianCnd( copulaparams, U, m );
elseif strcmpi(family,'t')
    Y = tsStudentCnd( copulaparams, U, m );
elseif strcmpi(family,'frank') || strcmpi(family,'gumbel') || strcmpi(family,'clayton')
    Y = tsArchimCnd(family, U, copulaparams.rho, m);

end

end

function [ Y ] = tsGaussianCnd( copulaparams, U, m )
%
%tsGaussianCnd Computation of Gaussian conditional distribution function.
%[ Y ] = tsGaussianCnd( copulaparams, U, m )
%           retuns a variable Y containing gaussian conditional distribution function

% originally obtained from https://github.com/mscavnicky/copula-matlab

%   References:
%       [1] Wang (2012) - Numerical approximations and goodness-of-fit of
%       copulas

X = norminv(U(:, 1:m-1));
y = norminv(U(:, m));

rho = copulaparams.rho(1:m, 1:m);
irho = inv(rho);
add = irho(m,m);

edges =  X * (irho(1:m-1,m) + irho(m,1:m-1)');

H = sqrt(add) * y + edges / (2 * sqrt(add));
subA = sum((X * irho(1:m-1,1:m-1)) .* X, 2);
A = subA - edges.^2 / (4 * add);

t1 = exp(-0.5 * A) .* normcdf(H);
t2 = (2*pi)^(0.5 * (m-1)) * det(rho)^0.5 * sqrt(add);

N = t1 ./ t2;
D = mvnpdf( X, 0, rho(1:m-1,1:m-1) );

Y = N ./ D;

end

function [ Y ] = tsStudentCnd( copulaparams, U, m )

%tsStudentCnd Computation of Student-t conditional distribution function.
%[ Y ] = tsStudentCnd( copulaparams, U, m )
%           retuns a variable Y containing student-t conditional distribution function

% originally obtained from https://github.com/mscavnicky/copula-matlab

%   References:
%       [1] Wang (2012) - Numerical approximations and goodness-of-fit of
%       copulas


nu = copulaparams.nu;
X = tinv(U(:, 1:m-1), nu);
y = tinv(U(:, m), nu);

rho = copulaparams.rho(1:m,1:m);
irho = inv(rho);
add = irho(m,m);

edges =  X * (irho(1:m-1,m) + irho(m,1:m-1)');

subA = sum((X * irho(1:m-1,1:m-1)) .* X, 2);
A = subA - edges.^2 / (4 * add);

H = @(z,i) (1 + (((sqrt(add)*z + edges(i)/(2*sqrt(add))).^2 + A(i)) / nu)).^(-0.5*(nu+m));

t1 = gamma((nu+m)/2) / (gamma(nu/2) * (nu*pi)^(m/2) * det(rho)^0.5);
t2 = arrayfun(@(i) integral(@(z) H(z,i), -Inf, y(i)), (1:size(y))');

N = t1 * t2;
D = mvtpdf( X, rho(1:m-1,1:m-1), nu );
Y = N ./ D;

end

function [ Y ] = tsArchimCnd( family, U, alpha, m )

% tsArchimCnd Conditional distribution function for Archimedean copulas.
%[ Y ] = tsArchimCnd( family, U, alpha, m )
%           retuns a variable Y containing archimedean type copula conditional distribution function
% this function Computes conditional CDF of d-dimensional copula, where m-th variable
%   is conditined upon the first m-1 variables.

% originally obtained from https://github.com/mscavnicky/copula-matlab


X1 = sum(tsGeneratorInverse(family, U(:,1:m), alpha), 2);
N = tsGeneratorDerivative(family, X1, alpha, m-1);

X2 = sum(tsGeneratorInverse(family, U(:,1:m-1), alpha), 2);
D = tsGeneratorDerivative(family, X2, alpha, m-1);

Y = N ./ D;

end





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

function [ Y ] = tsGeneratorInverse( family, X, p )

% tsGeneratorInverse Inverse of archimedean copula generator
%[ Y ] = tsGeneratorInverse( family, X, p )
% originally obtained from https://github.com/mscavnicky/copula-matlab
%   Reference:
%       [1] Nelsen. R, (2006) Introduction to Copulas, Second Edition, page 116

if strcmpi(family,'clayton')

    Y = ( X .^ -p ) - 1;
elseif strcmpi(family,'gumbel')
    Y = ( -log(X) ) .^ p;
elseif strcmpi(family,'frank')
    Y = -log( ( exp(-p * X) - 1 ) / ( exp(-p) - 1 ) );
else
    error('Copula family %s not recognized.', family);
end

end

function [ Y ] = tsGeneratorDerivative( family, X, alpha, m )

% tsGeneratorDerivative derivative of the generator of the Archimedean copula
%[ Y ] = tsGeneratorDerivative( family, X, alpha, m )
%           Compute values of the m-th derivative of the
%           generator of the Archimedean copula family using numerical methods.
% originally obtained from https://github.com/mscavnicky/copula-matlab
%   References:
%       [1] Hofert (2011) - Likelihood Inference for Archimedean Copulas
[n, d] = size(X);

if strcmpi(family,'clayton')

    Y = prod((0:m-1) + (1/alpha)) * (1 + X).^(-(m+1/alpha));
elseif strcmpi(family,'Gumbel')
    a = zeros(m, 1);
    for i=1:m
        for j=i:m
            a(i) = a(i) + alpha^(-j) * tsStirling1(m, j) * tsStirling2(j, i);
        end
        a(i) = (-1)^(m-i) * a(i);
    end

    P = zeros(n, d);
    for i=1:m
        P = P + a(i, 1) * X.^(i / alpha);
    end

    Y = tsGenerator('Gumbel', X, alpha) ./ (X.^m) .* P;
elseif strcmpi(family,'Frank')
    Y = (1/alpha) * tsNpolylog(-m+1, (1-exp(-alpha)) * exp(-X));
end

Y = Y * (-1)^m;

end

function [ s ] = tsStirling1( m, n )

% tsStirling1 stirling number
%[ s ] = tsStirling1( m, n )
%            Returns stirling number of the first kind.
%            Caches numbers up to 32 dimensions.
% originally obtained from https://github.com/mscavnicky/copula-matlab

persistent cache;
if isempty(cache)
    cache = NaN(32, 32);
    % Put ones on diagonal
    cache(eye(32) == 1) = 1;
    % Zeros in the first column
    cache(2:32,1) = 0;
end

s = cache(m+1,n+1);
if isnan(s)
    s = stirling1(m-1, n-1) - (m-1)*stirling1(m-1, n);
    cache(m+1,n+1) = s;
end

end

function [ s ] = tsStirling2(m, n)

% tsStirling2 stirling number
%[ s ] = tsStirling2( m, n )
%            Returns stirling number of the second kind.
%            Caches numbers up to 32 dimensions.
% originally obtained from https://github.com/mscavnicky/copula-matlab

persistent cache;
if isempty(cache)
    cache = NaN(32, 32);
    % Put ones on diagonal
    cache(eye(32) == 1) = 1;
    % Zeros in the first column
    cache(2:32,1) = 0;
    % Ones in the second column
    cache(2:32,2) = 1;
end

s = cache(m+1,n+1);
if isnan(s)
    s = stirling2(m-1, n-1) + n*stirling2(m-1, n);
    cache(m+1,n+1) = s;
end

end

function [ Y ] = tsNpolylog( n, X )

% tsNpolylog Polylogratihm function
%[ Y ] = tsNpolylog( n, X )
%            Polylogratihm function for n < 0.

% originally obtained from https://github.com/mscavnicky/copula-matlab

%   References:
%       [1] http://en.wikipedia.org/wiki/Polylogarithm#Particular_values

assert(n <= 0, 'Polylogarithm only implemented for non-positive numbers.');

n = abs(n);
W = X ./ (1 - X);

Y = zeros(size(X));
for i=0:n
    Y = Y + factorial(i) * (W .^ (i+1)) * tsStirling2(n+1, i+1);
end

end


function [ Y ] = tsGenerator( family, X, p )

% tsGenerator Archimedean copula function generator
%[ Y ] = tsGenerator( family, X, p )
% originally obtained from https://github.com/mscavnicky/copula-matlab
%   Reference:
%       Nelsen. R, (2006) Introduction to Copulas, Second Edition, page 116

switch family
    case 'Clayton'
        Y = ( 1 + X ) .^ ( -1 / p );
    case 'Gumbel'
        Y = exp( -X .^ ( 1 / p ) );
    case 'Frank'
        Y = ( -1 / p ) * log( 1 - ( 1 - exp(-p) ) .* exp(-X) );
    otherwise
        error 'Copula family not recognized.'
end

end

