function [gofStatistics] = tsCopulaUncertainty(resampProb,copulaAnalysis)
% a function to assess goodness-of-fit of the Monte-Carlo simulated
% values
% correlation parameters (Kendall, Spearman and Pearson) are calculated
% using the probability values

% resampProb is the resampled probability obtained from Monte-Carlo
% simulation

% copulaAnalysis is the copulaAnalysis structure file containing details of
% fitted copula and also information about marginals and sampled extremes

% for a time-varying copula, i.e., with resampProb in a cell data type,
% gofStatistics are calculated for all time segments 

% Reworked from https://github.com/mscavnicky/copula-matlab

% M. H. Bahmanpour, 2023

copulaFamily=copulaAnalysis.copulaParam.family;
copulaParam=copulaAnalysis.copulaParam;

U = resampProb;
if iscell(U) %time-varying copula
    [nc,dc]=cellfun(@(x) size(x),U);
    if strcmpi(copulaFamily, 'Gaussian')
        rho=copulaAnalysis.copulaParam.rho;
        Y=cellfun(@(x,y) copulapdf('Gaussian',x,y),U,rho,'UniformOutput',0);
        copulaParam.numParams = dc.*(dc-1) / 2; %number of parameters of a multivariate Gaussian copula
    elseif strcmpi(copulaFamily, 't')
        rho=copulaAnalysis.copulaParam.rho;
        nu=copulaAnalysis.copulaParam.nu;
        Y=cellfun(@(x,y,z) copulapdf('t',x,y,z),U,rho,nu,'UniformOutput',0);
        copulaParam.numParams = 1 + dc.*(dc-1) / 2;
    elseif strcmpi(copulaFamily, 'Gumbel') || strcmpi(copulaFamily, 'Clayton') || strcmpi(copulaFamily, 'Frank')
        alpha=copulaAnalysis.copulaParam.rho;
        copulaFamilyc=repmat({copulaFamily},1,size(U,2));
        Y=cellfun(@(x,y,z) copulapdf(z,x,y),U,alpha,copulaFamilyc,'UniformOutput',0);
        copulaParam.numParams = ones(1,size(U,2)); %Gumber, Clayton, and Frank copulas have one parameter in  case of d-dimensioanl
                                                   % copula future works,should proparly account for d-dimensional Archimedean copula
    end
    % else
    % Compute the log-likelihood

    ll=cellfun(@(x) sum(log(x)),Y);
 

    k = copulaParam.numParams;
    % Compute the AIC
    aic = -2*ll + (2*nc.*k)./(nc-k-1);
    % Compute the BIC
    bic = -2*ll + k.*log(nc);

    E = rosenblattTransform(copulaParam, U );
    % Produce vector with chi-square distribution
    
    C=cellfun(@(x) sum(norminv(x).^ 2,2),E,'UniformOutput',0);
    % Compute the AKS statistics
    aks=cellfun(@(x,y,y2) sum(abs(chi2cdf(x,y) - pseudoObservations(x))) / sqrt(y2),C,mat2cell(dc,1,ones(1,size(U,2))),mat2cell(nc,1,ones(1,size(U,2))),'UniformOutput',0);
   
    % Compute the SnC statistics

    snc=cellfun(@(x) sum((empirical(x) - prod(x, 2)) .^ 2),E,'UniformOutput',0);
    
    gofStatistics.snc=snc;  % since its a measure of departure, the smaller the better
    gofStatistics.aks=aks; % since its a measure of departure, the smaller the better
    gofStatistics.aic=aic;  % the smaller the better
    gofStatistics.bic=bic;  %the smaller the better
    gofStatistics.ll=ll; %the largest value represent the highest likelihood

    corrKendall=cellfun(@(x) corr(x,'type','Kendall'),U,'UniformOutput',0);
    corrSpearman=cellfun(@(x) corr(x,'type','Spearman'),U,'UniformOutput',0);
    corrPearson=cellfun(@(x) corr(x,'type','Pearson'),U,'UniformOutput',0);

    gofStatistics.corrKendall=corrKendall;
    gofStatistics.corrSpearman=corrSpearman;
    gofStatistics.corrPearson=corrPearson;
else
    [n,d]=size(U);
    if strcmpi(copulaFamily, 'Gaussian')

        rho=copulaAnalysis.copulaParam.rho;
        Y=copulapdf('Gaussian',U,rho);
        copulaParam.numParams = d.*(d-1) / 2; %number of parameters of a multivariate Gaussian copula

    elseif strcmpi(copulaFamily, 't')

        rho=copulaAnalysis.copulaParam.rho;
        nu=opulaAnalysis.copulaParam.nu;
        Y=copulapdf('t',U,rho,nu);
        copulaParam.numParams = 1 + d*(d-1) / 2;

    elseif strcmpi(copulaFamily, 'Gumbel') || strcmpi(copulaFamily, 'Clayton') || strcmpi(copulaFamily, 'Frank')

        alpha=copulaAnalysis.copulaParam.rho;
        Y=copulapdf(copulaFamily,U,alpha);
        copulaParam.numParams = 1;

    end

    % Compute the log-likelihood

    ll = sum(log(Y));

    k = copulaParam.numParams;
    % Compute the AIC
    aic = -2*ll + (2*n*k)/(n-k-1);
    % Compute the BIC
    bic = -2*ll + k*log(n);

    E = rosenblattTransform(copulaParam, U );
    % Produce vector with chi-square distribution
    C = sum( norminv( E ) .^ 2, 2 );

    % Compute the AKS statistics
    aks = sum(abs(chi2cdf(C, d) - pseudoObservations(C))) / sqrt(n);
    % Compute the SnC statistics

    snc = sum((empirical(E) - prod(E, 2)) .^ 2);

    gofStatistics.snc=snc;  % since its a measure of departure, the smaller the better
    gofStatistics.aks=aks; % since its a measure of departure, the smaller the better
    gofStatistics.aic=aic;  % the smaller the better
    gofStatistics.bic=bic;  %the smaller the better
    gofStatistics.ll=ll; %the largest value represent the highest likelihood

    corrKendall=corr(U,'type','Kendall');
    corrSpearman=corr(U,'type','Spearman');
    corrPearson=corr(U,'type','Pearson');

    gofStatistics.corrKendall=corrKendall;
    gofStatistics.corrSpearman=corrSpearman;
    gofStatistics.corrPearson=corrPearson;
end

end


function [ Tc ] = rosenblattTransform( copulaparams, U )
%rosenblattTransform Performs Rosenblatt's probability integral
%transformation under the null hypothesis that the data are generated using
%the given copula.
%
%   References:
%       [1] Breymann, Dependence Structures for Multivariate High-Frequency
%       Data in Finance, 2003
%reworked from https://github.com/mscavnicky/copula-matlab

if iscell(U)

    [nc,dc]=cellfun(@(x) size(x),U,'UniformOutput',0);
    
    Tc=cellfun(@(x,x1) zeros(x,x1),nc,dc,'UniformOutput',0);
    Tci=[Tc{:,:}];
    Ui=[U{:,:}];
    Tci(:,1:2:end)=Ui(:,1:2:end);

    Tc = mat2cell(Tci,size(Tci,1),2.*ones(1,size(U,2)));
   
    copulaparamsC={};

    for ij=1:size(U,2)
        rhox=copulaparams.rho;
        if strcmp(copulaparams.family,'t')
        nux=copulaparams.nu;
        end
        numparx=copulaparams.numParams;
        copulaparamx=copulaparams;
        copulaparamx.rho=rhox{ij};
        copulaparamx.numParams=numparx(ij);
         if strcmp(copulaparams.family,'t')
        copulaparamx.nu=nux{ij};
         end
        copulaparamsC=[copulaparamsC,{copulaparamx}];
        
    end

    d=unique([dc{:}]);

    for ii=2:d % for d-dimensional copula

        ii2=repmat({ii},1,size(U,2));
        Icc=cellfun(@(x,y,z) copulacnd(x,y,z),copulaparamsC,U,ii2,'UniformOutput',0);
     
        Tci=[Tc{:,:}];
        Icci=[Icc{:,:}];
        Tci(:,2:2:end)=Icci(:,1:1:end);
        Tc = mat2cell(Tci,size(Tci,1),2.*ones(1,size(U,2)));

    end
else

    [n,d]=size(U);

    Tc = zeros(n, d);

    Tc(:,1) = U(:,1);

    for i=2:d
        Tc(:,i) = copulacnd( copulaparams, U, i );
    end

end
end


function [ Y ] = copulacnd( copulaparams, U, m )
%copulacnd Conditional distribution function of different copula families.
%   Computes conditional CDF of d-dimensional copula, where m-th variable
%   is conditined upon the first m-1 variables.
% reworked from https://github.com/mscavnicky/copula-matlab
family = copulaparams.family;

switch family

    case 'independent'
        % Conditioning upon independent variables
        Y = U(:,m);
    case 'Gaussian'
        Y = gaussianCnd( copulaparams, U, m );
    case 't'
        Y = studentCnd( copulaparams, U, m );
    case {'Frank', 'Gumbel', 'Clayton'}
        Y = archimcnd(family, U, copulaparams.rho, m);

end

end

function [ Y ] = gaussianCnd( copulaparams, U, m )
%GAUSSIANCND Computation of Gaussian conditional distribution function.
%
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

function [ Y ] = studentCnd( copulaparams, U, m )
%STUDENTCND Computation of Student-t conditional distribution function.
%
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

function [ Y ] = archimcnd( family, U, alpha, m )
% ARCHIM.CND Conditional distribution function for Archimedean copulas.
%   Computes conditional CDF of d-dimensional copula, where m-th variable
%   is conditined upon the first m-1 variables.

X1 = sum(generatorInverse(family, U(:,1:m), alpha), 2);
N = generatorDerivative(family, X1, alpha, m-1);

X2 = sum(generatorInverse(family, U(:,1:m-1), alpha), 2);
D = generatorDerivative(family, X2, alpha, m-1);

Y = N ./ D;

end

function [ U ] = pseudoObservations( X )
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
    U(:,i) = rankmax(X(:,i)) / (n + 1);
end

end

function [ R ] = rankmax( X )
%RANKMAX Returns vector of one-based ranks for each element
%   For the groups of same element, the maximum rank is returned. This can
%   be viewed as number of elements smaller or equal than given number.

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

function [ C ] = empirical( U )
%COPULA.EMPIRICAL Empirical copula for d-dimensional data.
%   Works with uniform variates.
%
%   References:
%       [1] Genest, C. (2009) Goodness-of-fit tests for copulas: A review
%       and a power study

[n d] = size(U);

C = zeros(n, 1);
for i=1:n
    S = ones(n, 1);
    % Make AND using logical bitmaps of each column
    for j=1:d
        S = S .* (U(:, j) <= U(i, j));
    end
    C(i) = sum(S) / n;
end

end

function [ Y ] = generatorInverse( family, X, p )
%ARCHIM.GENERATORINVERSE Inverse of archimedean copula generator.
%   Please note that no parameter checking is done on this level.
%
%   Reference:
%       [1] Nelsen. R, (2006) Introduction to Copulas, Second Edition, page 116

switch family
    case 'Clayton'
        Y = ( X .^ -p ) - 1;
    case 'Gumbel'
        Y = ( -log(X) ) .^ p;
    case 'Frank'
        Y = -log( ( exp(-p * X) - 1 ) / ( exp(-p) - 1 ) );
    otherwise
        error('Copula family %s not recognized.', family);
end

end

function [ Y ] = generatorDerivative( family, X, alpha, m )
%ARCHIM.GENERATORDERIVATIVE Compute values of the m-th derivative of the
%generator of the Archimedean copula family using numerical methods.
%
%   References:
%       [1] Hofert (2011) - Likelihood Inference for Archimedean Copulas

[n, d] = size(X);

switch family
case 'Clayton' 
    Y = prod((0:m-1) + (1/alpha)) * (1 + X).^(-(m+1/alpha));
case 'Gumbel'
    a = zeros(m, 1);
    for i=1:m
       for j=i:m
          a(i) = a(i) + alpha^(-j) * stirling1(m, j) * stirling2(j, i);
       end        
       a(i) = (-1)^(m-i) * a(i); 
    end    
    
    P = zeros(n, d);
    for i=1:m
        P = P + a(i, 1) * X.^(i / alpha);
    end    
    
    Y = generator('Gumbel', X, alpha) ./ (X.^m) .* P;  
case 'Frank'
    Y = (1/alpha) * npolylog(-m+1, (1-exp(-alpha)) * exp(-X));
end

Y = Y * (-1)^m;

end

function [ s ] = stirling1( m, n )
%STIRLING1 Returns stirling number of the first kind.
%   Caches numbers up to 32 dimensions.

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

function [ s ] = stirling2(m, n)
%STIRLING2 Computes stirling number of the seconds kind.
%   Caches numbers up to 32 dimensions.

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

function [ Y ] = npolylog( n, X )
%NPOLYLOG Polylogratihm function for n < 0.
%
%   References:
%       [1] http://en.wikipedia.org/wiki/Polylogarithm#Particular_values

assert(n <= 0, 'Polylogarithm only implemented for non-positive numbers.');

n = abs(n);
W = X ./ (1 - X);

Y = zeros(size(X));
for i=0:n
    Y = Y + factorial(i) * (W .^ (i+1)) * stirling2(n+1, i+1);
end

end


function [ Y ] = generator( family, X, p )
%ARCHIM.GENERATOR Archimedean copula function generator.
%   Please note that no parameter checking is done in this function.
%
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