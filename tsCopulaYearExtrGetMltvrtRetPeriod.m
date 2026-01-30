function [returnPeriod, prob] = tsCopulaYearExtrGetMltvrtRetPeriod(randomSample, level)
  % computes the multivariate return period according to 
  % Slavadori and De Michele 2004, Salvadori et al. 2011, 
  % used by Zscheischler et al. 2017

  ndim = size(randomSample, 2);
  if ndim ~= size(level, 2)
    error('tsCopulaComputeMultivariateRetPeriod: randomSample and retLev should have the same number of columns');
  end
  
  nretLev = size(level, 1);
  returnPeriod = zeros([nretLev, 1])*nan;
  prob = zeros([nretLev, 1])*nan;
  for iretLev = 1:nretLev
    lvli = level(iretLev, :);
    lvlmtx = lvli(ones([size(randomSample, 1), 1]), :);
    cnd = sum(randomSample > lvlmtx, 2) == ndim;
    prob(iretLev) = sum(cnd(:))/numel(cnd);
    returnPeriod(iretLev) = 1/prob(iretLev);
  end

  