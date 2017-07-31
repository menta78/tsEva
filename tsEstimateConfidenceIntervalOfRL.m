function [ lowCI, highCI ] = tsEstimateConfidenceIntervalOfRL( rl, stdErr, p )
  %% tsEstimateConfidenceIntervalOfRL: 
  % estimates the confidence interval of the return level (rl), assuming that the
  % uncertainty distribution around rl is a lognormal with standard deviation
  % stdErr
  %%

  m = rl;
  vr = stdErr.^2;

  mu = log((m.^2)./sqrt(vr + m.^2));
  sigma = sqrt(log(vr./(m.^2) + 1));

  highCI = logninv(p, mu, sigma);
  lowCI = logninv((1 - p), mu, sigma);
  

end

