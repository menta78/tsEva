function [retPer, exceedProb] = tsEvaGetReturnPeriodOfLevelGEV( epsilon, sigma, mu, retLev)
  if epsilon ~= 0
    % GEV
    G = exp(-(1 + epsilon.*(retLev - mu)/sigma).^(-1/epsilon));
  else
    % Gumbel
    G = exp(-exp(-(retLev-mu)./sigma));
  end
  exceedProb = 1 - G;
  retPer = 1./exceedProb;
end
    