function [retPer, exceedProb] = tsEvaGetReturnPeriodOfLevelGPD( epsilon, sigma, pPeak, retLev)
  if epsilon ~= 0
    H = 1 - (1 + epsilon.*retLev./sigma).^(-1./epsilon);
  else
    H = 1 - exp(-retLev./sigma);
  end
  exceedProb = (1 - H)*pPeak;
  retPer = 1./exceedProb;
end
