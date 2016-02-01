function [ trendSeries, filledTimeStamps, filledSeries ] = runningMeanTrend( timeStamps, series, timeWindow)
  mint = min(timeStamps);
  maxt = max(timeStamps);
  mindt = min(diff(timeStamps));
  filledTimeStamps = (mint:mindt:maxt);
  filledSeries = interp1(timeStamps, series, filledTimeStamps, 'nearest');
  filledSeries = tsRemoveConstantSubseries(filledSeries, 4);
  
  nRunMn = ceil(timeWindow/mindt);
  trendSeries = nanRunningMean(filledSeries, nRunMn);
end

