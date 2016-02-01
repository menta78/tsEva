function [ trendSeries, filledTimeStamps, filledSeries, nRunMn ] = tsEvaRunningMeanTrend( timeStamps, series, timeWindow)
  indxs = ~isnan(series);
  timeStamps = timeStamps(indxs);
  series = series(indxs);

  mint = min(timeStamps);
  maxt = max(timeStamps);
  mindt = min(diff(timeStamps));
  filledTimeStamps = (mint:mindt:maxt)';
  filledSeries = interp1(timeStamps, series, filledTimeStamps, 'nearest');
  filledSeries = tsRemoveConstantSubseries(filledSeries, 4);
  
  nRunMn = ceil(timeWindow/mindt);
  trendSeries = tsEvaNanRunningMean(filledSeries, nRunMn);
  trendSeries = tsEvaNanRunningMean(trendSeries, ceil(nRunMn/2));
end

