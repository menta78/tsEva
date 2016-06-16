function [ trendSeries, filledTimeStamps, filledSeries, nRunMn ] = tsEvaRunningMeanTrend( timeStamps, series, timeWindow)
[ filledTimeStamps, filledSeries, dt ] = tsEvaFillSeries( timeStamps, series );  
  nRunMn = ceil(timeWindow/dt);
  trendSeries = tsEvaNanRunningMean(filledSeries, nRunMn);
  trendSeries = tsEvaNanRunningMean(trendSeries, ceil(nRunMn/2));
end

