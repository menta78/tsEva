addpath('../');

load('timeAndSeriesHebrides.mat');
timeAndSeries = timeAndSeriesHebrides;
timeWindow = 365.25*6; % 6 years
percent = 80;

timeStamps = timeAndSeries(:,1);
series = timeAndSeries(:,2);

[ trendSeries, filledTimeStamps, filledSeries, nRunMn ] = tsEvaRunningMeanTrend( timeStamps, series, timeWindow);
[ rnprcnt ] = tsEvaNanRunningPercentile( filledSeries, nRunMn, percent );
