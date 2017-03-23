function [trasfData] = tsEvaTransformSeriesToStationaryTrendOnly_ciPercentile( timeStamps, series, timeWindow, percentile, varargin )

disp('computing the trend ...');
[statSeries, trendSeries, filledTimeStamps, filledSeries, nRunMn] = tsEvaDetrendTimeSeries(timeStamps, series, timeWindow, varargin{:});

disp(['computing the slowly varying ' num2str(percentile) 'th percentile ...']);
[percentileSeries, stdErr] = tsEvaNanRunningPercentile(statSeries, nRunMn, percentile, varargin{:});
meanPerc = nanmean(percentileSeries);
%normalizing to standard deviation (just to be able to make acceptable graphs with the scripts of this library)
stdDev = nanstd(statSeries);
stdDevSeries = percentileSeries/meanPerc*stdDev;
stdDevError = stdErr/meanPerc*stdDev;

statSeries = statSeries./stdDevSeries;
[~, ~, statSer3Mom, statSer4Mom] = tsEvaNanRunningStatistics(statSeries, nRunMn);
statSer3Mom = tsEvaNanRunningMean(statSer3Mom, ceil(nRunMn));
statSer4Mom = tsEvaNanRunningMean(statSer4Mom, ceil(nRunMn));

% N is the size of each sample used to compute the average
N = nRunMn;
% the error on the trend is computed as the error on the average:
%  error(average) = stdDev/sqrt(N)
trendError = nanmean(stdDevSeries)/N^.5;

trasfData.runningStatsMulteplicity = nRunMn;
trasfData.stationarySeries = statSeries;
trasfData.trendSeries = trendSeries;
trasfData.trendSeriesNonSeasonal = trendSeries;
trasfData.trendError = trendError;
trasfData.stdDevSeries = stdDevSeries;
trasfData.stdDevSeriesNonSeasonal = stdDevSeries;
trasfData.stdDevError = stdDevError;
trasfData.timeStamps = filledTimeStamps;
trasfData.nonStatSeries = filledSeries;
trasfData.statSer3Mom = statSer3Mom;
trasfData.statSer4Mom = statSer4Mom;
end

