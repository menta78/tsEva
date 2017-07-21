function [trasfData] = tsEvaTransformSeriesToStationaryTrendOnly( timeStamps, series, timeWindow, varargin )

disp('computing the trend ...');
[statSeries, trendSeries, filledTimeStamps, filledSeries, nRunMn] = tsEvaDetrendTimeSeries(timeStamps, series, timeWindow, varargin{:});

disp('computing the slowly varying standard deviation ...');
varianceSeries = tsEvaNanRunningVariance(statSeries, nRunMn);
%further smoothing
varianceSeries = tsEvaNanRunningMean(varianceSeries, ceil(nRunMn/2));
%
stdDevSeries = varianceSeries.^.5;
statSeries = statSeries./stdDevSeries;
[~, ~, statSer3Mom, statSer4Mom] = tsEvaNanRunningStatistics(statSeries, nRunMn);
statSer3Mom = tsEvaNanRunningMean(statSer3Mom, ceil(nRunMn));
statSer4Mom = tsEvaNanRunningMean(statSer4Mom, ceil(nRunMn));

% N is the size of each sample used to compute the average
N = nRunMn;
% the error on the trend is computed as the error on the average:
%  error(average) = stdDev/sqrt(N)
trendError = nanmean(stdDevSeries)/N^.5;

% variance(stdDev) ~ 2 stdDev^4 / (n - 1)

% stdDevError is approximated as constant
avgStdDev = nanmean(stdDevSeries);
S = 2;
% computation of the error on the standard deviation explained in 
% Mentaschi et al 2016
stdDevError = avgStdDev*( 2*S^2/N^3 )^(1./4.);

trasfData.runningStatsMulteplicity = nRunMn;
trasfData.stationarySeries = statSeries;
trasfData.trendSeries = trendSeries;
trasfData.trendSeriesNonSeasonal = trendSeries;
trasfData.trendError = trendError;
trasfData.stdDevSeries = stdDevSeries;
trasfData.stdDevSeriesNonSeasonal = stdDevSeries;
trasfData.stdDevError = stdDevError*ones(size(stdDevSeries));
trasfData.timeStamps = filledTimeStamps;
trasfData.nonStatSeries = filledSeries;
trasfData.statSer3Mom = statSer3Mom;
trasfData.statSer4Mom = statSer4Mom;
end

