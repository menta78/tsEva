function [trasfData] = tsEvaTransformSeriesToStatSeasonal_ciPercentile( timeStamps, series, timeWindow, percentile, varargin  )
% this function decomposes the series into a season-dependent trend and a
% season-dependent standard deviation.
% The season-dependent standard deviation is given by a seasonal factor
% multiplied by a slowly varying standard deviation.

% transformation non stationary -> stationary
% x(t) = [y(t) - trend(t) - ssn_trend(t)]/[stdDev(t)*ssn_stdDev(t)]
% transformation stationary -> non stationary
% y(t) = stdDev(t)*ssn_stdDev(t)*x(t) + trend(t) + ssn_trend(t)

%trasfData.trendSeries = trend(t) + ssn_trend(t)
%trasfData.stdDevSeries = stdDev(t)*ssn_stdDev(t)

seasonalityTimeWindow = 2*30.4; % 2 months

disp('computing trend ...');
[statSeries, trendSeries, filledTimeStamps, filledSeries, nRunMn] = tsEvaDetrendTimeSeries(timeStamps, series, timeWindow, varargin{:});
disp('computing trend seasonality ...');
trendSeasonality = tsEstimateAverageSeasonality(filledTimeStamps, statSeries);
statSeries = statSeries - trendSeasonality;

disp(['computing the slowly varying ' num2str(percentile) 'th percentile ...']);
percentileSeries = tsEvaNanRunningPercentile(statSeries, nRunMn, percentile, varargin{:});
%
seasonalVarNRun = round(nRunMn/timeWindow*seasonalityTimeWindow);
%seasonalVarSeries is a moving variance computed on a short time
%window of 1-3 months, able to vary with season.
disp('computing standard deviation seasonality ...');
seasonalVarSeries = tsEvaNanRunningVariance(statSeries, seasonalVarNRun);
seasonalStdDevSeries = sqrt(seasonalVarSeries./percentileSeries);
seasonalStdDevSeries = tsEstimateAverageSeasonality(filledTimeStamps, seasonalStdDevSeries);
stdDevSeriesNonSeasonal = sqrt(percentileSeries);

statSeries = statSeries./(stdDevSeriesNonSeasonal.*seasonalStdDevSeries);
[~, ~, statSer3Mom, statSer4Mom] = tsEvaNanRunningStatistics(statSeries, nRunMn);
statSer3Mom = tsEvaNanRunningMean(statSer3Mom, ceil(nRunMn));
statSer4Mom = tsEvaNanRunningMean(statSer4Mom, ceil(nRunMn));

% N is the size of each sample used to compute the average
N = nRunMn;
% the error on the trend is computed as the error on the average:
%  error(average) = stdDev/sqrt(N)
trendNonSeasonalError = nanmean(stdDevSeriesNonSeasonal)/N^.5;

% computation of the error on the standard deviation explained in 
% Mentaschi et al 2016
S = 2;
avgStdDev = nanmean(stdDevSeriesNonSeasonal);
stdDevNonSeasonalError = avgStdDev*( 2*S^2/N^3 )^(1./4.);

Ntot = length(series);
trendSeasonalError = stdDevNonSeasonalError*sqrt(12/Ntot + 1/N);
stdDevSeasonalError = seasonalStdDevSeries*(288./Ntot^2/N)^(1./4.);

trendError = sqrt(trendNonSeasonalError^2 + trendSeasonalError^2);
stdDevError = sqrt( (stdDevSeriesNonSeasonal.*stdDevSeasonalError).^2 + (seasonalStdDevSeries.*stdDevNonSeasonalError).^2 );

trasfData.runningStatsMulteplicity = nRunMn;
trasfData.stationarySeries = statSeries;
trasfData.trendSeries = trendSeries + trendSeasonality;
trasfData.trendSeriesNonSeasonal = trendSeries;
trasfData.stdDevSeries = stdDevSeriesNonSeasonal.*seasonalStdDevSeries;
trasfData.stdDevSeriesNonSeasonal = stdDevSeriesNonSeasonal;

trasfData.trendNonSeasonalError = trendNonSeasonalError;
trasfData.stdDevNonSeasonalError = stdDevNonSeasonalError;
trasfData.trendSeasonalError = trendSeasonalError;
trasfData.stdDevSeasonalError = stdDevSeasonalError;

trasfData.trendError = trendError;
trasfData.stdDevError = stdDevError;

trasfData.timeStamps = filledTimeStamps;
trasfData.nonStatSeries = filledSeries;
trasfData.statSer3Mom = statSer3Mom;
trasfData.statSer4Mom = statSer4Mom;
end

