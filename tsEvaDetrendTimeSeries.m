function [ detrendSeries, trendSeries, filledTimeStamps, filledSeries, nRunMn ] = tsEvaDetrendTimeSeries( timeStamps, series, timeWindow, varargin )

args.extremeLowThreshold = -Inf;
args = tsEasyParseNamedArgs(varargin, args);
extremeLowThreshold = args.extremeLowThreshold;

[trendSeries, filledTimeStamps, filledSeries, nRunMn] = tsEvaRunningMeanTrend(timeStamps, series, timeWindow);
statSeries = filledSeries;
statSeries(statSeries < extremeLowThreshold) = nan;
detrendSeries = statSeries - trendSeries;

end

