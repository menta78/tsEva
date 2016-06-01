function phandles = tsEvaPlotSeriesTrendStdDevFromAnalyisObj( nonStationaryEvaParams, stationaryTransformData, varargin )
args.plotPercentile = -1;
args = tsEasyParseNamedArgs(varargin, args);
plotPercentile = args.plotPercentile;
  
timestamps = stationaryTransformData.timeStamps;
series = stationaryTransformData.nonStatSeries;
trend = stationaryTransformData.trendSeries;
stdDev = stationaryTransformData.stdDevSeries;

phandles = tsEvaPlotSeriesTrendStdDev(timestamps, series, trend, stdDev, varargin{:});

if plotPercentile ~= -1
  prcntile = tsEvaNanRunningPercentile(series, stationaryTransformData.runningStatsMulteplicity, plotPercentile);
  figure(phandles{1});
  hold on;
  hndl = plot(timestamps, prcntile);
  phandles{length(phandles) + 1} = hndl;
end
hold off;

end

