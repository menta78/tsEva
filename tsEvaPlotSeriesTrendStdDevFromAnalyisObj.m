function phandles = tsEvaPlotSeriesTrendStdDevFromAnalyisObj( nonStationaryEvaParams, stationaryTransformData, varargin )
  
timestamps = stationaryTransformData.timeStamps;
series = stationaryTransformData.nonStatSeries;
trend = stationaryTransformData.trendSeries;
stdDev = stationaryTransformData.stdDevSeries;

phandles = tsEvaPlotSeriesTrendStdDev(timestamps, series, trend, stdDev, varargin{:});

end

