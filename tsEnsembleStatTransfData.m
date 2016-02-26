function [ stationaryTransformDataEnsemble ] = tsEnsembleStatTransfData( stationaryTransformDataArray )
  % From a cell array of non stationaryTransformData computes the average
  % stationaryTransformDataArray is a cell array of stationaryTransformData, 
  % the transformation object type returned by tsEvaNonStationary
  % ASSUMING THAT ALL THE SERIES HAVE THE SAME LENGTH 
  % AND UNIFORM TIME STAMPS
 
  n = length(stationaryTransformDataArray);
  for ii=1:n
    td = stationaryTransformDataArray{ii};
    if ii == 1
      etd.timeStamps = td.timeStamps;

      etd.trendSeries = td.trendSeries;
      etd.trendSeriesNonSeasonal = td.trendSeriesNonSeasonal;
      etd.trendError = td.trendError;
      etd.stdDevSeries = td.stdDevSeries;
      etd.stdDevSeriesNonSeasonal = td.stdDevSeriesNonSeasonal;
      etd.stdDevError = td.stdDevError;
    else
      etd.trendSeries = etd.trendSeries + td.trendSeries;
      etd.trendSeriesNonSeasonal = etd.trendSeriesNonSeasonal + td.trendSeriesNonSeasonal;
      etd.trendError = etd.trendError + td.trendError;
      etd.stdDevSeries = etd.stdDevSeries + td.stdDevSeries;
      etd.stdDevSeriesNonSeasonal = etd.stdDevSeriesNonSeasonal + td.stdDevSeriesNonSeasonal;
      etd.stdDevError = etd.stdDevError + td.stdDevError;
    end
  end
  
  etd.trendSeries = etd.trendSeries/n;
  etd.trendSeriesNonSeasonal = etd.trendSeriesNonSeasonal/n;
  etd.trendError = etd.trendError/n;
  etd.stdDevSeries = etd.stdDevSeries/n;
  etd.stdDevSeriesNonSeasonal = etd.stdDevSeriesNonSeasonal/n;
  etd.stdDevError = etd.stdDevError/n;
  
  stationaryTransformDataEnsemble = etd;

end

