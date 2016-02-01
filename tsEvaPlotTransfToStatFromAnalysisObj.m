function phandles = tsEvaPlotTransfToStatFromAnalysisObj( nonStationaryEvaParams, stationaryTransformData, varargin )
  timestamps = stationaryTransformData.timeStamps;
  series = stationaryTransformData.stationarySeries;
  srmean = zeros(size(series));
  srstddev = ones(size(series));
%  std3mom = nthroot(stationaryTransformData.statSer3Mom, 3.);
%  std4mom = nthroot(stationaryTransformData.statSer4Mom, 4.);
  st3mom = stationaryTransformData.statSer3Mom;
  st4mom = stationaryTransformData.statSer4Mom;
  phandles = tsEvaPlotTransfToStat(timestamps, series, srmean, srstddev, st3mom, st4mom, varargin{:});
end

