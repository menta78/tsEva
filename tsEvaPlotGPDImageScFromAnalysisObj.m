function phandles = tsEvaPlotGPDImageScFromAnalysisObj( Y, nonStationaryEvaParams, stationaryTransformData, varargin )

  timeStamps = stationaryTransformData.timeStamps;
  epsilon = nonStationaryEvaParams(2).parameters.epsilon;
  sigma = nonStationaryEvaParams(2).parameters.sigma;
  threshold = nonStationaryEvaParams(2).parameters.threshold;

  phandles = tsEvaPlotGPDImageSc(Y, timeStamps, epsilon, sigma, threshold, varargin{:});
  
end

