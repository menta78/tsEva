function phandles = tsEvaPlotGEV3DFromAnalysisObj( X, nonStationaryEvaParams, stationaryTransformData, varargin )

  timeStamps = stationaryTransformData.timeStamps;
  epsilon = nonStationaryEvaParams(1).parameters.epsilon;
  sigma = nonStationaryEvaParams(1).parameters.sigma;
  mu = nonStationaryEvaParams(1).parameters.mu;

  phandles = tsEvaPlotGEV3D(X, timeStamps, epsilon, sigma, mu, varargin{:});
  
end

