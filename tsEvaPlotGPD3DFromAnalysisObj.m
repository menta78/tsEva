function phandles = plotGPD3DFromAnalysisObj( X, nonStationaryEvaParams, stationaryTransformData )

  timeStamps = stationaryTransformData.timeStamps;
  epsilon = nonStationaryEvaParams(2).parameters.epsilon;
  sigma = nonStationaryEvaParams(2).parameters.sigma;
  threshold = nonStationaryEvaParams(2).parameters.threshold;

  phandles = plotGPD3D(X, timeStamps, epsilon, sigma, threshold);
  
end

