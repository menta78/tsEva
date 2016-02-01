function phandles = tsEvaPlotReturnLevelsGEVFromAnalysisObj( nonStationaryEvaParams, timeIndex, varargin )
% nonStationaryEvaParams, stationaryTransformData are as these returned by function nonStationaryEvaJRCApproach
% timeIndex is the index at which the time varying analysis should be
% estimated.

epsilon = nonStationaryEvaParams(1).parameters.epsilon;
sigma = mean(nonStationaryEvaParams(1).parameters.sigma(timeIndex));
mu = mean(nonStationaryEvaParams(1).parameters.mu(timeIndex));
epsilonStdErr = nonStationaryEvaParams(1).paramErr.epsilonErr;
sigmaStdErr = mean(nonStationaryEvaParams(1).paramErr.sigmaErr(timeIndex));
muStdErr = mean(nonStationaryEvaParams(1).paramErr.muErr(timeIndex));

phandles = tsEvaPlotReturnLevelsGEV(epsilon, sigma, mu, epsilonStdErr, sigmaStdErr, muStdErr, varargin{:});
end

