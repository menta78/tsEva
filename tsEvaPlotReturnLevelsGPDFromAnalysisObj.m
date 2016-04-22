function [phandles, returnPeriods, returnLevels, retrunLevelsErrs] = tsEvaPlotReturnLevelsGPDFromAnalysisObj( nonStationaryEvaParams, timeIndex, varargin )
% nonStationaryEvaParams, stationaryTransformData are as these returned by function nonStationaryEvaJRCApproach
% timeIndex is the index at which the time varying analysis should be
% estimated.

epsilon = nonStationaryEvaParams(2).parameters.epsilon;
sigma = nonStationaryEvaParams(2).parameters.sigma(timeIndex);
threshold = nonStationaryEvaParams(2).parameters.threshold(timeIndex);
dtSampleYears = nonStationaryEvaParams(2).parameters.timeDeltaYears;
percentile = nonStationaryEvaParams(2).parameters.percentile;
epsilonStdErr = nonStationaryEvaParams(2).paramErr.epsilonErr;
sigmaStdErr = nonStationaryEvaParams(2).paramErr.sigmaErr(timeIndex);
thresholdStdErr = nonStationaryEvaParams(2).paramErr.thresholdErr(timeIndex);

[phandles, returnPeriods, returnLevels, retrunLevelsErrs] = tsEvaPlotReturnLevelsGPD...
  (epsilon, sigma, threshold, dtSampleYears, percentile, epsilonStdErr, sigmaStdErr, thresholdStdErr,...
   'dtSampleYears', dtSampleYears, varargin{:});
end

