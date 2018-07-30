function [phandles, returnPeriods, returnLevels, retrunLevelsErrs] = tsEvaPlotReturnLevelsGPDFromAnalysisObj( nonStationaryEvaParams, timeIndex, varargin )
% nonStationaryEvaParams, stationaryTransformData are as these returned by function nonStationaryEvaJRCApproach
% timeIndex is the index at which the time varying analysis should be
% estimated.

epsilon = nonStationaryEvaParams(2).parameters.epsilon;
sigma = nonStationaryEvaParams(2).parameters.sigma(timeIndex);
threshold = nonStationaryEvaParams(2).parameters.threshold(timeIndex);
thStart =  nonStationaryEvaParams(2).parameters.timeHorizonStart;
thEnd = nonStationaryEvaParams(2).parameters.timeHorizonEnd;
timeHorizonInYears = (thEnd - thStart)/365.2425;
nPeaks = nonStationaryEvaParams(2).parameters.nPeaks;

epsilonStdErr = nonStationaryEvaParams(2).paramErr.epsilonErr;
sigmaStdErr = nonStationaryEvaParams(2).paramErr.sigmaErr(timeIndex);
thresholdStdErr = nonStationaryEvaParams(2).paramErr.thresholdErr(timeIndex);

[phandles, returnPeriods, returnLevels, retrunLevelsErrs] = tsEvaPlotReturnLevelsGPD...
  (epsilon, sigma, threshold, epsilonStdErr, sigmaStdErr, thresholdStdErr,...
   nPeaks, timeHorizonInYears, varargin{:});
end

