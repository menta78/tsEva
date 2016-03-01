function [returnLevels, returnLevelsErr] = tsEvaComputeReturnLevelsGPDFromAnalysisObj(nonStationaryEvaParams, returnPeriodsInYears, varargin )
args.timeIndex = -1;
args = tsEasyParseNamedArgs(varargin, args);
timeIndex = args.timeIndex;

epsilon = nonStationaryEvaParams(2).parameters.epsilon;
epsilonStdErr = nonStationaryEvaParams(2).paramErr.epsilonErr;
percentile = nonStationaryEvaParams(2).parameters.percentile;
dtSample = nonStationaryEvaParams(2).parameters.timeDeltaYears;
if timeIndex > 0
  sigma = nonStationaryEvaParams(2).parameters.sigma(timeIndex);
  threshold = nonStationaryEvaParams(2).parameters.threshold(timeIndex);
  sigmaStdErr = nonStationaryEvaParams(2).paramErr.sigmaErr(timeIndex);
  thresholdStdErr = nonStationaryEvaParams(2).paramErr.thresholdErr(timeIndex);
else
  sigma = nonStationaryEvaParams(2).parameters.sigma;
  threshold = nonStationaryEvaParams(2).parameters.threshold;
  sigmaStdErr = nonStationaryEvaParams(2).paramErr.sigmaErr;
  thresholdStdErr = nonStationaryEvaParams(2).paramErr.thresholdErr;
end

[returnLevels, returnLevelsErr] = tsEvaComputeReturnLevelsGPD( epsilon, sigma, threshold, percentile, epsilonStdErr, sigmaStdErr, thresholdStdErr, dtSample, returnPeriodsInYears );