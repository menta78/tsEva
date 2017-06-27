function [returnLevels, returnLevelsErr, returnLevelsErrFit, returnLevelsErrTransf] = tsEvaComputeReturnLevelsGPDFromAnalysisObj(nonStationaryEvaParams, returnPeriodsInYears, varargin )
args.timeIndex = -1;
args = tsEasyParseNamedArgs(varargin, args);
timeIndex = args.timeIndex;

epsilon = nonStationaryEvaParams(2).parameters.epsilon;
epsilonStdErr = nonStationaryEvaParams(2).paramErr.epsilonErr;
epsilonStdErrFit = epsilonStdErr;
epsilonStdErrTransf = 0;
percentile = nonStationaryEvaParams(2).parameters.percentile;
dtSample = nonStationaryEvaParams(2).parameters.timeDeltaYears;
nonStationary = isfield(nonStationaryEvaParams(2).paramErr, 'sigmaErrTransf');
if timeIndex > 0
  sigma = nonStationaryEvaParams(2).parameters.sigma(timeIndex);
  threshold = nonStationaryEvaParams(2).parameters.threshold(timeIndex);
  sigmaStdErr = nonStationaryEvaParams(2).paramErr.sigmaErr(timeIndex);
  sigmaStdErrFit = nonStationaryEvaParams(2).paramErr.sigmaErrFit(timeIndex);
  sigmaStdErrTransf = nonStationaryEvaParams(2).paramErr.sigmaErrTransf(timeIndex);
  thresholdStdErr = nonStationaryEvaParams(2).paramErr.thresholdErr(timeIndex);
  thresholdStdErrFit = nonStationaryEvaParams(2).paramErr.thresholdErrFit(timeIndex);
  thresholdStdErrTransf = nonStationaryEvaParams(2).paramErr.thresholdErrTransf(timeIndex);
else
  sigma = nonStationaryEvaParams(2).parameters.sigma;
  threshold = nonStationaryEvaParams(2).parameters.threshold;
  sigmaStdErr = nonStationaryEvaParams(2).paramErr.sigmaErr;
  if nonStationary
    thresholdStdErr = nonStationaryEvaParams(2).paramErr.thresholdErr;
    sigmaStdErrFit = nonStationaryEvaParams(2).paramErr.sigmaErrFit;
    sigmaStdErrTransf = nonStationaryEvaParams(2).paramErr.sigmaErrTransf;
    thresholdStdErrFit = nonStationaryEvaParams(2).paramErr.thresholdErrFit;
    thresholdStdErrTransf = nonStationaryEvaParams(2).paramErr.thresholdErrTransf;
  else
    thresholdStdErr = 0;
  end
end

[returnLevels, returnLevelsErr] = tsEvaComputeReturnLevelsGPD( epsilon, sigma, threshold, percentile, epsilonStdErr, sigmaStdErr, thresholdStdErr, dtSample, returnPeriodsInYears );
if nonStationary
  [~, returnLevelsErrFit] = tsEvaComputeReturnLevelsGPD( epsilon, sigma, threshold, percentile, epsilonStdErrFit, sigmaStdErrFit, thresholdStdErrFit, dtSample, returnPeriodsInYears );
  [~, returnLevelsErrTransf] = tsEvaComputeReturnLevelsGPD( epsilon, sigma, threshold, percentile, epsilonStdErrTransf, sigmaStdErrTransf, thresholdStdErrTransf, dtSample, returnPeriodsInYears );
else
  returnLevelsErrFit = returnLevelsErr;
  returnLevelsErrTransf = zeros(size(returnLevelsErr));
end
