function [returnLevels, returnLevelsErr] = tsEvaComputeReturnLevelsGPDFromAnalysisObj(nonStationaryEvaParams, returnPeriodsInYears, timeIndex )

epsilon = nonStationaryEvaParams(2).parameters.epsilon;
sigma = nonStationaryEvaParams(2).parameters.sigma(timeIndex);
threshold = nonStationaryEvaParams(2).parameters.threshold(timeIndex);
epsilonStdErr = nonStationaryEvaParams(2).paramErr.epsilonErr;
sigmaStdErr = nonStationaryEvaParams(2).paramErr.sigmaErr(timeIndex);
thresholdStdErr = nonStationaryEvaParams(2).paramErr.thresholdErr(timeIndex);
percentile = nonStationaryEvaParams(2).paramErr.percentile;
dtSample = nonStationaryEvaParams(2).timeDeltaYears;

[returnLevels, returnLevelsErr] = tsEvaComputeReturnLevelsGPD( epsilon, sigma, threshold, percentile, epsilonStdErr, sigmaStdErr, thresholdStdErr, dtSample, returnPeriodsInYears );