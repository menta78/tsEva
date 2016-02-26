function [returnLevels, returnLevelsErr] = tsEvaComputeReturnLevelsGEVFromAnalysisObj(nonStationaryEvaParams, returnPeriodsInYears, timeIndex, varargin)

epsilon = nonStationaryEvaParams(1).parameters.epsilon;
sigma = nonStationaryEvaParams(1).parameters.sigma(timeIndex);
mu = nonStationaryEvaParams(1).parameters.mu(timeIndex);
epsilonStdErr = nonStationaryEvaParams(1).paramErr.epsilonErr;
sigmaStdErr = nonStationaryEvaParams(1).paramErr.sigmaErr(timeIndex);
muStdErr = nonStationaryEvaParams(1).paramErr.muErr(timeIndex);

[returnLevels, returnLevelsErr] = tsEvaComputeReturnLevelsGEV( epsilon, sigma, mu, epsilonStdErr, sigmaStdErr, muStdErr, returnPeriodsInYears, varargin{:} );