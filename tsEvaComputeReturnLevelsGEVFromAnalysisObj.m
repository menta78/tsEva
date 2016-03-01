function [returnLevels, returnLevelsErr] = tsEvaComputeReturnLevelsGEVFromAnalysisObj(nonStationaryEvaParams, returnPeriodsInYears, varargin)
args.timeIndex = -1;
args = tsEasyParseNamedArgs(varargin, args);
timeIndex = args.timeIndex;

epsilon = nonStationaryEvaParams(1).parameters.epsilon;
epsilonStdErr = nonStationaryEvaParams(1).paramErr.epsilonErr;
if timeIndex > 0
  sigma = nonStationaryEvaParams(1).parameters.sigma(timeIndex);
  mu = nonStationaryEvaParams(1).parameters.mu(timeIndex);
  sigmaStdErr = nonStationaryEvaParams(1).paramErr.sigmaErr(timeIndex);
  muStdErr = nonStationaryEvaParams(1).paramErr.muErr(timeIndex);
else
  sigma = nonStationaryEvaParams(1).parameters.sigma;
  mu = nonStationaryEvaParams(1).parameters.mu;
  sigmaStdErr = nonStationaryEvaParams(1).paramErr.sigmaErr;
  muStdErr = nonStationaryEvaParams(1).paramErr.muErr;
end
  
[returnLevels, returnLevelsErr] = tsEvaComputeReturnLevelsGEV( epsilon, sigma, mu, epsilonStdErr, sigmaStdErr, muStdErr, returnPeriodsInYears, varargin{:} );