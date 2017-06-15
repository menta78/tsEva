function [returnLevels, returnLevelsErr, returnLevelsErrFit, returnLevelsErrTransf] = tsEvaComputeReturnLevelsGEVFromAnalysisObj(nonStationaryEvaParams, returnPeriodsInYears, varargin)
args.timeIndex = -1;
args = tsEasyParseNamedArgs(varargin, args);
timeIndex = args.timeIndex;

epsilon = nonStationaryEvaParams(1).parameters.epsilon;
epsilonStdErr = nonStationaryEvaParams(1).paramErr.epsilonErr;
epsilonStdErrFit = epsilonStdErr;
epsilonStdErrTransf = 0;
if timeIndex > 0
  sigma = nonStationaryEvaParams(1).parameters.sigma(timeIndex);
  mu = nonStationaryEvaParams(1).parameters.mu(timeIndex);
  sigmaStdErr = nonStationaryEvaParams(1).paramErr.sigmaErr(timeIndex);
  sigmaStdErrFit = nonStationaryEvaParams(1).paramErr.sigmaErrFit(timeIndex);
  sigmaStdErrTransf = nonStationaryEvaParams(1).paramErr.sigmaErrTransf(timeIndex);
  muStdErr = nonStationaryEvaParams(1).paramErr.muErr(timeIndex);
  muStdErrFit = nonStationaryEvaParams(1).paramErr.muErrFit(timeIndex);
  muStdErrTransf = nonStationaryEvaParams(1).paramErr.muErrTransf(timeIndex);
else
  sigma = nonStationaryEvaParams(1).parameters.sigma;
  mu = nonStationaryEvaParams(1).parameters.mu;
  sigmaStdErrFit = nonStationaryEvaParams(1).paramErr.sigmaErrFit;
  sigmaStdErrTransf = nonStationaryEvaParams(1).paramErr.sigmaErrTransf;
  sigmaStdErr = nonStationaryEvaParams(1).paramErr.sigmaErr;
  muStdErr = nonStationaryEvaParams(1).paramErr.muErr;
  muStdErrFit = nonStationaryEvaParams(1).paramErr.muErrFit;
  muStdErrTransf = nonStationaryEvaParams(1).paramErr.muErrTransf;
end
  
[returnLevels, returnLevelsErr] = tsEvaComputeReturnLevelsGEV( epsilon, sigma, mu, epsilonStdErr, sigmaStdErr, muStdErr, returnPeriodsInYears, varargin{:} );
[~, returnLevelsErrFit] = tsEvaComputeReturnLevelsGEV( epsilon, sigma, mu, epsilonStdErrFit, sigmaStdErrFit, muStdErrFit, returnPeriodsInYears, varargin{:} );
[~, returnLevelsErrTransf] = tsEvaComputeReturnLevelsGEV( epsilon, sigma, mu, epsilonStdErrTransf, sigmaStdErrTransf, muStdErrTransf, returnPeriodsInYears, varargin{:} );

