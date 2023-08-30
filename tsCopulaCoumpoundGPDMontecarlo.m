function [resampleLevel, resampleProb, resampleRetPer] = tsCopulaCoumpoundGPDMontecarlo(timeStamp, ...
    marginalAnalysis, ...
    copulaParam, ...
    nResample, ...
    varargin)

%   args.logExtrapRetlev = true;
%   args = tsEasyParseNamedArgs(varargin, args);
%   logExtrapRetlev = args.logExtrapRetlev;

  copulaFamily = copulaParam.family;
  if strcmpi(copulaFamily, 'gaussian')
    resampleProb = copularnd(copulaFamily, copulaParam.rho, nResample);
  elseif strcmpi(copulaFamily, 't')
    resampleProb = copularnd(copulaFamily, copulaParam.rho, copulaParam.nu, nResample);
  elseif strcmpi(copulaFamily, 'gumbel') || strcmpi(copulaFamily, 'clayton') || strcmpi(copulaFamily, 'frank')
    resampleProb = copularnd(copulaFamily, copulaParam.theta, nResample);
  else
    error(['copulaFamily not supported: ' copulaFamily]);
  end

  resampleRetPer = 1./(1-resampleProb);

  % on the basis of the timestamp, find the non-stationary values of the
  % thresold and scale parameter
  nSeries = length(marginalAnalysis);
  ns
  for ivar = 1:nSeries
      itime = % compute the time index from timeStamp

      % getting the ns GPD parameters
      nonStatEvaParams = marginalAnalysis{ivar}(1);
      thrshld = nonStatEvaParams(2).parameters.threshold(itime);
      scaleParam = nonStatEvaParams(2).parameters.sigma(itime);
      shapeParam = nonStatEvaParams(2).parameters.epsilon;
      nPeaks = nonStatEvaParams(2).parameters.nPeaks;
      thStart =  nonStatEvaParams(2).parameters.timeHorizonStart;
      thEnd = nonStatEvaParams(2).parameters.timeHorizonEnd;
      timeHorizonInYears = (thEnd - thStart)/365.2425;

      rps = resampleRetPer(ivar);
      [lvls, ~] = tsEvaComputeReturnLevelsGPD(shapeParam, scaleParam, thrshld, 0, 0, 0, nPeaks, timeHorizonInYears, rps);
      resampleLevel(:,ivar) = lvls;
  end


