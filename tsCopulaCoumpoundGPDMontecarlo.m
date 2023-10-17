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
% x = gpinv(1-resampleProb,k,sigma,theta)
resampleRetPer = 1./(1-resampleProb);
% mx=nPeaks/timeHorizonInYears;


% on the basis of the timestamp, find the non-stationary values of the
% thresold and scale parameter
nSeries = length(marginalAnalysis);
%   ns
for ivar = 1:nSeries
    itime =find( marginalAnalysis{ivar}{2}.timeStamps==timeStamp);% compute the time index from timeStamp
    % itime=1:size(marginalAnalysis{ivar}{2}.timeStamps,1);
    % getting the ns GPD parameters
    nonStatEvaParams = marginalAnalysis{ivar}{1};
    thrshld = nonStatEvaParams(2).parameters.threshold(itime);
    scaleParam = nonStatEvaParams(2).parameters.sigma(itime);
    shapeParam = nonStatEvaParams(2).parameters.epsilon;
    nPeaks = nonStatEvaParams(2).parameters.nPeaks;
    thStart =  nonStatEvaParams(2).parameters.timeHorizonStart;
    thEnd = nonStatEvaParams(2).parameters.timeHorizonEnd;
    timeHorizonInYears = (thEnd - thStart)/365.2425;
    mx=nPeaks/timeHorizonInYears;
    resampleRetPer = (1./(resampleProb))/mx;

    rps = resampleRetPer(:,ivar);
    %       [lvls, ~] = tsEvaComputeReturnLevelsGPD(shapeParam, scaleParam, thrshld, 0, 0, 0, nPeaks, timeHorizonInYears, rps);
    [lvls, ~] = tsEvaComputeReturnLevelsGPD(shapeParam, scaleParam, thrshld, 0, 0, 0, nPeaks, timeHorizonInYears,rps');

    resampleLevel(:,ivar) = lvls;
end


