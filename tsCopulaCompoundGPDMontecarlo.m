function [resampleLevel, resampleProb, resampleRetPer] = tsCopulaCompoundGPDMontecarlo(copulaAnalysis,...
    nResample, ...
    varargin)


args.timeIndex = 1; %default time index to assess non-stationary return levels

args = tsEasyParseNamedArgs(varargin, args);

timeIndex = args.timeIndex;
copulaParam=copulaAnalysis.copulaParam;
marginalAnalysis=copulaAnalysis.marginalAnalysis;
copulaFamily = copulaParam.family;
timeVaryingCopula=copulaAnalysis.timeVaryingCopula;
switch timeVaryingCopula
    case false
if strcmpi(copulaFamily, 'gaussian')

    resampleProb = copularnd(copulaFamily, copulaParam.rho, nResample);
elseif strcmpi(copulaFamily, 't')
    resampleProb = copularnd(copulaFamily, copulaParam.rho, copulaParam.nu, nResample);
elseif strcmpi(copulaFamily, 'gumbel') || strcmpi(copulaFamily, 'clayton') || strcmpi(copulaFamily, 'frank')
    resampleProb = copularnd(copulaFamily, copulaParam.theta, nResample);
else
    error(['copulaFamily not supported: ' copulaFamily]);
end

    case true
        resampleProb={};
if strcmpi(copulaFamily, 'gaussian')
rhoCell=copulaParam.rho;
for ij=1:size(rhoCell,2)
    resampleProbT = copularnd(copulaFamily, rhoCell{ij}, nResample);
     resampleProb=[resampleProb,resampleProbT];
end
elseif strcmpi(copulaFamily, 't')
    resampleProb = copularnd(copulaFamily, copulaParam.rho, copulaParam.nu, nResample);
elseif strcmpi(copulaFamily, 'gumbel') || strcmpi(copulaFamily, 'clayton') || strcmpi(copulaFamily, 'frank')
    resampleProb = copularnd(copulaFamily, copulaParam.theta, nResample);
else
    error(['copulaFamily not supported: ' copulaFamily]);
end


end


% on the basis of the timeIndex, find the non-stationary values of the
% thresold and scale parameter
nSeries = length(marginalAnalysis);
switch timeVaryingCopula
    case false


for ivar = 1:nSeries

    nonStatEvaParams = marginalAnalysis{ivar}{1};
    thrshld = nonStatEvaParams(2).parameters.threshold(timeIndex);
    scaleParam = nonStatEvaParams(2).parameters.sigma(timeIndex);
    shapeParam = nonStatEvaParams(2).parameters.epsilon;
    nPeaks = nonStatEvaParams(2).parameters.nPeaks;
    thStart =  nonStatEvaParams(2).parameters.timeHorizonStart;
    thEnd = nonStatEvaParams(2).parameters.timeHorizonEnd;
    timeHorizonInYears = (thEnd - thStart)/365.2425;
    mx=nPeaks/timeHorizonInYears;

    resampleRetPer = (1./(resampleProb))/mx;

    rps = resampleRetPer(:,ivar);

    [lvls, ~] = tsEvaComputeReturnLevelsGPD(shapeParam, scaleParam, thrshld, 0, 0, 0, nPeaks, timeHorizonInYears,rps');

    resampleLevel(:,ivar) = lvls;
end
    case true
        resampleLevelCell={};
        resampleRetPerCell={};
for ij=1:size(rhoCell,2)
     resampleLevel=[];
     resampleRetPer=[];
for ivar = 1:nSeries

    nonStatEvaParams = marginalAnalysis{ivar}{1};
    thrshld = nonStatEvaParams(2).parameters.threshold(timeIndex);
    scaleParam = nonStatEvaParams(2).parameters.sigma(timeIndex);
    shapeParam = nonStatEvaParams(2).parameters.epsilon;
    nPeaks = nonStatEvaParams(2).parameters.nPeaks;
    thStart =  nonStatEvaParams(2).parameters.timeHorizonStart;
    thEnd = nonStatEvaParams(2).parameters.timeHorizonEnd;
    timeHorizonInYears = (thEnd - thStart)/365.2425;
    mx=nPeaks/timeHorizonInYears;

    resampleRetPer = (1./(resampleProb{ij}))/mx;

    rps = resampleRetPer(:,ivar);

    [lvls, ~] = tsEvaComputeReturnLevelsGPD(shapeParam, scaleParam, thrshld, 0, 0, 0, nPeaks, timeHorizonInYears,rps');

    resampleLevel(:,ivar) = lvls;
end
resampleLevelCell=[resampleLevelCell,resampleLevel];
resampleRetPerCell=[resampleRetPerCell,resampleRetPer];
end
resampleLevel=resampleLevelCell;
resampleRetPer=resampleRetPerCell;
end

