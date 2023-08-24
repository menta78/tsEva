function [copulaParam, jointExtremeMonovariateProb] = tsCopulaCompoundGPD(inputtimestamps, ...
    inputtimeseries, ...
    thresholdpercentiles, ...
    minpeakdistanceindays, ...
    maxdistancemultivariatepeaksindays, ...
    varargin)

args.copulaFamily = "gaussian";
args.transfType = "trendCiPercentile";
args = tsEasyParseNamedArgs(varargin, args);
transfType = args.transfType;
copulaFamily = args.copulaFamily;

nSeries = size(inputtimeseries, 2); % number of monovariate time series

marginalAnalysis = cell(nSeries);
for i = 1:nSeries
    thrshld = thresholdpercentiles(i);
    % invoking tsEvaNonStationary, GPD ONLY, must be:
    %    ciPercentile == thrshld
    %    evdType == "GPD"
    %    "potPercentiles" == thrshld

    [nonStatEvaParams, statTransfData] = tsEvaNonStationary("....", "transfType", transfType);
    [ nonStatEvaParams, statTransfData ] = tsEvaReduceOutputObjSize("...");
    marginalAnalysis{i} = [nonStatEvaParams, statTransfData];
end

% building the stationary input time series for joint GPD sampling
statInputTimeSeries = "...";
[jointextremes,thresholdsC,timestampstotal,pkstotal] =... 
    tsCopulaSampleJointPeaksMultiVariatePruning(statInputTimeSeries, ...
    thresholdpercentiles, ...
    minpeakdistanceindays, ...
    maxdistancemultivariatepeaksindays);

% translating the joint extremes into probabilities using the monovariate
% stationary distribution
jointExtremeMonovariateProb = size(jointextremes);
for i = 1:nSeries
    [nonStatEvaParams, statTransfData] = marginalAnalysis{i};
    % !!! the indices in stationaryParams.parameters must be the same as
    % used in tsEvaNonStationary !!!
    shapeParam = nonStatEvaParams(2).stationaryParams.parameters(2);
    scaleParam = nonStatEvaParams(2).stationaryParams.parameters(1);
    thrshldValue = nonStatEvaParams(2).stationaryParams.parameters(3);
    % COMPUTING THE PROBABILITY THAT X>jointExtreme.
    % this is given by 1 - F(jointExtreme), where F is the cumulative
    % distribution.
    jointExtremeMonovariateProb(i, :) = 1 - tsGPDCumulativeDistribution(thrshldValue, shapeParam, scaleParam, jointextremes(i,:));
end

% estimating the copula
copulaParam.family = copulaFamily;
copulaParam.familyId = tsCopulaGetFamilyId(copulaFamily);
copulaParam.nSeries = nSeries;
if strcmpi(copulaFamily, 'gaussian')
    % normal copula
    rho = copulafit(copulaFamily, jointExtremeMonovariateProb);
    copulaParam.rho = rho;
elseif strcmpi(copulaFamily, 't')
    % t copula
    [rho, nu] = copulafit(copulaFamily, jointExtremeMonovariateProb);
    copulaParam.rho = rho;
    copulaParam.nu = nu;
elseif strcmpi(copulaFamily, 'gumbel') || strcmpi(copulaFamily, 'clayton') || strcmpi(copulaFamily, 'frank')
    % one of the archimedean copulas
    [cprm, cci] = copulafit(copulaFamily, jointExtremeMonovariateProb);
    copulaParam.theta = cprm;
    copulaParam.cci = cci;
else
    error(['copulaFamily not supported: ' copulaFamily]);
end




            
