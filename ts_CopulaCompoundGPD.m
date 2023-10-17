function [copulaParam, jointExtremeMonovariateProb,marginalAnalysis,jointextremes,peaksjointidx,peaksjointidx2,thresholdsC] = ts_CopulaCompoundGPD(inputtimestamps, ...
    inputtimeseries, ...
    thresholdpercentiles, ...
    minpeakdistanceindays, ...
    maxdistancemultivariatepeaksindays, ...
    varargin)

args.copulaFamily = "gaussian";
args.timewindow = 1;
args.potPercentiles = [95:0.5:99];
args.transfType = "trendCiPercentile";
args.ciPercentile = 97;
args.minPeakDistanceInDays = 3;

args = tsEasyParseNamedArgs(varargin, args);

transfType = args.transfType;
copulaFamily = args.copulaFamily;
timewindow = args.timewindow;
ciPercentile = args.ciPercentile;
potPercentiles = args.potPercentiles;
minPeakDistanceInDays=args.minPeakDistanceInDays;
nSeries = size(inputtimeseries, 2); % number of monovariate time series

marginalAnalysis = cell(1,nSeries);
peakIndexes={};
for ii = 1:nSeries


    % outTmStmp = inputtimestamps;

    [nonStatEvaParams, statTransfData] = tsEvaNonStationary([inputtimestamps,inputtimeseries(:,ii)],timewindow,'transfType',transfType,...
        'ciPercentile',ciPercentile,'potPercentiles',potPercentiles,'minPeakDistanceInDays',minPeakDistanceInDays(ii));



    % [ nonStatEvaParams, statTransfData ] = tsEvaReduceOutputObjSize(nonStatEvaParams,statTransfData,outTmStmp,'maxTimeStepDist', 365.25);
    % [ nonStatEvaParams, statTransfData ] = tsEvaReduceOutputObjSize(nonStatEvaParams,statTransfData,outTmStmp);

    marginalAnalysis{ii} = [{nonStatEvaParams, statTransfData}];
end
statInputTimeSeries=cellfun(@(x) x{2}.stationarySeries,marginalAnalysis,'UniformOutput',0);%x{2}.stationarySeries

statInputTimeSeries=[statInputTimeSeries{:}];
timeStampSeries=cellfun(@(x) x{2}.timeStamps,marginalAnalysis,'UniformOutput',0);%x{2}.stationarySeries
timeStampSeries=[timeStampSeries{:}];
% building the stationary input time series for joint GPD sampling

[jointextremes,jointnonpeak,thresholdsC,peaksjointidx,peaksjointidx2] =...
    tsCopulaSampleJointPeaksMultiVariatePruning(inputtimestamps,statInputTimeSeries, ...
    thresholdpercentiles, ...
    minpeakdistanceindays, ...
    maxdistancemultivariatepeaksindays);
% jointextremes=cat(1,jointextremes,jointnonpeak);


% translating the joint extremes into probabilities using the monovariate
% stationary distribution

jointExtremeMonovariateProb = nan(size(jointextremes(:,:,1)));
for ii = 1:nSeries
    nonStatEvaParams = marginalAnalysis{ii}{1};
    statTransfData = marginalAnalysis{ii}{2};
    % !!! the indices in stationaryParams.parameters must be the same as
    % used in tsEvaNonStationary !!!
    shapeParam = nonStatEvaParams(2).stationaryParams.parameters(2);
    scaleParam = nonStatEvaParams(2).stationaryParams.parameters(1);
    thrshldValue = nonStatEvaParams(2).stationaryParams.parameters(3);
    % COMPUTING THE PROBABILITY THAT X>jointExtreme.
    % this is given by 1 - F(jointExtreme), where F is the cumulative
    % distribution.
    %1-(1-((1+(shapeParam*(jointextremes(:,ii,2)-thrshldValue))/scaleParam).^(-1/shapeParam)))
    jointExtremeMonovariateProb(:,ii) = (1- cdf('gp',jointextremes(:,ii,2),shapeParam,scaleParam,thrshldValue))*(1-thresholdpercentiles(ii)/100);
    % jointExtremeMonovariateProb(:,ii) = 1 - tsGPDCumulativeDistribution(thrshldValue, shapeParam, scaleParam, jointextremes(:,ii,2));
end
xx=jointExtremeMonovariateProb;
xx(xx==1)=nan;
mxx=max(xx,[],1);
for ijx=1:size(jointExtremeMonovariateProb,2)
    jointExtremeMonovariateProb(jointExtremeMonovariateProb(:,ijx)==1,ijx)= mxx(ijx);
end
[~,iss]=sort(jointextremes(:,1,2),1,'descend');
[~,iss2]=sort(jointextremes(:,2,2),1,'descend');
plot(jointextremes(iss,1,2),jointExtremeMonovariateProb(iss,1))
hold on
plot(jointextremes(iss2,2,2),jointExtremeMonovariateProb(iss2,2))
legend('series1','series2')
ylabel('probability')
xlabel('stationarized peak values')
title('P(X)>x')
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




            
