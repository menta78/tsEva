function [CopulaAnalysis] = tsCopulaCompoundGPD(inputtimestamps, ...
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

for ii = 1:nSeries

    [nonStatEvaParams, statTransfData] = tsEvaNonStationary([inputtimestamps,inputtimeseries(:,ii)],timewindow,'transfType',transfType,...
        'ciPercentile',ciPercentile,'potPercentiles',potPercentiles,'minPeakDistanceInDays',minPeakDistanceInDays(ii));

    marginalAnalysis{ii} = {nonStatEvaParams, statTransfData};
end


% building the stationary input time series for joint GPD sampling

statInputTimeSeries=cellfun(@(x) x{2}.stationarySeries,marginalAnalysis,'UniformOutput',0);

statInputTimeSeries=[statInputTimeSeries{:}];

% perform the sampling of stationarized series

[jointextremes,jointnonpeak,thresholdsStationary,jointExtremeIndices,peakIndicesAll] =...
    tsCopulaSampleJointPeaksMultiVariatePruning(inputtimestamps,statInputTimeSeries, ...
    thresholdpercentiles, ...
    minpeakdistanceindays, ...
    maxdistancemultivariatepeaksindays);



% translating the joint extremes into probabilities using the monovariate
% stationary distribution

jointExtremeMonovariateProb = nan(size(jointextremes(:,:,1)));
for ii = 1:nSeries
    nonStatEvaParams = marginalAnalysis{ii}{1};

    shapeParam = nonStatEvaParams(2).stationaryParams.parameters(2);
    scaleParam = nonStatEvaParams(2).stationaryParams.parameters(1);
    thrshldValue = nonStatEvaParams(2).stationaryParams.parameters(3);

    % COMPUTING THE PROBABILITY THAT X>jointExtreme.
    % this is given by 1 - F(jointExtreme), where F is the cumulative
    % distribution.
    %  Pr⁡〖{X>x}= ξ_u [1+ξ((x-u)/σ)]^(-1⁄ξ)  〗
    %   ξ_u=k/n    see Coles, 2001, pp. 81-82


    jointExtremeMonovariateProb(:,ii) = (1- cdf('gp',jointextremes(:,ii,2),shapeParam,scaleParam,thrshldValue))*(1-thresholdpercentiles(ii)/100);
    % jointExtremeMonovariateProb(:,ii) = (1- cdf('gp',jointextremes(:,ii,2),shapeParam,scaleParam,thrshldValue))*(npeak/nsample);

end

% plot P(X) > x : the probability of exceeding x (where x is jointextremes)
[~,iss]=sort(jointextremes(:,1,2),1,'descend');
[~,iss2]=sort(jointextremes(:,2,2),1,'descend');
plot(jointextremes(iss,1,2),jointExtremeMonovariateProb(iss,1))
hold on
plot(jointextremes(iss2,2,2),jointExtremeMonovariateProb(iss2,2))
legend('Series1','Series2')
ylabel('Probability')
xlabel('Stationarized peak values')
title('P(X) > x')

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

trendSeries=cellfun(@(x) x{2}.trendSeries,marginalAnalysis,'UniformOutput',0);
stdDevSeries=cellfun(@(x) x{2}.stdDevSeries,marginalAnalysis,'UniformOutput',0);
thresholdPotNS=cellfun(@(x,y,z) y*z+x,trendSeries,stdDevSeries,mat2cell(thresholdsStationary,1,ones(1,length(thresholdsStationary))),'UniformOutput',0)

jointExtremes=jointextremes(:,:,2);
jointExtremesNS=[];
for ii=1:size(jointExtremes,2)
jointExtremesNS=[jointExtremesNS,inputtimeseries(jointExtremeIndices(:,ii),ii)];
end
CopulaAnalysis.copulaParam=copulaParam;
CopulaAnalysis.jointExtremeMonovariateProb=jointExtremeMonovariateProb;
CopulaAnalysis.marginalAnalysis=marginalAnalysis;
CopulaAnalysis.jointExtremes=jointExtremesNS;
CopulaAnalysis.jointExtremeTimeStamps=jointextremes(:,:,1);
CopulaAnalysis.jointExtremeIndices=jointExtremeIndices;
CopulaAnalysis.peakIndicesAll=peakIndicesAll;
CopulaAnalysis.stationaryThresholdSampling=thresholdsStationary;
CopulaAnalysis.thresholdPotNS=[thresholdPotNS{:}];




