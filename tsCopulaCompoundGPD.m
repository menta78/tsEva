function [CopulaAnalysis] = tsCopulaCompoundGPD(inputtimestamps, ...
    inputtimeseries, ...
    thresholdpercentiles, ...
    minpeakdistanceindays, ...
    maxdistancemultivariatepeaksindays, ...
    varargin)
% This function estimates multivariate copula based on compound (joint)
% extremes sampled from the inputtimeseries in accordance with sampling
% criteria. Marginal distributions are assessed using tsEvaNonStationary.
% Copula functions supported include "gaussian" , "t", "frank", "clayton",
% and "gumbel". Also supports the calculation of a time-varying copula.
% Bahmanpour, M. H., 2023

% input data:
% - inputtimestamps: 1d array with length nt, time stamps for the input
%   time series. must be the same for all the time series
% - inputtimeseries: data of each time series. 2d array with size [nt x n],
% where n is the number of variables.
% - thresholdpercentiles: percentile threshold to be used for each time
% series. 1d array with length n, where n is the number of variables to be
%   considered.
% - minpeakdistanceindays: minimum time distance among peaks of the same
%   variable. 1d array with length n, where n is the number of variables to
%   be considered.
% - maxdistancemultivariatepeaksindays: maximum time distance among peaks
%   of different variables for the peaks to be considered joint. 1d array
%   with length size(nchoosek([1:n],2),1), where n is the number of
%   variables to be considered. nchoosek([1:n],2) shows the format in which
%   maxdistancemultivariatepeaksindays will be interpreted.

% output data:
% - CopulaAnalysis: structure containing:
%          copulaParam: parameter(s) of copula estimated from copulafit
%          jointExtremeMonovariateProb: probability of joint extremes
%          marginalAnalysis: a cell array of marginal distributions
%          jointExtremes: extremes on the original input time series
%          jointExtremeTimeStamps: time stamp of extreme events
%          jointExtremeIndices: indices of the joint extremes
%          peakIndicesAll: indices of all peaks (not just joint ones)
%          stationaryThresholdSampling: threshold level on stationarized series
%          thresholdPotNS: non-stationary threshold level
%          methodology: type of marginal distribution
%          timeVaryingCopula: 1 if timevarying, 0 if stationary copula
args.copulaFamily = "Gaussian";
args.marginalDistributions = "gpd";
args.timewindow = 30*365;
args.potPercentiles = [95:0.5:99];
args.transfType = "trendCiPercentile";
args.ciPercentile = 97;
args.minPeakDistanceInDays = 3;
args.timeSlide=365;
args = tsEasyParseNamedArgs(varargin, args);

transfType = args.transfType;
copulaFamily = args.copulaFamily;
timewindow = args.timewindow;
ciPercentile = args.ciPercentile;
potPercentiles = args.potPercentiles;
minPeakDistanceInDays=args.minPeakDistanceInDays;
timeSlide=args.timeSlide;
marginalDistributions=args.marginalDistributions;

nSeries = size(inputtimeseries, 2); % number of monovariate time series
marginalAnalysis = cell(1,nSeries);
durationSeriesInYears=(inputtimestamps(end)-inputtimestamps(1))/365;
if timewindow/365<durationSeriesInYears
    timeVaryingCopula=1;
else
    timeVaryingCopula=0;
end

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
    npeak=nonStatEvaParams(2).parameters.nPeaks;
    nSample=size(inputtimeseries, 1);
    % COMPUTING THE PROBABILITY THAT X>jointExtreme.
    % this is given by 1 - F(jointExtreme), where F is the cumulative
    % distribution.
    %  Pr⁡〖{X>x}= ξ_u [1+ξ((x-u)/σ)]^(-1⁄ξ)  〗
    %   ξ_u=k/n    see Coles, 2001, pp. 81-82


    % jointExtremeMonovariateProb(:,ii) = (1- cdf('gp',jointextremes(:,ii,2),shapeParam,scaleParam,thrshldValue))*(1-thresholdpercentiles(ii)/100);
    jointExtremeMonovariateProb(:,ii) = (1- cdf('gp',jointextremes(:,ii,2),shapeParam,scaleParam,thrshldValue))*(npeak/nSample);

end

% estimating the copula
copulaParam.family = copulaFamily;
copulaParam.familyId = tsCopulaGetFamilyId(copulaFamily);
copulaParam.nSeries = nSeries;
switch timeVaryingCopula
    case false
        if strcmpi(copulaFamily, 'Gaussian')
            % normal copula
            rho = copulafit(copulaFamily, jointExtremeMonovariateProb);
            copulaParam.rho = rho;
        elseif strcmpi(copulaFamily, 't')
            % t copula
            [rho, nu] = copulafit(copulaFamily, jointExtremeMonovariateProb);
            copulaParam.rho = rho;
            copulaParam.nu = nu;
        elseif strcmpi(copulaFamily, 'Gumbel') || strcmpi(copulaFamily, 'Clayton') || strcmpi(copulaFamily, 'Frank')
            % one of the archimedean copulas
            [cprm, cci] = copulafit(copulaFamily, jointExtremeMonovariateProb);
            copulaParam.rho = cprm;
            copulaParam.nu = cci;
        else
            error(['copulaFamily not supported: ' copulaFamily]);
        end

    case true

        timePeaks=jointextremes(:,:,1);
        rhoTotal={};
        nuTotal={};
        cprmTotal={};cciTotal={};
        dt = tsEvaGetTimeStep(inputtimestamps);
        timeWindowIndices = round(timewindow/dt);
        timeSlideIndices = round(timeSlide/dt);
        beginIndex=0;
        jointExtremeMonovariateProbCell={};
        inputtimestampsWindowCell={};
        IndexWindowCell={};
        timePeaksCell={};

        while beginIndex+timeWindowIndices<=length(inputtimestamps)

            inputtimestampsWindow=inputtimestamps(beginIndex+1:beginIndex+timeWindowIndices,:);
            [Lia,Locb] = ismember(timePeaks,inputtimestampsWindow);
            inputtimestampsWindowCell=[inputtimestampsWindowCell,inputtimestampsWindow];
            WindowIndex=all(Lia,2);
            [~,Locb2] = ismember(timePeaks(WindowIndex,:),inputtimestamps);
            timePeaksCell=[timePeaksCell,timePeaks(WindowIndex,:)];
            IndexWindowCell=[IndexWindowCell,Locb2];
            jointExtremeMonovariateProbWindow=jointExtremeMonovariateProb(WindowIndex,:);
            jointExtremeMonovariateProbCell=[jointExtremeMonovariateProbCell,jointExtremeMonovariateProbWindow];
            beginIndex=beginIndex+timeSlideIndices;

            if strcmpi(copulaFamily, 'Gaussian')
                rho = copulafit(copulaFamily, jointExtremeMonovariateProbWindow);
                rhoTotal=[rhoTotal,rho];
                copulaParam.rho = rhoTotal;
                copulaParam.inputtimestampsWindowCell=inputtimestampsWindowCell;

            elseif strcmpi(copulaFamily, 't')
                [rho, nu] = copulafit(copulaFamily, jointExtremeMonovariateProbWindow);
                rhoTotal=[rhoTotal,rho];
                nuTotal=[nuTotal,nu];
                copulaParam.rho = rhoTotal;
                copulaParam.nu = nuTotal;
                copulaParam.inputtimestampsWindowCell=inputtimestampsWindowCell;

            elseif strcmpi(copulaFamily, 'Gumbel') || strcmpi(copulaFamily, 'Clayton') || strcmpi(copulaFamily, 'Frank')
                % one of the archimedean copulas
                [cprm, cci] = copulafit(copulaFamily, jointExtremeMonovariateProbWindow);
                cprmTotal=[cprmTotal,cprm];
                cciTotal=[cciTotal,cci];
                copulaParam.rho = cprmTotal;
                copulaParam.nu = cciTotal;
                copulaParam.inputtimestampsWindowCell=inputtimestampsWindowCell;

            else
                error(['copulaFamily not supported: ' copulaFamily]);
            end
        end

end

trendSeries=cellfun(@(x) x{2}.trendSeries,marginalAnalysis,'UniformOutput',0);
stdDevSeries=cellfun(@(x) x{2}.stdDevSeries,marginalAnalysis,'UniformOutput',0);
thresholdPotNS=cellfun(@(x,y,z) y*z+x,trendSeries,stdDevSeries,mat2cell(thresholdsStationary,1,ones(1,length(thresholdsStationary))),'UniformOutput',0);

jointExtremes=jointextremes(:,:,2);
jointExtremesNS=[];

for ii=1:size(jointExtremes,2)
    jointExtremesNS=[jointExtremesNS,inputtimeseries(jointExtremeIndices(:,ii),ii)];
end

switch timeVaryingCopula
    case false

        CopulaAnalysis.copulaParam=copulaParam;
        CopulaAnalysis.jointExtremeMonovariateProb=jointExtremeMonovariateProb;
        CopulaAnalysis.marginalAnalysis=marginalAnalysis;
        CopulaAnalysis.jointExtremes=jointExtremesNS;
        timeIndex=1;
        jointExtremeMonovariateProbNS=[];
        for ijx=1:size(marginalAnalysis,2)

            nonStatEvaParams =marginalAnalysis{ijx}{1};
            % timeIndex=length(nonStatEvaParams(2).parameters.threshold);
            thrshld = nonStatEvaParams(2).parameters.threshold(timeIndex);
            scaleParam = nonStatEvaParams(2).parameters.sigma(timeIndex);
            shapeParam = nonStatEvaParams(2).parameters.epsilon;
            npeak=nonStatEvaParams(2).parameters.nPeaks;
            nSample=size(inputtimeseries, 1);
            jointExtremeMonovariateProbNS(:,ijx)=1-((1-cdf('gp',jointExtremesNS(:,ijx),shapeParam, scaleParam, thrshld))*(npeak/nSample));
    

        end
        CopulaAnalysis.jointExtremeTimeStamps=jointextremes(:,:,1);
        CopulaAnalysis.jointExtremeIndices=jointExtremeIndices;
        CopulaAnalysis.peakIndicesAll=peakIndicesAll;
        CopulaAnalysis.stationaryThresholdSampling=thresholdsStationary;
        CopulaAnalysis.thresholdPotNS=[thresholdPotNS{:}];
        CopulaAnalysis.methodology=marginalDistributions;
        CopulaAnalysis.timeVaryingCopula=timeVaryingCopula;
        CopulaAnalysis.jointExtremeMonovariateProbNS=jointExtremeMonovariateProbNS;
    case true

        jointExtremesNS={};
        for ij=1:size(timePeaksCell,2)
            Indices=IndexWindowCell{ij};
            temporarySeries=[];
            for kk=1:size(Indices,2)
                temporarySeries(:,kk)=inputtimeseries(Indices(:,kk),kk);

            end
            jointExtremesNS=[jointExtremesNS,temporarySeries];
        end

        CopulaAnalysis.copulaParam=copulaParam;
        CopulaAnalysis.jointExtremeMonovariateProb=jointExtremeMonovariateProbCell;
        CopulaAnalysis.marginalAnalysis=marginalAnalysis;
        CopulaAnalysis.jointExtremes=jointExtremesNS;
        CopulaAnalysis.jointExtremeTimeStamps=timePeaksCell;
        CopulaAnalysis.jointExtremeIndices=IndexWindowCell;
        CopulaAnalysis.peakIndicesAll=peakIndicesAll;
        CopulaAnalysis.stationaryThresholdSampling=thresholdsStationary;
        CopulaAnalysis.thresholdPotNS=[thresholdPotNS{:}];
        CopulaAnalysis.methodology=marginalDistributions;
        CopulaAnalysis.timeVaryingCopula=timeVaryingCopula;
      
end



