function [CopulaAnalysis] = tsCopulaExtremes(inputtimestamps,inputtimeseries, varargin)

% tsCopulaExtremes joint distribution of non-stationary compound extremes

% [CopulaAnalysis] = tsCopulaExtremes(inputtimestamps,inputtimeseries, varargin)
% returns a variable of type structure containing various parameters
% related to the joint distribution of the non-stationary compound
% extremes

% Copula functions supported include some MATLAB built-in copula functions:

% "gaussian"
% "frank"
% "gumbel"

% The compound (joint) extremes are sampled using the stationarized series.
% Transformation of the non-stationary series to stationarized series and
% the calculation of marginal distributions are performed using method of
% Mentaschi, et al., 2016 [1], applied on each margin separately.

% input:
%  inputtimestamps                           - 1d array with length nt, time stamps for the input
%                                              time series. must be the same for all the time series
%  inputtimeseries                           - 2d array with size [nt x n], where n is the number of variables.

% other (optional) inputs:


%  minPeakDistanceInDaysMonovarSampling      - 1d array with length n, minimum time distance (in days) among peaks of the same
%                                              variable used for sampling
%  maxPeakDistanceInDaysMultivarSampling     - maximum time distance (in days) among peaks of different variables
%                                              for the peaks to be considered joint. 1d array with
%                                              length(maxPeakDistanceInDaysMultivarSampling) =size(nchoosek([1:n],2),1),
%                                              where nchoosek([1:n],2) shows the format in which maxPeakDistanceInDaysMultivarSampling
%                                              will be interpreted. Alternatively, it can be 1d array with length 1
%  copulaFamily                              - A string indicating the type of copula function to be applied
%  marginalDistributions                     - supports marginal distributions of type "gev" or "gpd"
%  timewindow                                - scalar indicating time window (in days). if timewindow is less than duration of inputtimeseries, a time-varying
%                                              copula would be adopted; in case of a transfType parameter other than
%                                              "trendlinear", this parameter is also used within  tsEvaNonStationary; type help tsEvaNonStationary
%  potPercentiles                            - A 1d cell array of length n, with each cell either taking a scalar or 1d array of variable length, indicating
%                                              peak-over-threshold percentile levels used as input for tsEvaNonStationary for sampling;
%                                              type  help tsEvaNonStationary
%  transfType                                - a string indicating transformation type (from non-stationary to stationary series)
%                                              ;type  help tsEvaNonStationary
%  ciPercentile                              - 1d array of length n, indicating percentile level used for assessing amplitude of the confidence interval
%                                              type help tsEvaNonStationary
%  samplingOrder                             - 1d array of length n, indicating the order of precedence of peaks of univariate series



% output:
%  CopulaAnalysis:                           - A variable of type structure containing:
%                                               copulaParam                         -- Parameter(s) of copula estimated from copulafit
%                                               jointExtremeMonovariateProb         -- Exceedance monovariate probility of joint extremes
%                                                                                      used for fitting copula functions
%                                               marginalAnalysis                    -- A cell array of marginal distributions data
%                                               jointExtremes                       -- Joint extremes on the original input time series
%                                               jointExtremeTimeStamps              -- Time stamps of extreme events
%                                               jointExtremeIndices                 -- Indices of the joint extremes
%                                               peakIndicesAll                      -- Indices of all extremes (joint and non-joint extremes)
%                                               stationaryThresholdSampling         -- Threshold level on stationarized series
%                                               thresholdPotNS                      -- Non-stationary threshold level
%                                               methodology                         -- Type of the marginal distribution adopted
%                                               timeVaryingCopula                   -- [1] if timevarying, [0] if stationary copula
%                                               jointExtremeMonovariateProbNS       -- Monovariate probility of joint extremes
%

% M.H.Bahmanpour, 2025

%REFERENCES

% [1] Bahmanpour, M.H., Mentaschi, L., Tilloy, A., Vousdoukas, M.,
%     Federico, I., Coppini, G., and Feyen, L., 2025,
%     Transformed-Stationary EVA 2.0: A Generalized Framework for
%     Non-stationary Joint Extreme Analysis (submitted to Hydrology and
%     Earth System Sciences; Feb 2025)
% [2] Mentaschi, L., Vousdoukas, M. I., Voukouvalas, E., Sartini, L.,
%     Feyen, L., Besio, G., & Alfieri, L. (2016). The
%     transformed-stationary approach: a generic and simplified methodology
%     for non-stationary extreme value analysis. Hydrology and Earth System
%     Sciences, 20(9), 3527–3547. https://doi.org/10.5194/hess-20-3527-2016

%%%%%%%%%%%%%%%%%%%%%%

% setting the default parameters
args.copulaFamily = 'gaussian';
args.marginalDistributions = 'gpd'; % can be gev or gpd
args.timewindow = 100*365.25; %long enough that in most cases a time-invariant copula be adopted as default
args.potPercentiles = cell(1,size(inputtimeseries,2));
args.transfType = 'trendCiPercentile';
args.ciPercentile = [99,99];
args.timeSlide=365.25; % 1 year
args.minPeakDistanceInDaysMonovarSampling=[3,3];
args.maxPeakDistanceInDaysMultivarSampling=3;
args.peakType='allexceedthreshold';
args.samplingOrder=0;
args.smoothInd=-1;
args.timeVaryingCopula = 1;
% parsing of input parameters, overrides if different with the default
args = tsEasyParseNamedArgs(varargin, args);

copulaFamily = args.copulaFamily;
marginalDistributions=args.marginalDistributions;
timewindow = args.timewindow;
potPercentiles = args.potPercentiles;
transfType = args.transfType;
ciPercentile = args.ciPercentile;
timeSlide=args.timeSlide;
minPeakDistanceInDaysMonovarSampling=args.minPeakDistanceInDaysMonovarSampling;
maxPeakDistanceInDaysMultivarSampling=args.maxPeakDistanceInDaysMultivarSampling;
peakType=args.peakType;
samplingOrder=args.samplingOrder;
smoothInd=args.smoothInd;
timeVaryingCopula = args.timeVaryingCopula;
if smoothInd == -1
    smoothInd = ceil(timewindow/365.23/4);
end

% number of monovariate time series
nSeries = size(inputtimeseries, 2);

% determine whether or not a time-varying copula should be adopted
durationSeriesInYears=(inputtimestamps(end)-inputtimestamps(1))/365.25;

if timeVaryingCopula
    copulaTimeWindow = timewindow;
else
    copulaTimeWindow = durationSeriesInYears*365.25*2;
end

% perform transformation (from non-stationary to stationary) and obtain
% marginal distribution data
samplingThresholdPrct = zeros(nSeries, 1);
marginalAnalysis = cell(1,nSeries);

for ii = 1:nSeries

    [nonStatEvaParams, statTransfData] = tsEvaNonStationary([inputtimestamps,...
        inputtimeseries(:,ii)],timewindow,'transfType',transfType,...
        'ciPercentile',ciPercentile(ii),'potPercentiles',potPercentiles{ii},...
        'minPeakDistanceInDays',minPeakDistanceInDaysMonovarSampling(ii), ...
        'evdType', args.marginalDistributions);
    marginalAnalysis{ii} = {nonStatEvaParams, statTransfData};
    samplingThresholdPrct(ii) = nonStatEvaParams(2).parameters.percentile;
end


% building the stationary input time series for joint sampling

statInputTimeSeries=cellfun(@(x) x{2}.stationarySeries,marginalAnalysis,'UniformOutput',0);

statInputTimeSeries=[statInputTimeSeries{:}];

% it is possible that inputtimestamps changes slighltly following application of
% tsEvaNonStationary; so we extract them again

inputtimestamps=cellfun(@(x) x{2}.timeStamps,marginalAnalysis,'UniformOutput',0);

inputtimestamps=[inputtimestamps{:}];

inputtimestamps=inputtimestamps(:,1);

inputtimeseries=cellfun(@(x) x{2}.nonStatSeries,marginalAnalysis,'UniformOutput',0);

inputtimeseries=[inputtimeseries{:}];

% perform sampling of joint extremes from stationarized series
if strcmpi(marginalDistributions,'gpd')
    [samplingAnalysis] =...
        tsCopulaSampleJointPeaksMultiVariatePruning(inputtimestamps,statInputTimeSeries, ...
        'samplingThresholdPrct',samplingThresholdPrct, ...
        'minPeakDistanceInDaysMonovarSampling',minPeakDistanceInDaysMonovarSampling, ...
        'maxPeakDistanceInDaysMultivarSampling',maxPeakDistanceInDaysMultivarSampling,...
        'marginalAnalysis',marginalAnalysis,'peakType',peakType,'samplingOrder',samplingOrder);


    if strcmpi(peakType,'anyexceedthreshold')
        jointextremes=samplingAnalysis.jointextremes;
        jointextremes2=samplingAnalysis.jointextremes2;
        jointextremes=cat(1,jointextremes,jointextremes2);
    elseif strcmpi(peakType,'allexceedthreshold')
        jointextremes=samplingAnalysis.jointextremes;
    end

     if nSeries>3        
        nBivarComb=nchoosek(1:nSeries,2);
        indicesCell=cell(size(nBivarComb,1),1);
        thresholdsC=samplingAnalysis.thresholdsC;
        for in=1:size(nBivarComb,1)
            thr=thresholdsC(nBivarComb(in,:));
            jointextr=jointextremes(:,nBivarComb(in,:),2);
            indices=find(jointextr(:,1)>=thr(1)&jointextr(:,2)>=thr(2));
            indicesCell{in}=indices;
        end
    end
    % translating the joint extremes into probabilities using the monovariate
    % stationary distribution
    %pre-allocation
    gpdCDFCopula = nan(size(jointextremes(:,:,1)));

    for ii = 1:nSeries
        nonStatEvaParams = marginalAnalysis{ii}{1};

        shapeParam = nonStatEvaParams(2).stationaryParams.parameters(2);
        scaleParam = nonStatEvaParams(2).stationaryParams.parameters(1);
        thrshldValue = nonStatEvaParams(2).stationaryParams.parameters(3);

        gpdCdf=((cdf('gp',jointextremes(:,ii,2),shapeParam, scaleParam, thrshldValue)));
        gpdCdf(gpdCdf == 0) = 1e-7;
        gpdCDFCopula(:,ii) = gpdCdf;

    end
elseif  strcmpi(marginalDistributions,'gev')
    jointextremes=[];
    jointExtremeIndices=[];
    jointExtremesNS=[];
    for ii = 1:nSeries

        [annualMax, annualMaxTimeStamp, annualMaxIndexes] = tsEvaComputeAnnualMaxima([inputtimestamps,statInputTimeSeries(:,ii)]);

        tmpA=cat(3,annualMaxTimeStamp,annualMax);
        jointExtremeIndices=[jointExtremeIndices,annualMaxIndexes];
        jointextremes=[jointextremes,tmpA];
        jointExtremesNS=[jointExtremesNS,inputtimeseries(annualMaxIndexes,ii)];

    end

    gpdCDFCopula = nan(size(jointextremes(:,:,1)));

    for ii = 1:nSeries
        nonStatEvaParams = marginalAnalysis{ii}{1};

        scaleParam = nonStatEvaParams(1).stationaryParams.parameters(2);
        shapeParam = nonStatEvaParams(1).stationaryParams.parameters(1);
        locationParam = nonStatEvaParams(1).stationaryParams.parameters(3);

        gpdCdf=((cdf('gev',jointextremes(:,ii,2),shapeParam, scaleParam, locationParam)));
        gpdCdf(gpdCdf == 0) = 1e-7;
        gpdCDFCopula(:,ii) = gpdCdf;

    end
end

% estimating the copula
copulaParam.family = copulaFamily;
copulaParam.nSeries = nSeries;

%apply a non-stationary copula
if nSeries<=3
monovarProbJointExtrCell={};
timeStampsByTimeWindow={};
IndexWindowCell={};
timePeaksCell={};
rhoTotal={}; % handling with a cell, because we don't know in advance what each copula needs

beginIndex=0;

timePeaks=jointextremes(:,:,1);
dt = tsEvaGetTimeStep(inputtimestamps);
timeWindowIndices = min(round(copulaTimeWindow/dt), length(inputtimestamps));
timeSlideIndices = round(timeSlide/dt);

while beginIndex+timeWindowIndices<=length(inputtimestamps)
    % select portion that falls in each window and store it in a cell array
    % for the last timeWindow, the duration is changed to cover
    % until the end of the inputtimestamps (this ensures no peak is
    % left behind)
    if beginIndex+timeSlideIndices+timeWindowIndices>length(inputtimestamps)
        inputtimestampsWindow=inputtimestamps(beginIndex+1:end,:);
    else
        inputtimestampsWindow=inputtimestamps(beginIndex+1:beginIndex+timeWindowIndices,:);
    end

    timeStampsByTimeWindow=[timeStampsByTimeWindow,inputtimestampsWindow];
    % of all joint peaks, find ones that fall within the time
    % window, use this index to also select probabilities that
    % belong to this window
    [Lia,~] = ismember(timePeaks,inputtimestampsWindow);
    WindowIndex=all(Lia,2);
    timePeaksCell=[timePeaksCell,timePeaks(WindowIndex,:)];           

    %keep probability of extremes (i.e., CDF in an unscaled manner)
    monovarProbJointExtrWindow=gpdCDFCopula(WindowIndex,:);
    monovarProbJointExtrCell=[monovarProbJointExtrCell,monovarProbJointExtrWindow];
    % global indexing of the window in the inputtimestamps
    [~,Locb2] = ismember(timePeaks(WindowIndex,:),inputtimestamps);

    IndexWindowCell=[IndexWindowCell,Locb2];

    %increase the beginIndex which controls the while loop
    beginIndex=beginIndex+timeSlideIndices;

    rho = tsCopulaFit(copulaFamily, monovarProbJointExtrWindow);

    rhoTotal=[rhoTotal,rho];

end
copulaParam.timeStampsByTimeWindow=timeStampsByTimeWindow;
copulaParam.rhoTimeStamps = linspace(timeStampsByTimeWindow{1}(1), ...
    timeStampsByTimeWindow{end}(end), ...
    length(timeStampsByTimeWindow));

inputtimeseriesC=repmat({inputtimeseries},1,size(IndexWindowCell,2));
jointExtremesNS = cellfun(@(x, y) ...
    cell2mat(arrayfun(@(k) x(y(:,k), k), 1:nSeries, 'UniformOutput', false)), ...
    inputtimeseriesC, IndexWindowCell, 'UniformOutput', false);

rhoTotalRaw=rhoTotal;

% smoothing
N = length(rhoTotal); % Number of NxN cell arrays
for iSeries1 = 1:nSeries
    for iSeries2 = iSeries1+1:nSeries
        cmpPrm = ones(nSeries);
        comp = zeros(1, N);
        for it = 1:N
            comp(it) = rhoTotal{it}(iSeries1,iSeries2); % Extract the component
        end
        comp = smoothdata(comp,'movmean',smoothInd);
        for it = 1:N
            rhoTotal{it}(iSeries1,iSeries2) = comp(it); % Extract the component
            rhoTotal{it}(iSeries2,iSeries1) = comp(it); % Extract the component
        end

    end
end

copulaParam.rho = rhoTotal;                   
copulaParam.rhoRaw = rhoTotalRaw;
copulaParam.smoothInd = smoothInd;

cellTimePeaks=vertcat(timePeaksCell{:});
[yMax,iB,~] = unique(vertcat(jointExtremesNS{:}),'stable','rows');
tMax=cellTimePeaks(iB,:);

CopulaAnalysis.copulaParam=copulaParam;
CopulaAnalysis.marginalAnalysis=marginalAnalysis;
CopulaAnalysis.methodology=marginalDistributions;
CopulaAnalysis.timeVaryingCopula=timeVaryingCopula;

CopulaAnalysis.jointExtremes=jointExtremesNS;
CopulaAnalysis.jointExtremeTimeStamps=timePeaksCell;
CopulaAnalysis.jointExtremeIndices=IndexWindowCell;
CopulaAnalysis.jointExtremeMonovariateProbNS=monovarProbJointExtrCell;
CopulaAnalysis.yMax=yMax;
CopulaAnalysis.tMax=tMax;
CopulaAnalysis.timeWindow=timewindow;
if strcmpi(marginalDistributions,'gpd')
    CopulaAnalysis.peakIndicesAll=samplingAnalysis.peakIndicesAll;
    CopulaAnalysis.stationaryThresholdSampling=samplingAnalysis.thresholdsC;
    CopulaAnalysis.thresholdPotNS=[samplingAnalysis.thresholdsNonStation{:}];
elseif strcmpi(marginalDistributions,'gev')
    CopulaAnalysis.peakIndicesAll=nan(1,nSeries);
    CopulaAnalysis.stationaryThresholdSampling=nan(1,nSeries);
    CopulaAnalysis.thresholdPotNS=nan(1,nSeries);
end
return
end
if nSeries>3

   
monovarProbJointExtrCell={};
IndexWindowCell={};
timePeaksCell={};
rhoTotal={}; % handling with a cell, because we don't know in advance what each copula needs

beginIndex=0;
for ijx=1:length(indicesCell)
timePeaks=jointextremes(indicesCell{ijx},nBivarComb(ijx,:),1);
gpdCDFCopulax=gpdCDFCopula(indicesCell{ijx},nBivarComb(ijx,:));
dt = tsEvaGetTimeStep(inputtimestamps);
timeWindowIndices = min(round(copulaTimeWindow/dt), length(inputtimestamps));
timeSlideIndices = round(timeSlide/dt);

timeStampsByTimeWindow={};

while beginIndex+timeWindowIndices<=length(inputtimestamps)
    % select portion that falls in each window and store it in a cell array
    % for the last timeWindow, the duration is changed to cover
    % until the end of the inputtimestamps (this ensures no peak is
    % left behind)
    if beginIndex+timeSlideIndices+timeWindowIndices>length(inputtimestamps)
        inputtimestampsWindow=inputtimestamps(beginIndex+1:end,:);
    else
        inputtimestampsWindow=inputtimestamps(beginIndex+1:beginIndex+timeWindowIndices,:);
    end

    timeStampsByTimeWindow=[timeStampsByTimeWindow,inputtimestampsWindow];
    % of all joint peaks, find ones that fall within the time
    % window, use this index to also select probabilities that
    % belong to this window
    [Lia,~] = ismember(timePeaks,inputtimestampsWindow);
    WindowIndex=all(Lia,2);
    timePeaksCell=[timePeaksCell,timePeaks(WindowIndex,:)];           

    %keep probability of extremes (i.e., CDF in an unscaled manner)
    monovarProbJointExtrWindow=gpdCDFCopulax(WindowIndex,:);
    monovarProbJointExtrCell=[monovarProbJointExtrCell,monovarProbJointExtrWindow];
    % global indexing of the window in the inputtimestamps
    [~,Locb2] = ismember(timePeaks(WindowIndex,:),inputtimestamps);

    IndexWindowCell=[IndexWindowCell,Locb2];

    %increase the beginIndex which controls the while loop
    beginIndex=beginIndex+timeSlideIndices;
    if isempty(monovarProbJointExtrWindow) || size(monovarProbJointExtrWindow,1)<2
        rho = ones(2,2);
    else
        rho = tsCopulaFit(copulaFamily, monovarProbJointExtrWindow);
    end
    rhoTotal=[rhoTotal,rho];

end
beginIndex=0;
end

rhoTotal2=reshape(rhoTotal,length(rhoTotal)/length(nBivarComb),length(nBivarComb));
rhoTotal2=rhoTotal2';
rhoTotal3=repmat({ones(nSeries)},1,size(rhoTotal2,2));
xx=cellfun(@(x) x(find(triu(x,1))),rhoTotal2,'UniformOutput',1);
nBivarCell=repmat({nBivarComb},1,size(rhoTotal2,2));
C = mat2cell(xx, size(xx,1), ones(1,size(xx,2)));
rhoTotal = cellfun(@(x,y,z) setelem(x,y,z),rhoTotal3,nBivarCell,C,'UniformOutput', 0);

end
copulaParam.timeStampsByTimeWindow=timeStampsByTimeWindow;
copulaParam.rhoTimeStamps = linspace(timeStampsByTimeWindow{1}(1), ...
    timeStampsByTimeWindow{end}(end), ...
    length(timeStampsByTimeWindow));

nBivarCombCC=repmat(mat2cell(nBivarComb,ones(size(nBivarComb,1),1),size(nBivarComb,2)),1,size(timeStampsByTimeWindow,2));
nBivarComMatch=reshape(nBivarCombCC',1,size(nBivarCombCC,1)*size(nBivarCombCC,2));

%
inputtimeseriesC=repmat({inputtimeseries},1,size(IndexWindowCell,2));
inputtimeseriesCM=cellfun(@(x,y) x(:,y),inputtimeseriesC,nBivarComMatch,'UniformOutput',0);


jointExtremesNS = cellfun(@(x, y) ...
    cell2mat(arrayfun(@(k) x(y(:,k), k), 1:2, 'UniformOutput', false)), ...
    inputtimeseriesCM, IndexWindowCell, 'UniformOutput', false);

rhoTotalRaw=rhoTotal;

% smoothing
N = length(rhoTotal); % Number of NxN cell arrays
for iSeries1 = 1:nSeries
    for iSeries2 = iSeries1+1:nSeries
        cmpPrm = ones(nSeries);
        comp = zeros(1, N);
        for it = 1:N
            comp(it) = rhoTotal{it}(iSeries1,iSeries2); % Extract the component
        end
        comp = smoothdata(comp,'movmean',smoothInd);
        for it = 1:N
            rhoTotal{it}(iSeries1,iSeries2) = comp(it); % Extract the component
            rhoTotal{it}(iSeries2,iSeries1) = comp(it); % Extract the component
        end

    end
end

copulaParam.rho = rhoTotal;                   
copulaParam.rhoRaw = rhoTotalRaw;
copulaParam.smoothInd = smoothInd;

timePeaksCellReshape=reshape(timePeaksCell,size(timeStampsByTimeWindow,2),size(nBivarComb,1));
cellTimePeaks = cellfun(@(col) vertcat(col{:}), num2cell(timePeaksCellReshape, 1), 'UniformOutput', false);

jointExtremesNSReshape=reshape(jointExtremesNS,size(timeStampsByTimeWindow,2),size(nBivarComb,1));
jointExtremesNSCat = cellfun(@(col) vertcat(col{:}), num2cell(jointExtremesNSReshape, 1), 'UniformOutput', false);
[yMax,iB,~]=cellfun(@(x) unique(x,'stable','rows'),jointExtremesNSCat,'UniformOutput',0);

tMax=cellfun(@(x,y) x(y,:),cellTimePeaks,iB,'UniformOutput',0);

CopulaAnalysis.copulaParam=copulaParam;
CopulaAnalysis.marginalAnalysis=marginalAnalysis;
CopulaAnalysis.methodology=marginalDistributions;
CopulaAnalysis.timeVaryingCopula=timeVaryingCopula;


CopulaAnalysis.jointExtremes=jointExtremesNSReshape';
CopulaAnalysis.jointExtremeTimeStamps=timePeaksCellReshape';
CopulaAnalysis.jointExtremeIndices=reshape(IndexWindowCell,size(timeStampsByTimeWindow,2),size(nBivarComb,1))';

CopulaAnalysis.jointExtremeMonovariateProbNS=reshape(monovarProbJointExtrCell,size(timeStampsByTimeWindow,2),size(nBivarComb,1))';

CopulaAnalysis.yMax=yMax;
CopulaAnalysis.tMax=tMax;
CopulaAnalysis.timeWindow=timewindow;
if strcmpi(marginalDistributions,'gpd')
    CopulaAnalysis.peakIndicesAll=samplingAnalysis.peakIndicesAll;
    CopulaAnalysis.stationaryThresholdSampling=samplingAnalysis.thresholdsC;
    CopulaAnalysis.thresholdPotNS=[samplingAnalysis.thresholdsNonStation{:}];
elseif strcmpi(marginalDistributions,'gev')
    CopulaAnalysis.peakIndicesAll=nan(1,nSeries);
    CopulaAnalysis.stationaryThresholdSampling=nan(1,nSeries);
    CopulaAnalysis.thresholdPotNS=nan(1,nSeries);
end
end


function y = setelem(x, idx, val)
lin_idx = sub2ind(size(x), idx(:,1), idx(:,2));

% Extract values
x(lin_idx)=val;
   x = triu(x) + triu(x,1)'; 
    y = x;
end


