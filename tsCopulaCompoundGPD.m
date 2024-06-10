function [CopulaAnalysis] = tsCopulaCompoundGPD(inputtimestamps,inputtimeseries, varargin)

%tsCopulaCompoundGPD joint distribution of non-stationary compound extremes
% [CopulaAnalysis] = tsCopulaCompoundGPD(inputtimestamps,inputtimeseries,varargin)
%                     returns a variable of type structure containing various parameters 
%                     related to the joint distribution of the non-stationary extremes 

% Copula functions supported include MATLAB built-in copula functions, i.e.,:
% "gaussian"
% "t"
% "frank"
% "clayton"
% "gumbel"

% The compound (joint) extremes are sampled using the stationarized series. Transformation of the non-stationary series to 
% stationarized series and the calculation of marginal distributions are performed using method of Mentaschi, et al., 2016 [1],
% applied on each margin separately.

% input:
%  inputtimestamps                           - 1d array with length nt, time stamps for the input
%                                              time series. must be the same for all the time series
%  inputtimeseries                           - 2d array with size [nt x n], where n is the number of variables.

% other (optional) inputs:

%  samplingThresholdPrct                     - 1d array with length n, percentile threshold to be used for each time
%                                              series. 
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
%  timeSlide                                 - a scalar parameter (in number of days) where each timewindow is slided throughout the duration of inputtimeseries
%                                              to calculate time-varying copula; default value is 365


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

% M.H.Bahmanpour, 2024

%   References:
%       [1] Mentaschi, L., Vousdoukas, M., Voukouvalas, E., Sartini, L., Feyen, L., Besio, G., and Alfieri, L.: 
%           The transformed-stationary approach: a generic and simplified methodology for non-stationary extreme value analysis,
%           Hydrol. Earth Syst. Sci., 20, 3527–3547, https://doi.org/10.5194/hess-20-3527-2016, 2016

%%%%%%%%%%%%%%%%%%%%%%

% setting the default parameters
args.copulaFamily = 'gaussian';
args.marginalDistributions = 'gpd';
args.timewindow = 100*365; %long enough that in most cases a time-invariant copula be adopted as default
args.potPercentiles = {99,99}; %for bivariate case; 
args.transfType = 'trendCiPercentile';
args.ciPercentile = [99,99];
args.timeSlide=365.25; % 1 year
args.minPeakDistanceInDaysMonovarSampling=[3,3];%minpeakdistanceindays
args.maxPeakDistanceInDaysMultivarSampling=3;%maxdistancemultivariatepeaksindays
args.peakType='anyexceedthreshold';
% parsing of input parameters, overrides if different with the default
args = tsEasyParseNamedArgs(varargin, args);
copulaFamily = args.copulaFamily;
marginalDistributions=args.marginalDistributions;
timewindow = args.timewindow;
potPercentiles = args.potPercentiles;
transfType = args.transfType;
ciPercentile = args.ciPercentile;
timeSlide=args.timeSlide;
minPeakDistanceInDaysMonovarSampling=args.minPeakDistanceInDaysMonovarSampling;%minpeakdistanceindays
maxPeakDistanceInDaysMultivarSampling=args.maxPeakDistanceInDaysMultivarSampling;%maxdistancemultivariatepeaksindays
peakType=args.peakType;
% number of monovariate time series
nSeries = size(inputtimeseries, 2);
marginalAnalysis = cell(1,nSeries);

% determine whether or not a time-varying copula should be adopted
durationSeriesInYears=(inputtimestamps(end)-inputtimestamps(1))/365;

if timewindow/365<durationSeriesInYears
    timeVaryingCopula=1;
else
    timeVaryingCopula=0;
end

% perform transformation (from non-stationary to stationary) and obtain
% marginal distribution data
samplingThresholdPrct = zeros(nSeries, 1);
for ii = 1:nSeries

    [nonStatEvaParams, statTransfData] = tsEvaNonStationary([inputtimestamps,inputtimeseries(:,ii)],timewindow,'transfType',transfType,...
        'ciPercentile',ciPercentile(ii),'potPercentiles',potPercentiles{ii},'minPeakDistanceInDays',minPeakDistanceInDaysMonovarSampling(ii));

    marginalAnalysis{ii} = {nonStatEvaParams, statTransfData};
    samplingThresholdPrct(ii) = nonStatEvaParams(2).parameters.percentile;
end


% building the stationary input time series for joint GPD sampling

statInputTimeSeries=cellfun(@(x) x{2}.stationarySeries,marginalAnalysis,'UniformOutput',0);

statInputTimeSeries=[statInputTimeSeries{:}];

% it is possible that inputtimestamps changes following application of tsEvaNonStationary

inputtimestamps=cellfun(@(x) x{2}.timeStamps,marginalAnalysis,'UniformOutput',0);
inputtimeseries=cellfun(@(x) x{2}.nonStatSeries,marginalAnalysis,'UniformOutput',0);
inputtimeseries=[inputtimeseries{:}];
inputtimestamps=[inputtimestamps{:}];
inputtimestamps=inputtimestamps(:,1);
% perform the sampling of stationarized series

[samplingAnalysis] =...
    tsCopulaSampleJointPeaksMultiVariatePruning(inputtimestamps,statInputTimeSeries, ...
    'samplingThresholdPrct',samplingThresholdPrct, ...
    'minPeakDistanceInDaysMonovarSampling',minPeakDistanceInDaysMonovarSampling, ...
    'maxPeakDistanceInDaysMultivarSampling',maxPeakDistanceInDaysMultivarSampling,...
    'marginalAnalysis',marginalAnalysis,'peakType',peakType);

% translating the joint extremes into probabilities using the monovariate
% stationary distribution
if strcmpi(peakType,'anyexceedthreshold')
    jointextremes=samplingAnalysis.jointextremes;
    jointextremes2=samplingAnalysis.jointextremes2;
    jointextremes=cat(1,jointextremes,jointextremes2);
elseif strcmpi(peakType,'allexceedthreshold')
    jointextremes=samplingAnalysis.jointextremes;
end
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

% estimating the copula
copulaParam.family = copulaFamily;
copulaParam.familyId = tsCopulaGetFamilyId(copulaFamily);
copulaParam.nSeries = nSeries;

switch timeVaryingCopula
    case false
        %apply a stationary copula
        if strcmpi(copulaFamily, 'gaussian')
            % normal copula
            rho = copulafit(copulaFamily, gpdCDFCopula);
            copulaParam.rho = rho;
        elseif strcmpi(copulaFamily, 't')
            % t copula
            [rho, nu] = copulafit(copulaFamily, gpdCDFCopula);
            copulaParam.rho = rho;
            copulaParam.nu = nu;
        elseif strcmpi(copulaFamily, 'gumbel') || strcmpi(copulaFamily, 'clayton') || strcmpi(copulaFamily, 'frank')
            % one of the archimedean copulas supported by MATLAB , where only bivariate case is supported by MATLAB copulafit function
            [cprm, cci] = copulafit(copulaFamily, gpdCDFCopula);
            copulaParam.rho = cprm;
            copulaParam.nu = cci;
        else
            error(['copulaFamily not supported: ' copulaFamily]);
        end

    case true
        %apply a non-stationary copula
        jointExtremeMonovariateProbCell={};
        monovarProbJointExtrCell={};
        inputtimestampsWindowCell={};
        IndexWindowCell={};
        timePeaksCell={};
        rhoTotal={};
        nuTotal={};
        cprmTotal={};cciTotal={};
        beginIndex=0;

        timePeaks=jointextremes(:,:,1);
        dt = tsEvaGetTimeStep(inputtimestamps);
        timeWindowIndices = round(timewindow/dt);
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
           
            inputtimestampsWindowCell=[inputtimestampsWindowCell,inputtimestampsWindow];
            % of all joint peaks, find ones that fall within the time
            % window, use this index to also select probabilities that
            % belong to this window
            [Lia,~] = ismember(timePeaks,inputtimestampsWindow);
            WindowIndex=all(Lia,2);
            timePeaksCell=[timePeaksCell,timePeaks(WindowIndex,:)];
            jointExtremeMonovariateProbWindow=gpdCDFCopula(WindowIndex,:);%jointExtremeMonovariateProb(WindowIndex,:);
            jointExtremeMonovariateProbCell=[jointExtremeMonovariateProbCell,jointExtremeMonovariateProbWindow];
            
            %keep probability of extremes (i.e., CDF in an unscaled manner)
            monovarProbJointExtrWindow=gpdCDFCopula(WindowIndex,:);
            monovarProbJointExtrCell=[monovarProbJointExtrCell,monovarProbJointExtrWindow];
            % global indexing of the window in the inputtimestamps
            [~,Locb2] = ismember(timePeaks(WindowIndex,:),inputtimestamps);

            IndexWindowCell=[IndexWindowCell,Locb2];

            %increase the beginIndex which controls the while loop
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


switch timeVaryingCopula
    case false
 
        CopulaAnalysis.copulaParam=copulaParam;
        CopulaAnalysis.jointExtremeMonovariateProb=gpdCDFCopula;
        CopulaAnalysis.marginalAnalysis=marginalAnalysis;
        CopulaAnalysis.jointExtremes=samplingAnalysis.jointExtremesNS;
        CopulaAnalysis.jointExtremeTimeStamps=jointextremes(:,:,1);
        CopulaAnalysis.jointExtremeIndices=samplingAnalysis.jointExtremeIndices;
        CopulaAnalysis.peakIndicesAll=samplingAnalysis.peakIndicesAll;
        CopulaAnalysis.stationaryThresholdSampling=samplingAnalysis.thresholdsC;
        CopulaAnalysis.thresholdPotNS=[samplingAnalysis.thresholdsNonStation{:}];
        CopulaAnalysis.methodology=marginalDistributions;
        CopulaAnalysis.timeVaryingCopula=timeVaryingCopula;
        CopulaAnalysis.jointExtremeMonovariateProbNS=gpdCDFCopula;
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
        CopulaAnalysis.peakIndicesAll=samplingAnalysis.peakIndicesAll;
        CopulaAnalysis.stationaryThresholdSampling=samplingAnalysis.thresholdsC;
        CopulaAnalysis.thresholdPotNS=[samplingAnalysis.thresholdsNonStation{:}];%[thresholdPotNS{:}];
        CopulaAnalysis.methodology=marginalDistributions;
        CopulaAnalysis.timeVaryingCopula=timeVaryingCopula;
        CopulaAnalysis.jointExtremeMonovariateProbNS=monovarProbJointExtrCell;
end



