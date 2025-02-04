function [CopulaAnalysis] = tsCopulaExtremes(inputtimestamps,inputtimeseries, varargin)

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
args.timewindow = 100*365.25; %long enough that in most cases a time-invariant copula be adopted as default
args.potPercentiles = cell(1,size(inputtimeseries,2)); %with empty pot percentiles, pot selected automatically 
args.transfType = 'trendCiPercentile';
args.ciPercentile = [99,99];
args.timeSlide=365.25; % 1 year
args.minPeakDistanceInDaysMonovarSampling=[3,3];%minpeakdistanceindays
args.maxPeakDistanceInDaysMultivarSampling=3;%maxdistancemultivariatepeaksindays
args.peakType='anyexceedthreshold';
args.samplingOrder=0;
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
samplingOrder=args.samplingOrder;
% number of monovariate time series
nSeries = size(inputtimeseries, 2);

% determine whether or not a time-varying copula should be adopted
durationSeriesInYears=(inputtimestamps(end)-inputtimestamps(1))/365.25;

if timewindow/365.25<durationSeriesInYears
    timeVaryingCopula=1;
else
    timeVaryingCopula=0;
end

% perform transformation (from non-stationary to stationary) and obtain
% marginal distribution data
samplingThresholdPrct = zeros(nSeries, 1);
marginalAnalysis = cell(1,nSeries);
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

inputtimestamps=[inputtimestamps{:}];


inputtimestamps=inputtimestamps(:,1);

dt0 = tsEvaGetTimeStep(inputtimestamps);
minPeakDistance = minPeakDistanceInDaysMonovarSampling/dt0;


inputtimeseries=cellfun(@(x) x{2}.nonStatSeries,marginalAnalysis,'UniformOutput',0);
% inputtimeseries= cellfun(@(x) fillmissing(x, 'constant', 0), inputtimeseries, 'UniformOutput', false); t
threshLevel=cellfun(@(x,y) prctile(x,y),inputtimeseries,potPercentiles,'UniformOutput',0);

% peaks are sampled from non-stationary series, disregarding any knowledge
% of stationary series that's why we didn't use information on the peaks
% from marginalAnalysis because those peaks are sampled based on stationary
% series; 
[pksCell,~]=cellfun(@(x,y,z) findpeaks(x,'MinPeakHeight',y,'MinPeakDistance',z),inputtimeseries,threshLevel,num2cell(minPeakDistance),'UniformOutput',0);

inputtimeseries=[inputtimeseries{:}];

% perform the sampling of stationarized series
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
    % translating the joint extremes into probabilities using the monovariate
    % stationary distribution
    %pre-allocation
    gpdCDFCopula = nan(size(jointextremes(:,:,1)));
    aicNonStat = nan(1,nSeries);
     aicStat = nan(1,nSeries);
    for ii = 1:nSeries
        nonStatEvaParams = marginalAnalysis{ii}{1};

        shapeParam = nonStatEvaParams(2).stationaryParams.parameters(2);
        scaleParam = nonStatEvaParams(2).stationaryParams.parameters(1);
        thrshldValue = nonStatEvaParams(2).stationaryParams.parameters(3);

        gpdCdf=((cdf('gp',jointextremes(:,ii,2),shapeParam, scaleParam, thrshldValue)));
        gpdCdf(gpdCdf == 0) = 1e-7;
        gpdCDFCopula(:,ii) = gpdCdf;

        % AIC of marginals in the non-stationary model
        gpdPdf=((pdf('gp',jointextremes(:,ii,2),shapeParam, scaleParam, thrshldValue)));
        k=2;
        n=length(gpdPdf);
        ll=sum(log(gpdPdf));
        aic = -2*ll + (2*n*k)/(n-k-1);
        aicNonStat(ii)=aic;
        % AIC of marginals disregarding presence of non-stationarity
        parmhat = gpfit(pksCell{ii}-threshLevel{ii});
        gpdPdf2=((pdf('gp',pksCell{ii},parmhat(1), parmhat(2), threshLevel{ii})));
        k2=2;
        n2=length(gpdPdf2);
        ll2=sum(log(gpdPdf2));
        aic2 = -2*ll2 + (2*n2*k2)/(n2-k2-1);
        aicStat(ii)=aic2;
    end
elseif  strcmpi(marginalDistributions,'gev')
    jointextremes=[];
    jointExtremeIndices=[];
    jointExtremesNS=[];
    for ii = 1:nSeries

        [annualMax, annualMaxTimeStamp, annualMaxIndexes] = tsEvaComputeAnnualMaxima([inputtimestamps,statInputTimeSeries(:,ii)]);
         [annualMax2, ~, ~] = tsEvaComputeAnnualMaxima([inputtimestamps,inputtimeseries(:,ii)]);
        tmpA=cat(3,annualMaxTimeStamp,annualMax);
        jointExtremeIndices=[jointExtremeIndices,annualMaxIndexes];
        jointextremes=[jointextremes,tmpA];
        jointExtremesNS=[jointExtremesNS,inputtimeseries(annualMaxIndexes,ii)];
    
    end
    
     gpdCDFCopula = nan(size(jointextremes(:,:,1)));
     aicNonStat = nan(1,nSeries);
     aicStat = nan(1,nSeries);
    for ii = 1:nSeries
        nonStatEvaParams = marginalAnalysis{ii}{1};

        scaleParam = nonStatEvaParams(1).stationaryParams.parameters(2);
        shapeParam = nonStatEvaParams(1).stationaryParams.parameters(1);
        locationParam = nonStatEvaParams(1).stationaryParams.parameters(3);

        gpdCdf=((cdf('gev',jointextremes(:,ii,2),shapeParam, scaleParam, locationParam)));
         gpdCdf(gpdCdf == 0) = 1e-7;
        gpdCDFCopula(:,ii) = gpdCdf;

        gpdPdf=((pdf('gev',jointextremes(:,ii,2),shapeParam, scaleParam, locationParam)));
          k=3;
          n=length(gpdPdf);
          ll=sum(log(gpdPdf));
          aicNonStat(ii) = -2*ll + (2*n*k)/(n-k-1);
      
       
         parmhat = gevfit(jointExtremesNS(:,ii));
        gpdPdf2=((pdf('gev',jointExtremesNS(:,ii),parmhat(1), parmhat(2), parmhat(3))));
        k2=3;
        n2=length(gpdPdf2);
        ll2=sum(log(gpdPdf2));
        aic2 = -2*ll2 + (2*n2*k2)/(n2-k2-1);
        aicStat(ii)=aic2;
    end
end


% estimating the copula
copulaParam.family = copulaFamily;
% copulaParam.familyId = tsCopulaGetFamilyId(copulaFamily);
copulaParam.nSeries = nSeries;

switch timeVaryingCopula
    case false
        %apply a stationary copula
        rhoC=cell(1,length(copulaFamily));
        nuC=cell(1,length(copulaFamily));
         aicStatCopula=nan(1,length(copulaFamily));
        for iFamily=1:length(copulaFamily)
            if strcmpi(copulaFamily{iFamily}, 'gaussian')
                % normal copula
                % rho = copulafit('gaussian', gpdCDFCopula);
                rho = corr(gpdCDFCopula, 'type', 'Spearman'); % apparently works bettern than copulafit
                %rho = copulaparam('Gaussian',corr(gpdCDFCopula, 'type', 'kendall'));
                y = copulapdf('Gaussian',gpdCDFCopula,rho);
                k=nSeries*(nSeries-1)/2;
                n=length(y);
                ll=sum(log(y));
                aic = -2*ll + (2*n*k)/(n-k-1);
                aicStatCopula(iFamily)=aic;
                rhoC{iFamily}=rho;

            elseif strcmpi(copulaFamily{iFamily}, 't')
                % t copula
                [rho, nu] = copulafit('t', gpdCDFCopula);           
                rhoC{iFamily}=rho;
                nuC{iFamily}=nu;
            elseif strcmpi(copulaFamily{iFamily}, 'gumbel') || strcmpi(copulaFamily{iFamily}, 'clayton') || strcmpi(copulaFamily{iFamily}, 'frank')
                % Bi-variate Archimedean copulas supported by MATLAB 
      
                [alphaHat] = copulafit(copulaFamily{iFamily}, gpdCDFCopula);
                rhoC{iFamily}=alphaHat;
                y = copulapdf(copulaFamily{iFamily},gpdCDFCopula,alphaHat);
                k=1;
                n=length(y);
                ll=sum(log(y));
                aic = -2*ll + (2*n*k)/(n-k-1);
                aicStatCopula(iFamily)=aic;
            else
                error(['copulaFamily not supported: ' copulaFamily{iFamily}]);
            end

        end
        copulaParam.rho=rhoC;
        copulaParam.nu=nuC;
    case true
        %apply a non-stationary copula
        jointExtremeMonovariateProbCell={};
        monovarProbJointExtrCell={};
        inputtimestampsWindowCell={};
        IndexWindowCell={};
        timePeaksCell={};
        rhoTotal={};
        nuTotal={};
        rho0=cell(length(copulaFamily),1);
         nu0=cell(length(copulaFamily),1);
        cprmTotal={};cciTotal={};
        beginIndex=0;

        timePeaks=jointextremes(:,:,1);
        dt = tsEvaGetTimeStep(inputtimestamps);
        timeWindowIndices = round(timewindow/dt);
        timeSlideIndices = round(timeSlide/dt);
        gpdPDFCopulaReconst=nan(size(gpdCDFCopula,1),1); %reconstructed PDF copula with non-overlapping windows for purpose of AIC calculation
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
              copulaParam.inputtimestampsWindowCell=inputtimestampsWindowCell;
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

          
            for iFamily=1:length(copulaFamily)
                if strcmpi(copulaFamily{iFamily}, 'Gaussian')
                    % rho = copulafit('gaussian', jointExtremeMonovariateProbWindow);
                    rho = corr(jointExtremeMonovariateProbWindow, 'type', 'Spearman'); % apparently works bettern than copulafit
                    rho0{iFamily}=rho;

                    % rhoTotal=[rhoTotal,rho];
                    % copulaParam.rho = rhoTotal;
                    % copulaParam.inputtimestampsWindowCell=inputtimestampsWindowCell;

                elseif strcmpi(copulaFamily{iFamily}, 't')
                    [rho, nu] = copulafit('t', jointExtremeMonovariateProbWindow);
                    rho0{iFamily}=rho;
                    nu0{iFamily}=nu;
                    % rhoTotal=[rhoTotal,rho];
                    % nuTotal=[nuTotal,nu];
                    % copulaParam.rho = rhoTotal;
                    % copulaParam.nu = nuTotal;
                    % copulaParam.inputtimestampsWindowCell=inputtimestampsWindowCell;

                elseif strcmpi(copulaFamily{iFamily}, 'Gumbel') || strcmpi(copulaFamily{iFamily}, 'Clayton') || strcmpi(copulaFamily{iFamily}, 'Frank')
                    % one of the archimedean copulas
                    [alpha, ~] = copulafit(copulaFamily{iFamily}, jointExtremeMonovariateProbWindow);
                    rho0{iFamily}=alpha;
                    copulapdfNonStat=copulapdf(copulaFamily{iFamily}, jointExtremeMonovariateProbWindow,alpha);
                    if all(isnan(gpdPDFCopulaReconst(WindowIndex)))
                        gpdPDFCopulaReconst(WindowIndex)=copulapdfNonStat;
                    end
                    % cprmTotal=[cprmTotal,cprm];
                    % cciTotal=[cciTotal,cci];
                    % copulaParam.rho = cprmTotal;
                    % copulaParam.nu = cciTotal;
                    % copulaParam.inputtimestampsWindowCell=inputtimestampsWindowCell;

                else
                    error(['copulaFamily not supported: ' copulaFamily{iFamily}]);
                end
                
            end
            rhoTotal=[rhoTotal,rho0];
            nuTotal=[nuTotal,nu0];
        end

        ttrho=cellfun(@(x) mean(x),inputtimestampsWindowCell);
        xxrho=[rhoTotal{:}];
        xxrho2=[rhoTotal{:}];
        if strcmpi(copulaFamily{1},'gaussian')
            xxrhorshp=reshape(xxrho,nSeries.^2,size(rhoTotal,2));
            xxrho=[];
            xxrhorshp(xxrhorshp==1)=nan;
            for ix=1:size(xxrhorshp,2)
                xxrho=[xxrho,mean(xxrhorshp(:,ix),'omitmissing')];
            end
            timePeakss=mean(timePeaks,2);
            rhoInterp2=interp1(ttrho,xxrho,mean(timePeaks,2),'nearest','extrap');
            ypdf=nan(size(gpdCDFCopula,1),1);
            aicNonStatCopula=nan(1,length(copulaFamily));
            for ii=1:size(gpdCDFCopula,1)
                [~,immin]=min(abs(ttrho-timePeakss(ii)));
                crx=xxrho2(:,(immin-1)*nSeries+1:(immin-1)*nSeries+3);
                ypdf(ii) = copulapdf(copulaFamily{iFamily},gpdCDFCopula(ii,:),crx);
            end
        else
            rhoInterp2=interp1(ttrho,xxrho,mean(timePeaks,2),'nearest','extrap');
            ypdf=nan(size(gpdCDFCopula,1),1);
            aicNonStatCopula=nan(1,length(copulaFamily));
           
            for ii=1:size(gpdCDFCopula,1)
                ypdf(ii) = copulapdf(copulaFamily{iFamily},gpdCDFCopula(ii,:),rhoInterp2(ii));
            end
        end
        
        
        if strcmpi(copulaFamily{1},'gaussian')
            k=nSeries*(nSeries-1)/2;
        else

            k=1;
        end
                n=length(ypdf);
                ll=sum(log(ypdf));
                aic = -2*ll + (2*n*k)/(n-k-1);
                aicNonStatCopula(iFamily)=aic;
        copulaParam.rho=rhoTotal;
        copulaParam.nu=nuTotal;

end


switch timeVaryingCopula
    case false
 
        CopulaAnalysis.copulaParam=copulaParam;
        CopulaAnalysis.jointExtremeMonovariateProb=gpdCDFCopula;
        CopulaAnalysis.marginalAnalysis=marginalAnalysis;
        CopulaAnalysis.jointExtremeTimeStamps=jointextremes(:,:,1);
        CopulaAnalysis.jointExtremesStationary=jointextremes(:,:,2);
        if strcmpi(marginalDistributions,'gpd')
        CopulaAnalysis.jointExtremes=samplingAnalysis.jointExtremesNS;
        
        CopulaAnalysis.jointExtremeIndices=samplingAnalysis.jointExtremeIndices;
        CopulaAnalysis.peakIndicesAll=samplingAnalysis.peakIndicesAll;
        CopulaAnalysis.stationaryThresholdSampling=samplingAnalysis.thresholdsC;
        CopulaAnalysis.thresholdPotNS=[samplingAnalysis.thresholdsNonStation{:}];
        elseif strcmpi(marginalDistributions,'gev')
            CopulaAnalysis.jointExtremes=jointExtremesNS;
        
        CopulaAnalysis.jointExtremeIndices=jointExtremeIndices;
        CopulaAnalysis.peakIndicesAll=nan(1,nSeries);
        CopulaAnalysis.stationaryThresholdSampling=nan(1,nSeries);
        CopulaAnalysis.thresholdPotNS=nan(1,nSeries);
        end
        CopulaAnalysis.methodology=marginalDistributions;
        CopulaAnalysis.timeVaryingCopula=timeVaryingCopula;
        CopulaAnalysis.jointExtremeMonovariateProbNS=gpdCDFCopula;
         CopulaAnalysis.aicStat=aicStat;
        CopulaAnalysis.aicNonStat=aicNonStat;
         CopulaAnalysis.aicStatCopula=aicStatCopula;
    case true
         % if strcmpi(marginalDistributions,'gpd')
        jointExtremesNS={};
       
        for ij=1:size(timePeaksCell,2)
            Indices=IndexWindowCell{ij};
            temporarySeries=[];
            for kk=1:size(Indices,2)
                temporarySeries(:,kk)=inputtimeseries(Indices(:,kk),kk);

            end
            jointExtremesNS=[jointExtremesNS,temporarySeries];
          
        end

         % elseif strcmpi(marginalDistributions,'gev')
         % 
         %     jointExtremesNS=[];
         %       jointExtremesNS2=[];
         %     for ii = 1:nSeries
         % 
         % 
         %          jointExtremesNS=[jointExtremesNS2,inputtimeseries(jointExtremeIndices(:,ii),ii)];
         % 
         %     end
         % end
        CopulaAnalysis.copulaParam=copulaParam;
        CopulaAnalysis.jointExtremeMonovariateProb=jointExtremeMonovariateProbCell;
        CopulaAnalysis.marginalAnalysis=marginalAnalysis;
        CopulaAnalysis.jointExtremes=jointExtremesNS;
        CopulaAnalysis.jointExtremeTimeStamps=timePeaksCell;
        CopulaAnalysis.jointExtremeIndices=IndexWindowCell;
        if strcmpi(marginalDistributions,'gpd')
        CopulaAnalysis.peakIndicesAll=samplingAnalysis.peakIndicesAll;
        CopulaAnalysis.stationaryThresholdSampling=samplingAnalysis.thresholdsC;
        CopulaAnalysis.thresholdPotNS=[samplingAnalysis.thresholdsNonStation{:}];%[thresholdPotNS{:}];
        elseif strcmpi(marginalDistributions,'gev')
            CopulaAnalysis.peakIndicesAll=nan(1,nSeries);
        CopulaAnalysis.stationaryThresholdSampling=nan(1,nSeries);
        CopulaAnalysis.thresholdPotNS=nan(1,nSeries);
        end
        CopulaAnalysis.methodology=marginalDistributions;
        CopulaAnalysis.timeVaryingCopula=timeVaryingCopula;
        CopulaAnalysis.jointExtremeMonovariateProbNS=monovarProbJointExtrCell;
        CopulaAnalysis.aicStat=aicStat;
        CopulaAnalysis.aicNonStat=aicNonStat;
        CopulaAnalysis.aicNonStatCopula=aicNonStatCopula;
end



