% An example script related to case study 03 of Bahmanpour, et al., 2025
% (under review)

% The third case study examined the relationship between surface
% temperature and the 6-month Standardized Precipitation-Evapotranspiration
% Index (SPEI) in a region south of Milan, Italy (9.25E, 45.25N). Hourly
% surface temperature data from the ECMWF ERA-5 dataset (1959–2023) was
% paired with monthly SPEI data (1959–2022) (Zhang, 2023). To better
% capture heatwave dynamics, the temperature data were smoothed using a
% 10-day running mean. The analysis was restricted to the period from April
% to September, which aligns with the growing season when drought impacts
% are at their peak and heatwaves present a significant hazard. This case
% study demonstrated a scenario where block-maxima sampling is a valid and
% simpler alternative to the POT method for analyzing extremes. A 35-year
% non-stationary time window was used in this analysis.


%REFERENCES

% [1] Bahmanpour, M.H., Mentaschi, L., Tilloy, A., Vousdoukas, M.,
%     Federico, I., Coppini, G., and Feyen, L., 2025,
%     Transformed-Stationary EVA 2.0: A Generalized Framework for
%     Non-stationary Joint Extreme Analysis (under review)
% [2] Mentaschi, L., Vousdoukas, M. I., Voukouvalas, E., Sartini, L.,
%     Feyen, L., Besio, G., & Alfieri, L. (2016). The
%     transformed-stationary approach: a generic and simplified methodology
%     for non-stationary extreme value analysis. Hydrology and Earth System
%     Sciences, 20(9), 3527–3547. https://doi.org/10.5194/hess-20-3527-2016
% [3] Zhang, Xuanze (2023). A dataset of monthly SPI and SPEI derived from
%     ERA5 over 1959-2022. figshare. Dataset.
%     https://doi.org/10.6084/m9.figshare.24485389.v1

% M. H. Bahmanpour, 2025


close all;
clc
clearvars;
addpath('../');

load caseStudy03_data

%percentile levels of univariate series (used for transformation)

ciPercentile = [99,99]; 

% peak-over-threshold levels used for sampling of univariate series; it
% should be a 1-d cell array where each cell can have different size
potPercentiles=[{75},{97}]; 

%Non-stationary time window (in days) used for time-varying joint
%distribution
timeWindowNonStat=365*35;

%minimum distance (in days) between univariate peaks 
 minDeltaUnivarSampli=[30,30];  
  %maximum distance (in days) between multivariate peaks; can either take
 %one value or has to have a format and size matching
 %size(nchoosek([1:numvar],2),1) where numvar is 2 in bivariate case, 3 in
 %trivariate case, and so on
 maxDeltaMultivarSampli=12*30; %12

 %copula family; Gumbel, gaussian and Frank are possible choices
copulaFamily={'gumbel'};  

%methodology to perform univariate transformation from non-stationary to
%stationary
transfType='trendlinear';
peakType='allExceedThreshold';
marginalDistributions='gev';
[copulaAnalysis] = tsCopulaExtremes(timeAndSeries1(:,1), ...
    [timeAndSeries2(:,2),timeAndSeries1(:,2)], ...
    'minPeakDistanceInDaysMonovarSampling',minDeltaUnivarSampli, ...
    'maxPeakDistanceInDaysMultivarSampling',maxDeltaMultivarSampli, ...
    'copulaFamily',copulaFamily,...
    'transfType',transfType,'timeWindow',timeWindowNonStat,...
    'ciPercentile',ciPercentile,'potPercentiles',potPercentiles,...
    'peakType',peakType, ...
    'marginalDistributions',marginalDistributions);

[monteCarloAnalysis1] = tsCopulaCompoundGPDMontecarlo(copulaAnalysis,...
    'nResample',10000,'timeIndex','middle','nonStationarity','margins');%'nonStationarity','margins'
[monteCarloAnalysis2] = tsCopulaCompoundGPDMontecarlo(copulaAnalysis,...
    'nResample',1000,'timeIndex','middle','nonStationarity','margins');%'nonStationarity','margins'

% append monteCarloAnalysis to copulaAnalysis 
fields = fieldnames(monteCarloAnalysis1); %

for ii = 1:numel(fields)
    copulaAnalysis.(fields{ii}) = monteCarloAnalysis1.(fields{ii}); 
end

[gofStatistics,iToCopula] = tsCopulaGOFNonStat(copulaAnalysis);

fields = fieldnames(monteCarloAnalysis2); %

for ii = 1:numel(fields)
    copulaAnalysis.(fields{ii}) = monteCarloAnalysis2.(fields{ii}); 
end

if length(copulaFamily)>1
    copulaAnalysis.copulaParam.family=copulaFamily{iToCopula};
    if copulaAnalysis.timeVaryingCopula==1
        copulaAnalysis.copulaParam.rho=copulaAnalysis.copulaParam.rho(iToCopula,:);
        copulaAnalysis.copulaParam.rhoMean=copulaAnalysis.copulaParam.rhoMean(iToCopula,:);
    else
        copulaAnalysis.copulaParam.rho=copulaAnalysis.copulaParam.rho(iToCopula);
         copulaAnalysis.copulaParam.rhoMean=copulaAnalysis.copulaParam.rhoMean(iToCopula,:);

    end
    
    for ii = 1:numel(fields)
        if strcmpi(fields{ii},'resamplelevel') || strcmpi(fields{ii},'resampleprob')
            if copulaAnalysis.timeVaryingCopula==1
                copulaAnalysis.(fields{ii}) = monteCarloAnalysis2.(fields{ii})(iToCopula,:);
            else
                copulaAnalysis.(fields{ii}) = monteCarloAnalysis2.(fields{ii})(iToCopula);
            end

        end
    end
end

axxArray = tsCopulaPlotBivariate(copulaAnalysis,gofStatistics, ...
    'ylbl', {'River discharge (m^3s^{-1})','SWH (m)'},'smoothInd',10);

[rpAnalysis,~]=tsCopulaComputeandPlotBivarRP(copulaAnalysis,'axxArray',axxArray);


