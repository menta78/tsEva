% An example script related to case study 01 of Bahmanpour, et al., 2025
% (under review)

%This case study explored the evolving relationship between river discharge
%and significant wave height (SWH) near the coast over time. The focus was
%on the mouth of the La Liane River in France, a fast-responding river
%influenced by precipitation. Wave data comprised 3-hourly SWH records from
%a high-resolution global wave model (Mentaschi et al., 2023) with
%nearshore resolutions of 2–4 km, covering the period 1950–2020. River
%discharge data was obtained from the HERA hydrological reanalysis (Tilloy
%et la., 2025). The dataset, generated with the OS LISFLOOD model (Burek et
%al., 2013), provides high resolution (approx.. 1.5 km) simulation of river
%discharge for every river with an upstream area >100km2 across Europe.).
%The data comes at six-hourly records over the same time frame.

% The analysis used the Generalized Pareto Distribution (GPD) for
% univariate margins and a time-varying Gumbel copula to model dependence.
% Univariate non-stationary margins were treated based on the method of
% Mentaschi, et al., 2016 Univariate peaks were selected with a minimum
% separation of 30 days, while joint extremes were defined as events
% occurring within a maximum interval of 45 days, ensuring that high river
% discharge followed extreme wave events. Non-stationarity was assessed
% within a 40-year moving window, with thresholds set at the 95th
% percentile for river discharge and the 99th percentile for wave height.


%REFERENCES

% [1] Bahmanpour, M.H., Mentaschi, L., Tilloy, A., Vousdoukas, M.,
%     Federico, I., Coppini, G., and Feyen, L., 2025,
%     Transformed-Stationary EVA 2.0: A Generalized Framework for
%     Non-stationary Joint Extreme Analysis (under review)
% [2] Tilloy, A., Paprotny, D., Grimaldi, S., Gomes, G., Bianchi, A.,
%     Lange,
%     S., Beck, H., & Feyen, L. (2024). HERA: a high-resolution
%     pan-European hydrological reanalysis (1950-2020). Earth Syst. Sci.
%     Data Discuss., 2024, 1–38. https://doi.org/10.5194/essd-2024-41
% [3] Mentaschi, L., Vousdoukas, M. I., Mentaschi, L., García-Sánchez, G.,
%     Fernández-Montblanc, T., Roland, A., Voukouvalas, E., Federico, I.,
%     Abdolali, A., Zhang, Y. J., & Feyen, L. (2023). A global
%     unstructured, coupled, high-resolution hindcast of waves and storm
%     surge. Frontiers in Marine Science, 10.
%     https://www.frontiersin.org/journals/marine-science/articles/10.3389/fmars.2023.1233679
% [4] Mentaschi, L., Vousdoukas, M. I., Voukouvalas, E., Sartini, L.,
%     Feyen, L., Besio, G., & Alfieri, L. (2016). The
%     transformed-stationary approach: a generic and simplified methodology
%     for non-stationary extreme value analysis. Hydrology and Earth System
%     Sciences, 20(9), 3527–3547. https://doi.org/10.5194/hess-20-3527-2016

% M. H. Bahmanpour, 2025

close all;
clc
clearvars;
addpath('../');

load caseStudy01_data

% find the overlapping part of both data sources; define a 3-hourly time
% frame
timeCommon=(max(timeSWH(1),timeRiverDisch(1))):3/24:(min(timeSWH(end),timeRiverDisch(end)));

indexGoodData=find(~isnan(riverineDischarge));
indexGoodDataW=find(~isnan(SWH));

%interpolate river and SWH data using common time frame
riverineDischarge_=interp1(timeRiverDisch(indexGoodData),riverineDischarge(indexGoodData),timeCommon);
timeSeriesRiver=[timeCommon',riverineDischarge_'];

SWH_=interp1(timeSWH(indexGoodDataW),SWH(indexGoodDataW),timeCommon);
timeSeriesSWH=[timeCommon',SWH_'];

%percentile levels of univariate series (used for transformation)
ciPercentile = [99,99];  

% peak-over-threshold levels used for sampling of univariate series; it
% should be a 1-d cell array where each cell can have different size
potPercentiles=[{95},{99}];   %95 99                           

%Non-stationary time window (in days) used for time-varying joint
%distribution
timeWindowJointDist = 365.25*40;   
%minimum distance (in days) between univariate peaks 
minDeltaUnivarSampli=[30,30]; %30, 30

 %maximum distance (in days) between multivariate peaks; can either take
 %one value or has to have a format and size matching
 %size(nchoosek([1:numvar],2),1) where numvar is 2 in bivariate case, 3 in
 %trivariate case, and so on
maxDeltaMultivarSampli=45; %45   

%copula family; Gumbel, gaussian and Frank are possible choices
copulaFamily={'gumbel'};  

%methodology to perform univariate transformation from non-stationary to
%stationary
transfType='trendlinear';

marginalDistributions='gpd';
samplingOrder=[2,1];
%

[copulaAnalysis] = tsCopulaExtremes(timeSeriesRiver(:,1), ...
    [timeSeriesRiver(:,2),timeSeriesSWH(:,2)], ...
    'minPeakDistanceInDaysMonovarSampling',minDeltaUnivarSampli, ...
    'maxPeakDistanceInDaysMultivarSampling',maxDeltaMultivarSampli, ...
    'copulaFamily',copulaFamily,...
    'transfType',transfType,'timeWindow',timeWindowJointDist,...
    'ciPercentile',ciPercentile,'potPercentiles',potPercentiles,...
    'marginalDistributions',marginalDistributions,'samplingOrder',samplingOrder);

[monteCarloAnalysis] = tsCopulaCompoundGPDMontecarlo(copulaAnalysis,...
    'nResample',1000,'timeIndex','middle');

% append monteCarloAnalysis to copulaAnalysis 
fields = fieldnames(monteCarloAnalysis); %

for ii = 1:numel(fields)
    copulaAnalysis.(fields{ii}) = monteCarloAnalysis.(fields{ii}); 
end

[gofStatistics] = tsCopulaGOFNonStat(copulaAnalysis);



axxArray = tsCopulaPlotBivariate(copulaAnalysis,gofStatistics, ...
    'ylbl', {'River discharge (m^3s^{-1})','SWH (m)'});
    
[rpAnalysis,~]=tsCopulaComputeandPlotBivarRP(copulaAnalysis,'axxArray',axxArray);

% saveas(gcf, ['testcase01','.png'], 'png');

