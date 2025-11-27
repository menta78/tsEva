% An example script related to case study 02 of Bahmanpour, et al., 2025
% (under review)

% The second case study evaluated the spatial relationship of SWH across
% three locations scattered around the Marshall Islands. Wave data
% comprised 3-hourly SWH records from a high-resolution global wave model
% (Mentaschi et al., 2023) with nearshore resolutions of 2–4 km, covering
% the period 1950–2020. This trivariate analysis highlighted spatial
% dependencies, employing a non-stationary Gaussian copula with
% non-stationary margins, modeled with GPD. Each variable was sampled at
% the 99th percentile, with univariate peaks spaced a minimum of 12 hours
% apart and a maximum allowable distance of 12 hours for multivariate
% peaks.


%REFERENCES

% [1] Bahmanpour, M.H., Mentaschi, L., Tilloy, A., Vousdoukas, M.,
%     Federico, I., Coppini, G., and Feyen, L., 2025, Transformed-Stationary
%     EVA 2.0: A Generalized Framework for Non-stationary Joint Extreme
%     Analysis (under review)
% [3] Mentaschi, L., Vousdoukas, M. I., Mentaschi, L., García-Sánchez, G.,
%     Fernández-Montblanc, T., Roland, A., Voukouvalas, E., Federico, I.,
%     Abdolali, A., Zhang, Y. J., & Feyen, L. (2023). A global unstructured,
%     coupled, high-resolution hindcast of waves and storm surge. Frontiers in
%     Marine Science, 10.
%     https://www.frontiersin.org/journals/marine-science/articles/10.3389/fmars.2023.1233679
% [4] Mentaschi, L., Vousdoukas, M. I., Voukouvalas, E., Sartini, L.,
%     Feyen, L., Besio, G., & Alfieri, L. (2016). The transformed-stationary approach:
%     a generic and simplified methodology for non-stationary extreme value
%     analysis. Hydrology and Earth System Sciences, 20(9), 3527–3547.
%     https://doi.org/10.5194/hess-20-3527-2016

% M. H. Bahmanpour, 2025

close all;
clc
clearvars;
addpath('../');

%load data
load caseStudy02_data

%percentile levels of univariate series (used for transformation)

ciPercentile = [99,99,99];     

% peak-over-threshold levels used for sampling of univariate series; it
% should be a 1-d cell array where each cell can have different size
potPercentiles=[{99},{99},{99}];  

%Non-stationary time window (in days) used for time-varying joint
%distribution
timeWindowNonStat=365*40;
%minimum distance (in days) between univariate peaks 
minDeltaUnivarSampli=[0.5,0.5,0.5];

 %maximum distance (in days) between multivariate peaks; can either take
 %one value or has to have a format and size matching
 %size(nchoosek([1:numvar],2),1) where numvar is 2 in bivariate case, 3 in
 %trivariate case, and so on
maxDeltaMultivarSampli=0.5; 

%copula family; can be gaussian or gumbel
%copulaFamily = 'Gaussian';
copulaFamily = {'Gaussian'};  

%methodology to perform univariate transformation from non-stationary to
%stationary
transfType='trendlinear';
peakType='allExceedThreshold';
[copulaAnalysis] = tsCopulaExtremes(timeAndSeries1(:,1), ...
    [timeAndSeries1(:,2),timeAndSeries2(:,2),timeAndSeries3(:,2)], ...
    'minPeakDistanceInDaysMonovarSampling',minDeltaUnivarSampli, ...
    'maxPeakDistanceInDaysMultivarSampling',maxDeltaMultivarSampli, ...
    'copulaFamily',copulaFamily,...
    'transfType',transfType,'timeWindow',timeWindowNonStat,...
    'ciPercentile',ciPercentile,'potPercentiles',potPercentiles,...
    'peakType',peakType);

[monteCarloAnalysis1] = tsCopulaMontecarlo(copulaAnalysis,...
    'nResample',10000,'timeIndex','middle'); % large montecarlo good for statistics
[monteCarloAnalysis2] = tsCopulaMontecarlo(copulaAnalysis,...
    'nResample',300,'timeIndex','middle'); % smaller montecarlo good for plotting

[gofStatistics] = tsCopulaGOFNonStat(copulaAnalysis, monteCarloAnalysis1, 'smoothInd',10);

axxArray = tsCopulaPlotTrivariate(copulaAnalysis, monteCarloAnalysis2, ...
    'gofStatistics', gofStatistics, ...
    'varLabels', {'{Loc 1}_{SWH (m)}','{Loc 2}_{SWH (m)}','{Loc 3}_{SWH (m)}'});




 
  

   
