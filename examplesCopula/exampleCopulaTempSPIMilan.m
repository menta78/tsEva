% an example script to assess joint distribution of temperature and SPEI
% index in Milan; hourly ECMWF data first converted to daily averages;
% original data has a temporal span of 1959 - 2023; Second variables is Standardized 
% Precipitation (Evapotranspiration) Index; longitude range of 8.75 - 9.75 and
% latitude range of 45.25 - 45.75 was used. ECMWF spatial resolution was 0.25 X 0.25 
% hence yielding 15 points in this square; for each 15 points, nearest point
% from SPEI-3 or SPI-3 was chosen. monthly SPEI-3 or SPI-3 was then converted to daily 
% values. both (SPEI or SPI) were multiplied by -1. SPEI and SPI dataset has a spatial
% resolution of 0.5 X 0.5 and a time coverage of 1959 - 2022. 

%Two dataset exist: SPIAprSep where only months of April to September were
%included and SPI were all months are included

% M. H. Bahmanpour, 2024

close all;
clc
clearvars;
addpath('../');


% copyfile('/Users/hadi/Downloads/remote-work-italy/tsEva_dvlp-copula_git9/SPI/temp.mat',pwd)
% copyfile('/Users/hadi/Downloads/remote-work-italy/tsEva_dvlp-copula_git9/SPI/SPI.mat',pwd)
% copyfile('/Users/hadi/Downloads/remote-work-italy/tsEva_dvlp-copula_git9/SPI/SPIAprSep.mat',pwd)
% copyfile('/Users/hadi/Downloads/remote-work-italy/tsEva_dvlp-copula_git9/SPI/tempAprSep.mat',pwd)
load temp
load SPI

%first match both series in terms of time duration
[Ai,Bi]=ismember(TIMESPI,timeT2MMmonthly);
TIMESPI=TIMESPI(Ai);
SPEI=SPEI(Ai,:);
SPI=SPI(Ai,:);
timeT2MMmonthly=timeT2MMmonthly(Bi);
T2MMmonthly=T2MMmonthly(Bi,:);



timeAndSeries1=[timeT2MMmonthly,T2MMmonthly(:,1)]; %1 - 8 - 10 SPEI ; 4 - 6 - 7 - 8
timeAndSeries2=[TIMESPI',SPEI(:,1)];
samplingThresholdPrct=[95,80]; %percentile thresholds in each series for sampling of data
ciPercentile = [95,95];  %to be used in tsEvaNonStationary
potPercentiles=[{95},{80}]; %o be used in tsEvaNonStationary; better be set to only one value but a range is also possible

timeWindow = 365*80; %if time window (in days) is equal or larger than
% length of time series, copula will be of stationary type; a very small
% time window may result in error (e.g., less than 3 years)

 minPeakDistanceInDaysMonovarSampling=[1,3];    %minimum distance between monovariate peaks of each series
 maxPeakDistanceInDaysMultivarSampling=[12*30]; %can either take one value or has to have a format and size matching size(nchoosek([1:numvar],2),1) where numvar is 2 in bivariate case, 3 in trivariate case, and so on
minPeakDistanceInDaysMonovarDistribution=[1,3]; %used for sampling peaks inside tsevavnonstationary

[copulaAnalysis] = tsCopulaCompoundGPD(timeAndSeries1(:,1), ...
    [timeAndSeries1(:,2),timeAndSeries2(:,2)], ...
    'samplingThresholdPrct', samplingThresholdPrct,...
    'minPeakDistanceInDaysMonovarSampling',minPeakDistanceInDaysMonovarSampling, ...
    'maxPeakDistanceInDaysMultivarSampling',maxPeakDistanceInDaysMultivarSampling, ...
    'copulaFamily','gumbel',...
    'transfType','trendlinear','timeWindow',timeWindow,...
    'ciPercentile',ciPercentile,'potPercentiles',potPercentiles,...
    'minPeakDistanceInDaysMonovarDistribution',minPeakDistanceInDaysMonovarDistribution,...
    'peakType','anyExceedThreshold');



[copulaAnalysis] = tsCopulaCompoundGPDMontecarlo(copulaAnalysis,...
    'nResample',1000,'timeIndex',23194);

[gofStatistics] = tsCopulaUncertainty(copulaAnalysis);

    figHnd = tsCopulaTimeVaryingPlot(copulaAnalysis, ...
        'xlbl', 'Temp °K', 'ylbl', '- SPEI');
      [copulaAnalysis,figHnd2]=tsCopulaJointReturnPeriod(copulaAnalysis,'plotType','OR',...
          'xlbl', 'Temp °K', 'ylbl', '- SPEI');
%%

load tempAprSep
load SPIAprSep

%first match both series in terms of time duration
[Ai,Bi]=ismember(TIMESPI,timeT2MMmonthly);
TIMESPI=TIMESPI(Ai);
SPEI=SPEI(Ai,:);
SPI=SPI(Ai,:);
timeT2MMmonthly=timeT2MMmonthly(Bi);
T2MMmonthly=T2MMmonthly(Bi,:);



timeAndSeries1=[timeT2MMmonthly,T2MMmonthly(:,1)]; %1 - 8 - 10 SPEI ; 4 - 6 - 7 - 8
timeAndSeries2=[TIMESPI',SPEI(:,1)];
samplingThresholdPrct=[95,80]; %percentile thresholds in each series for sampling of data
ciPercentile = [95,95];  %to be used in tsEvaNonStationary
potPercentiles=[{95},{80}]; %o be used in tsEvaNonStationary; better be set to only one value but a range is also possible

timeWindow = 365*80; %if time window (in days) is equal or larger than
% length of time series, copula will be of stationary type; a very small
% time window may result in error (e.g., less than 3 years)

 minPeakDistanceInDaysMonovarSampling=[1,3];    %minimum distance between monovariate peaks of each series
 maxPeakDistanceInDaysMultivarSampling=[12*30]; %can either take one value or has to have a format and size matching size(nchoosek([1:numvar],2),1) where numvar is 2 in bivariate case, 3 in trivariate case, and so on
minPeakDistanceInDaysMonovarDistribution=[1,3]; %used for sampling peaks inside tsevavnonstationary

[copulaAnalysis] = tsCopulaCompoundGPD(timeAndSeries1(:,1), ...
    [timeAndSeries1(:,2),timeAndSeries2(:,2)], ...
    'samplingThresholdPrct', samplingThresholdPrct,...
    'minPeakDistanceInDaysMonovarSampling',minPeakDistanceInDaysMonovarSampling, ...
    'maxPeakDistanceInDaysMultivarSampling',maxPeakDistanceInDaysMultivarSampling, ...
    'copulaFamily','gumbel',...
    'transfType','trendlinear','timeWindow',timeWindow,...
    'ciPercentile',ciPercentile,'potPercentiles',potPercentiles,...
    'minPeakDistanceInDaysMonovarDistribution',minPeakDistanceInDaysMonovarDistribution,...
    'peakType','anyExceedThreshold');



[copulaAnalysis] = tsCopulaCompoundGPDMontecarlo(copulaAnalysis,...
    'nResample',1000,'timeIndex',23194);

[gofStatistics] = tsCopulaUncertainty(copulaAnalysis);

    figHnd = tsCopulaTimeVaryingPlot(copulaAnalysis, ...
        'xlbl', 'Temp °K - AprSep', 'ylbl', '- SPEI - AprSep');
      [copulaAnalysis,figHnd2]=tsCopulaJointReturnPeriod(copulaAnalysis,'plotType','OR',...
          'xlbl', 'Temp °K - AprSep', 'ylbl', '- SPEI - AprSep');

