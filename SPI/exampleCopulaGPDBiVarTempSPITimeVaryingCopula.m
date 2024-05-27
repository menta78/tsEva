% an example script to assess joint distribution of coastal sea surface
% height against riverine discharge at various locations
% if timeWindow is selected larger than duration of time series, 
% a stationary (or time-invariant) copula will be used (i.e, one
% copula only). If the timeWindow is selected to be less than duration of
% time series, a time-varying copula will be plotted
% Copula could be of type Gaussian, t, Frank, Gumbel and Clayton

% M. H. Bahmanpour, 2024

close all;
clc
clearvars;
addpath('../');

% load('testCopulaCompoundGPDbivariate')  %data is SS (storm surge) from Marshall Islan in two different locations
% load /Users/hadi/Downloads/remote-work-italy/tsEva_dvlp-copula_git9/Alois_collaboration/lorenzodata
% load /Users/hadi/Downloads/remote-work-italy/tsEva_dvlp-copula_git9/Alois_collaboration/alloisdata
%index 2 is discharge
load temp
load SPI

[Ai,Bi]=ismember(TIMESPI,timeT2MMmonthly);
TIMESPI=TIMESPI(Ai);
SPEI=SPEI(Ai,:);
SPI=SPI(Ai,:);
timeT2MMmonthly=timeT2MMmonthly(Bi);
T2MMmonthly=T2MMmonthly(Bi,:);

% tn=max(TIMESPI(1),timeT2MMmonthly(1)):30:min(TIMESPI(end),timeT2MMmonthly(end));
% SPI=interp1(TIMESPI,SPI,tn);
% SPEI=interp1(TIMESPI,SPEI,tn);
% T2MMmonthly=interp1(timeT2MMmonthly,T2MMmonthly,tn);

% timeT2MMmonthly=timeT2MMmonthly(3:end-12);
% T2MMmonthly=T2MMmonthly(3:end-12,:);
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
    'nResample',1000);

[gofStatistics] = tsCopulaUncertainty(copulaAnalysis);

    figHnd = tsCopulaTimeVaryingPlot(copulaAnalysis, ...
        'xlbl', 'Temp °K', 'ylbl', '- SPEI');
      [copulaAnalysis,figHnd2]=tsCopulaJointReturnPeriod(copulaAnalysis,'plotType','OR',...
          'xlbl', 'Temp °K', 'ylbl', '- SPEI');

