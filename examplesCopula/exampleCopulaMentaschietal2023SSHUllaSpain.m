% an example script to assess joint distribution of wave and storm surge
% based on Mentaschi, et al., 2023; dataset covers a period of 1950 - 2023
% witha temporal resolution of 3 hours 

%spatially, data is for 5 different locations: 

% 10.4636 E 43.2851 N : off Cecina, western Italy, in front of the river mouth;
% -1.2514 E 45.5389 N a point off the river Garonne in the western coasts
% of France off Les Huttes
% 1.5368 E 50.7543 N: a point in the northern parts of France in the strait
% of Dover (in the English Channel) off the port of Port of Boulogne-sur-Mer  
% and in front of Liane river mouth
% -8.8036 E 42.5819 N: a point off north western parts of spain in front of 
% the Rio Ulla river mouth; this point is inside a natural harbour
% -9.0268 E 42.4651 N : a point off the Rio Ulla river mouth but outisde the 
% harbor and very close to the illa de salvora island

% Mentaschi et al., 2023 data includes SSH,SWH,T01,MWD,MWL

% M. H. Bahmanpour, 2024

close all;
clc
clearvars;
addpath('../');


% copyfile('/Users/hadi/Downloads/remote-work-italy/tsEva_dvlp-copula_git9/Alois_collaboration/lorenzodata.mat',pwd)

load lorenzodata

timeAndSeries1=[tt,xx(:,16)]; 
timeAndSeries2=[tt,xx(:,21)];


samplingThresholdPrct=[99,99]; %percentile thresholds in each series for sampling of data
ciPercentile = [99,99];  %to be used in tsEvaNonStationary
potPercentiles=[{99},{99}]; %o be used in tsEvaNonStationary; better be set to only one value but a range is also possible

timeWindow = 365*80; %if time window (in days) is equal or larger than
% length of time series, copula will be of stationary type; a very small
% time window may result in error (e.g., less than 3 years)

 minPeakDistanceInDaysMonovarSampling=[3,3];    %minimum distance between monovariate peaks of each series
 maxPeakDistanceInDaysMultivarSampling=[14]; %can either take one value or has to have a format and size matching size(nchoosek([1:numvar],2),1) where numvar is 2 in bivariate case, 3 in trivariate case, and so on
minPeakDistanceInDaysMonovarDistribution=[3,3]; %used for sampling peaks inside tsevavnonstationary

[copulaAnalysis] = tsCopulaCompoundGPD(timeAndSeries1(:,1), ...
    [timeAndSeries1(:,2),timeAndSeries2(:,2)], ...
    'samplingThresholdPrct', samplingThresholdPrct,...
    'minPeakDistanceInDaysMonovarSampling',minPeakDistanceInDaysMonovarSampling, ...
    'maxPeakDistanceInDaysMultivarSampling',maxPeakDistanceInDaysMultivarSampling, ...
    'copulaFamily','gumbel',...
    'transfType','trendlinear','timeWindow',timeWindow,...
    'ciPercentile',ciPercentile,'potPercentiles',potPercentiles,...
    'minPeakDistanceInDaysMonovarDistribution',minPeakDistanceInDaysMonovarDistribution,...
    'peakType','allExceedThreshold');



[copulaAnalysis] = tsCopulaCompoundGPDMontecarlo(copulaAnalysis,...
    'nResample',1000,'timeIndex',213280);

[gofStatistics] = tsCopulaUncertainty(copulaAnalysis);

    figHnd = tsCopulaTimeVaryingPlot(copulaAnalysis, ...
        'xlbl', 'SSH onshore (m)', 'ylbl', 'SSH offshore (m)');
      [copulaAnalysis,figHnd2]=tsCopulaJointReturnPeriod(copulaAnalysis,'plotType','OR',...
          'xlbl', 'SSH onshore (m)', 'ylbl', 'SSH offshore (m)');
%%

