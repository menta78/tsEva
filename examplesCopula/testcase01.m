% an example script to assess joint distribution of wave and storm surge
% based on Mentaschi, et al., 2023; dataset covers a period of 1950 - 2023
% witha temporal resolution of 3 hours 

%spatially, data is for 5 different locations: 

% (1)       10.4636 E 43.2851 N : off Cecina, western Italy, in front of the river mouth;
% (2)       -1.2514 E 45.5389 N a point off the river Garonne in the western coasts
%            of France off Les Huttes
% (3)       1.5368 E 50.7543 N: a point in the northern parts of France in the strait
%           of Dover (in the English Channel) off the port of Port of Boulogne-sur-Mer  
%            and in front of Liane river mouth
% (4)       -8.8036 E 42.5819 N: a point off north western parts of spain in front of 
%           the Rio Ulla river mouth; this point is inside a natural harbour
% (5)       -9.0268 E 42.4651 N : a point off the Rio Ulla river mouth but outisde the 
%           harbor and very close to the illa de salvora island

% 6-hourly river discharge data covering period of 1950 - 2020
% Cecina, Garonne, Liane, and Ulla
% Mentaschi et al., 2023 data includes SSH,SWH,T01,MWD,MWL

% M. H. Bahmanpour, 2024

close all;
clc
clearvars;
addpath('../');


% copyfile('/Users/hadi/Downloads/remote-work-italy/tsEva_dvlp-copula_git9/Alois_collaboration/lorenzodata.mat',pwd)
% copyfile('/Users/hadi/Downloads/remote-work-italy/tsEva_dvlp-copula_git9/Alois_collaboration/alloisdata.mat',pwd)

load alloisdata.mat
load lorenzodata

timeRiverDisch=tt2(:,1); %all are unique
timeSSHSWH=tt(:,1);

%River 1     2    3      4     4
%SSH:  1     6    11    16    21
%SWH:  2     7    12    17    22

riverineDischarge=xx2(:,3);  %Liane river %3 for liane
SWH=xx(:,12); %Wave data  %12 for wave   %4  17   


timeCommon=(max(timeSSHSWH(1),timeRiverDisch(1))):3/24:(min(timeSSHSWH(end),timeRiverDisch(end)));

ig=find(~isnan(riverineDischarge));
ig2=find(~isnan(SWH));

timeAndSeries1x=interp1(timeRiverDisch(ig),riverineDischarge(ig),timeCommon);
timeAndSeries1=[timeCommon',timeAndSeries1x'];

timeAndSeries3x=interp1(timeSSHSWH(ig2),SWH(ig2),timeCommon);
timeAndSeries3=[timeCommon',timeAndSeries3x'];


ciPercentile = [99,99];                                  % used for transformation from non-stationary to stationary
potPercentiles=[{95},{99}];                              % percentile level used for sampling; could either be a single value, or a 1-d array (in each cell having different size)

timeWindowNonStat = 365.25*40;                           % if timeWindow (in days) >= duration of input time series,
                                                          % a stationary copula would be chosen; a very small
                                                               % time window may result in error (e.g., less than 3 years)
timeWindowStat = 365.25*80; 
minPeakDistanceInDaysMonovarSampling=[3,3];                 %minimum distance between monovariate peaks of each series
maxPeakDistanceInDaysMultivarSampling=[90];                  %can either take one value or has to have a format and size matching size(nchoosek([1:numvar],2),1) where numvar is 2 in bivariate case, 3 in trivariate case, and so on

copulaFamily={'gaussian','gumbel','clayton','frank'};           %clayton, t and gaussian are possible choices

% first assessing stationary copula
[copulaAnalysisStat] = tsCopulaExtremes(timeAndSeries1(:,1), ...
    [timeAndSeries1(:,2),timeAndSeries3(:,2)], ...
    'minPeakDistanceInDaysMonovarSampling',minPeakDistanceInDaysMonovarSampling, ...
    'maxPeakDistanceInDaysMultivarSampling',maxPeakDistanceInDaysMultivarSampling, ...
    'copulaFamily',copulaFamily,...
    'transfType','trendlinear','timeWindow',timeWindowStat,...
    'ciPercentile',ciPercentile,'potPercentiles',potPercentiles,...
    'peakType','allExceedThreshold', ...
    'marginalDistributions','gpd');
[copulaAnalysisStat] = tsCopulaCompoundGPDMontecarlo(copulaAnalysisStat,...
    'nResample',1000,'timeIndex','middle');


[gofStatistics] = tsCopulaGOF(copulaAnalysisStat);

%%
copulaFamily={'gumbel'};
[copulaAnalysisNonStat] = tsCopulaExtremes(timeAndSeries1(:,1), ...
    [timeAndSeries1(:,2),timeAndSeries3(:,2)], ...
    'minPeakDistanceInDaysMonovarSampling',minPeakDistanceInDaysMonovarSampling, ...
    'maxPeakDistanceInDaysMultivarSampling',maxPeakDistanceInDaysMultivarSampling, ...
    'copulaFamily',copulaFamily,...
    'transfType','trendlinear','timeWindow',timeWindowNonStat,...
    'ciPercentile',ciPercentile,'potPercentiles',potPercentiles,...
    'peakType','allExceedThreshold', ...
    'marginalDistributions','gpd');



[copulaAnalysisNonStat] = tsCopulaCompoundGPDMontecarlo(copulaAnalysisNonStat,...
    'nResample',1000);
figHnd = tsCopulaPlotBivariate(copulaAnalysisNonStat,copulaAnalysisStat,gofStatistics, ...
    'ylbl', {'River discharge (m^3s^{-1})','SWH (m)'});
    
 

