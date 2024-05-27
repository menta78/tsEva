close all;
clc
clearvars;
addpath('../');



load('testCopulaCompoundGPDTrivariate')  %data is storm surge at three locations: Adriatic-North, Adriatic-South, and Carcavelos with identical timestamps



ciPercentile = [97.5, 97.5, 97.5];  %to be used in tsEvaNonStationary
potPercentiles={99, 99, 99}; %o be used in tsEvaNonStationary; better be set to only one value but a range is also possible

timeWindow = 365.25*15; %just for sake of consistancy, this value will not be used if trendlinear is the method of nonstationary analysis

maxDistanceMultivariatePeaksInDays=[7,7,7]; %can either take one value or has to have a format and size matching size(nchoosek([1:numvar],2),1) where numvar is 2 in bivariate case, 3 in trivariate case, and so on
minPeakDistanceInDaysMonovarSampling=[3,3,3]; %used for sampling peaks inside tsevavnonstationary

[copulaAnalysis] = tsCopulaCompoundGPD(timeAndSeries1(:,1), ...
    [timeAndSeries1(:,2),timeAndSeries2(:,2),timeAndSeries3(:,2)], ...
    'minPeakDistanceInDaysMonovarSampling', minPeakDistanceInDaysMonovarSampling,... 
    'maxDistanceMultivariatePeaksInDays', maxDistanceMultivariatePeaksInDays, ...
    'copulaFamily','gaussian',...
    'transfType','trendlinear',...
    'timeWindow',timeWindow,...
    'ciPercentile',ciPercentile,...
    'potPercentiles',potPercentiles);



nResample=100;

[copulaAnalysis] = tsCopulaCompoundGPDMontecarlo(copulaAnalysis,...
    'nResample',nResample,'timeIndex',1);


yMax=copulaAnalysis.jointExtremes; %non-stationary peaks


% plotting the resampling at the begining of the time series
plot3(copulaAnalysis.resampleLevel{1}(:,1),...
      copulaAnalysis.resampleLevel{1}(:,2),...
      copulaAnalysis.resampleLevel{1}(:,3),'.r')
hold on
plot3(yMax{1}(:,1),yMax{1}(:,2),yMax{1}(:,3),'.b')

view(-33,31)
grid on
xlabel('Adr.-North')
ylabel('Adr.-South')
zlabel('Carcavelos')
title('Storm surge at three locations')