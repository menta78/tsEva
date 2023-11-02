close all;
clc
clearvars;
addpath('../');



load('testCopulaCompoundGPDTrivariate')  %data is storm surge at three locations: Adriatic-North, Adriatic-South, and Carcavelos with identical timestamps



thresholdPercentiles=[97.5,97.5,97.5]; %threshold levels in each series for sampling of data
ciPercentile = 97.5;  %to be used in tsEvaNonStationary
potPercentiles=[97.5]; %o be used in tsEvaNonStationary; better be set to only one value but a range is also possible

timeWindow = 365.25*15; %just for sake of consistancy, this value will not be used if trendlinear is the method of nonstationary analysis

minPeakDistanceInDaysforjointpeaks=[3,3,3]; %minimum distance between monovariate peaks of each series
maxDistanceMultivariatePeaksInDays=[7]; %can either take one value or has to have a format and size matching size(nchoosek([1:numvar],2),1) where numvar is 2 in bivariate case, 3 in trivariate case, and so on
minPeakDistanceInDaysfornonstationary=[3,3,3]; %used for sampling peaks inside tsevavnonstationary

[copulaAnalysis] = tsCopulaCompoundGPD(timeAndSeries1(:,1), ...
    [timeAndSeries1(:,2),timeAndSeries2(:,2),timeAndSeries3(:,2)], ...
    thresholdPercentiles, ...
    minPeakDistanceInDaysforjointpeaks, ...
    maxDistanceMultivariatePeaksInDays, ...
    'copulaFamily','gaussian',...
    'transfType','trendlinear','timeWindow',timeWindow,...
    'ciPercentile',ciPercentile,'potPercentiles',potPercentiles,'minPeakDistanceInDays',minPeakDistanceInDaysfornonstationary);



nResample=1000;

[resampleLevel, resampleProb, resampleRetPer] = tsCopulaCompoundGPDMontecarlo(copulaAnalysis,...
    nResample,'timeIndex',1);


yMax=copulaAnalysis.jointExtremes; %non-stationary peaks


plot3(resampleLevel(:,1),resampleLevel(:,2),resampleLevel(:,3),'.r')
hold on
plot3(yMax(:,1),yMax(:,2),yMax(:,3),'.b')

view(-33,31)
grid on
xlabel('Adr.-North')
ylabel('Adr.-South')
zlabel('Carcavelos')
title('Storm surge at three locations')