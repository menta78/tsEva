close all;
clc
clearvars;
addpath('../');



load('testCopulaCompoundGPDbivariate')  %data is storm surge at two locations: Marshall-North, and Marshall-South - needs to have identical timestamps



thresholdPercentiles=[99,99]; %threshold levels in each series for sampling of data
ciPercentile = 99;  %to be used in tsEvaNonStationary
potPercentiles=[99]; %o be used in tsEvaNonStationary; better be set to only one value but a range is also possible

timeWindow = 365.25*15; %just for sake of consistancy, this value will not be used if trendlinear is the method of nonstationary analysis

minPeakDistanceInDaysforjointpeaks=[3,3]; %minimum distance between monovariate peaks of each series
maxDistanceMultivariatePeaksInDays=[1]; %can either take one value or has to have a format and size matching size(nchoosek([1:numvar],2),1) where numvar is 2 in bivariate case, 3 in trivariate case, and so on
minPeakDistanceInDaysfornonstationary=[3,3]; %used for sampling peaks inside tsevavnonstationary

[copulaAnalysis] = tsCopulaCompoundGPD(timeAndSeries1(:,1), ...
    [timeAndSeries1(:,2),timeAndSeries2(:,2)], ...
    thresholdPercentiles, ...
    minPeakDistanceInDaysforjointpeaks, ...
    maxDistanceMultivariatePeaksInDays, ...
    'copulaFamily','gaussian',...
    'transfType','trendlinear','timeWindow',timeWindow,...
    'ciPercentile',ciPercentile,'potPercentiles',potPercentiles,'minPeakDistanceInDays',minPeakDistanceInDaysfornonstationary);


tsCopulaPlotJointReturnPeriod(copulaAnalysis,'marginalDistributions','gp','plotType','AND')% plot type: AND, OR

