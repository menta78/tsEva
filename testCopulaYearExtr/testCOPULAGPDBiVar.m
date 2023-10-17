close all;
clc
clearvars;
addpath('../');



load('testCopulaCompoundGPDbivariate')  %data is storm surge at two locations: Marshall-North, and Marshall-South - needs to have identical timestamps
% load testESL_closePts_bivariateNonStationary 


thresholdPercentiles=[99,99]; %threshold levels in each series for sampling of data
ciPercentile = 99;  %to be used in tsEvaNonStationary
potPercentiles=[99]; %o be used in tsEvaNonStationary; better be set to only one value 

timeWindow = 365.25*15; %just for sake of consistancy, this value will not be used if trendlinear is the method of nonstationary analysis

minPeakDistanceInDaysforjointpeaks=[3,3]; %minimum distance between monovariate peaks of each series
maxDistanceMultivariatePeaksInDays=[1]; %can either take one value or has to have a format and size matching size(nchoosek([1:numvar],2),1) where numvar is 2 in bivariate case, 3 in trivariate case, and so on
minPeakDistanceInDaysfornonstationary=[3,3]; %used for sampling peaks inside tsevavnonstationary

[copulaParam, jointExtremeMonovariateProb,marginalAnalysis,jointextremes,peaksjointidx,peaksjointidx2,thresholdsC] = ts_CopulaCompoundGPD(timeAndSeries1(:,1), ...
    [timeAndSeries1(:,2),timeAndSeries2(:,2)], ...
    thresholdPercentiles, ...
    minPeakDistanceInDaysforjointpeaks, ...
    maxDistanceMultivariatePeaksInDays, ...
    'copulaFamily','gaussian',...
    'transfType','trendlinear','timeWindow',timeWindow,...
    'ciPercentile',ciPercentile,'potPercentiles',potPercentiles,'minPeakDistanceInDays',minPeakDistanceInDaysfornonstationary);


nResample=1000;
timeStamp=cellfun(@(x) x{2}.timeStamps(1),marginalAnalysis,'UniformOutput',1);
timeStamp=timeStamp(1);
[resampleLevel, resampleProb, resampleRetPer] = tsCopulaCoumpoundGPDMontecarlo(timeStamp, ...
    marginalAnalysis, ...
    copulaParam, ...
    nResample);



yMax=[timeAndSeries1(peaksjointidx(:,1),2),timeAndSeries2(peaksjointidx(:,2),2)];
tMax=[timeAndSeries1(peaksjointidx(:,1),1),timeAndSeries2(peaksjointidx(:,2),1)];


trendSeries=cellfun(@(x) (x{2}.trendSeries),marginalAnalysis,'UniformOutput',0);%x{2}.stationarySeries
 stdDevSeries=cellfun(@(x) x{2}.stdDevSeries,marginalAnalysis,'UniformOutput',0);%x{2}.stationarySeries


figure

subplot(2,1,1)
plot(datetime(datevec(timeAndSeries1(:,1))),timeAndSeries1(:,2))
hold on
plot(datetime(datevec(timeAndSeries2(:,1))),timeAndSeries2(:,2))

thresholdC=cellfun(@(x) x{1}(2).parameters.threshold,marginalAnalysis,'UniformOutput',0);%x{2}.stationarySeries
plot(datetime(datevec(timeAndSeries1(:,1))),thresholdC{1})
plot(datetime(datevec(timeAndSeries2(:,1))),thresholdC{2})
plot(datetime(datevec(timeAndSeries1(:,1))),trendSeries{1})
plot(datetime(datevec(timeAndSeries2(:,1))),trendSeries{2})
plot(datetime(datevec(timeAndSeries1(:,1))),stdDevSeries{1})
plot(datetime(datevec(timeAndSeries2(:,1))),stdDevSeries{2})

plot(datetime(datevec(tMax(:,1))),yMax(:,1),'hr')
plot(datetime(datevec(tMax(:,2))),yMax(:,2),'.k')
legend('S1','S2','Location-Param1','Location-param2','Trend1','Trend2','A1','A2','P1','P2')
title('Non-stationary series')
subplot(2,1,2)
yMax2=jointextremes(:,:,2);
tMax2=jointextremes(:,:,1);
stationarySeries=cellfun(@(x) x{2}.stationarySeries,marginalAnalysis,'UniformOutput',0);%x{2}.stationarySeries
timeStampSeries=cellfun(@(x) x{2}.timeStamps,marginalAnalysis,'UniformOutput',0);%x{2}.stationarySeries
plot(datetime(datevec(timeStampSeries{1})),stationarySeries{1})
hold on
plot(datetime(datevec(timeStampSeries{2})),stationarySeries{2})
plot(datetime(datevec(timeStampSeries{1})),thresholdsC(1).*ones(1,length(timeStampSeries{1})))
plot(datetime(datevec(timeStampSeries{2})),thresholdsC(2).*ones(1,length(timeStampSeries{1})))
plot(datetime(datevec(tMax2(:,1))),yMax2(:,1),'hr')
plot(datetime(datevec(tMax2(:,2))),yMax2(:,2),'.k')
legend('S1','S2','99th percentile-1','99th percentile-2','P1','P2')
title('Stationary series')
figHnd = tsCopulaPeakExtrPlotSctrBivar(resampleLevel, yMax, 'xlbl', 'Marshall-north', 'ylbl', 'Marshall-south');



