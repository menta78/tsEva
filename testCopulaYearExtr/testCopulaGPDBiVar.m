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



nResample=1000;

[resampleLevel, resampleProb, resampleRetPer] = tsCopulaCoumpoundGPDMontecarlo(copulaAnalysis,...
    nResample,'timeIndex',1);


peakIndices=copulaAnalysis.jointExtremeIndices;
yMax=[timeAndSeries1(peakIndices(:,1),2),timeAndSeries2(peakIndices(:,2),2)]; %non-stationary peaks, not same as copulaAnalysis.jointExtremes
tMax=copulaAnalysis.jointExtremeTimeStamps;

figure 

subplot(2,1,1) % the non-stationary series and threshold parameter of each 
% series is plotted to see if the peaks that were sampled from stationarized
% series also coincide with peaks in the non-stationary ones

plot(datetime(datevec(timeAndSeries1(:,1))),timeAndSeries1(:,2))
hold on
plot(datetime(datevec(timeAndSeries2(:,1))),timeAndSeries2(:,2))
marginalAnalysis=copulaAnalysis.marginalAnalysis;
thresholdC=cellfun(@(x) x{1}(2).parameters.threshold,marginalAnalysis,'UniformOutput',0);%x{2}.stationarySeries
plot(datetime(datevec(timeAndSeries1(:,1))),thresholdC{1})
plot(datetime(datevec(timeAndSeries2(:,1))),thresholdC{2})
plot(datetime(datevec(tMax(:,1))),yMax(:,1),'hr')
plot(datetime(datevec(tMax(:,2))),yMax(:,2),'.k')
legend('Series1','Series2','Thresh-Param1','Threshold-param2','Peaks1','Peaks2')
title('Non-stationary series')


subplot(2,1,2) % the stationarized series and jointExtremes are plotted to
% check if the sampling of peaks was done correctly 
yMax2=copulaAnalysis.jointExtremes;
tMax2=copulaAnalysis.jointExtremeTimeStamps;
stationarySeries=cellfun(@(x) x{2}.stationarySeries,marginalAnalysis,'UniformOutput',0);%x{2}.stationarySeries
thresholdsC=copulaAnalysis.thresholdSampling;
plot(datetime(datevec(timeAndSeries1(:,1))),stationarySeries{1})
hold on
plot(datetime(datevec(timeAndSeries1(:,1))),stationarySeries{2})
plot(datetime(datevec(timeAndSeries1(:,1))),thresholdsC(1).*ones(1,length(timeAndSeries1(:,1))))
plot(datetime(datevec(timeAndSeries1(:,1))),thresholdsC(2).*ones(1,length(timeAndSeries1(:,1))))
plot(datetime(datevec(tMax2(:,1))),yMax2(:,1),'hr')
plot(datetime(datevec(tMax2(:,2))),yMax2(:,2),'.k')
legend('Series1','Series2','Percentile-1','Percentile-2','Peaks-1','Peaks-2')
title('Stationary series')

% plotting of jointExtremes and simulated return levels from the copula 
figHnd = tsCopulaPeakExtrPlotSctrBivar(resampleLevel, yMax, 'xlbl', 'Marshall-north', 'ylbl', 'Marshall-south');



