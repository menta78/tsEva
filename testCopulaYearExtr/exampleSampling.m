close all;
clc
clearvars;

addpath('../');

load('testESL_closePts_bivariateNonStationary.mat');  %some data to test this script


if any(timeAndSeries1(:,1)~=timeAndSeries2(:,1))
    disp('time stamps of each series are not indentical-an interpolation needs to be done')
    
end


thresholdPercentiles=[99,99];
minPeakDistanceInDays=[3,3];
maxDistanceMultivariatePeaksInDays=[7,7];
timeWindow = 365.25*15;
ciPercentile = 98;

gdidx=find(timeAndSeries1(:,1)<=datenum(2004,12,20));
timeAndSeries1=timeAndSeries1(gdidx,:);
timeAndSeries2=timeAndSeries2(gdidx,:);


[~,isort1]=sort(timeAndSeries1(:,1),'ascend');
[~,isort2]=sort(timeAndSeries2(:,1),'ascend');
timeAndSeries1=timeAndSeries1(isort1,:);
timeAndSeries2=timeAndSeries2(isort2,:);
[nonStatEvaParams1, statTransfData1] = tsEvaNonStationary(timeAndSeries1, timeWindow, 'transfType', 'trendlinear',...
    'ciPercentile', ciPercentile,'potPercentiles',[97:0.5:99], 'minPeakDistanceInDays', 14);
plot(datetime(datevec(timeAndSeries1(:,1))),statTransfData1.stationarySeries,'--')
[nonStatEvaParams2, statTransfData2] = tsEvaNonStationary(timeAndSeries2, timeWindow, 'transfType', 'trendlinear',...
    'ciPercentile', ciPercentile,'potPercentiles',[97:0.5:99], 'minPeakDistanceInDays', 14);
hold on
plot(datetime(datevec(timeAndSeries1(:,1))),statTransfData2.stationarySeries,'--')

[jointExtremes, thresholds] = tsCopulaSampleJointExtremes(timeAndSeries1(:,1), ...
    [statTransfData1.stationarySeries,statTransfData2.stationarySeries], ...
    thresholdPercentiles, ...
    minPeakDistanceInDays, ...
    maxDistanceMultivariatePeaksInDays);

plot(datetime(datevec(jointExtremes(:,1,1))),jointExtremes(:,1,2),'ks')

plot(datetime(datevec(jointExtremes(:,2,1))),jointExtremes(:,2,2),'r+')


legend('Timeseries1','TimeSeries2','Pks1','Pks2')



