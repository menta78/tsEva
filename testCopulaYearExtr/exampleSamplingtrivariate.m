close all;
clc
clearvars;

addpath('../');

path(path,'/Users/bahmanpour/Downloads/tsEva-dvlp-copula_git5/testCopulaYearExtr/data')   %netcdf files need to be here
filenames={'jd_twl_ADRIATIC_NORTH.nc','jd_twl_ADRIATIC_SOUTH.nc','jd_twl_DUTCH_COAST.nc','jd_twl_CARCAVELOS.nc'};  %use this for case of trivariate example, only the first three filenames are used
 filenames={'jd_twl_ADRIATIC_NORTH.nc','jd_twl_ADRIATIC_SOUTH.nc'};   %use this for case of bivariate example


for ii=1:length(filenames)

    var0=datenum(1950,1,1)+double(ncread(filenames{ii},'date'));
    var1=double(ncread(filenames{ii},'TWL'));

    if ii==1
        timeAndSeries1=[var0,var1];
    elseif ii==2
        timeAndSeries2=[var0,var1];
    elseif ii==3
        timeAndSeries3=[var0,var1];

    else

        timeAndSeries4=[var0,var1]; %this will not be used
    end
end
if any(timeAndSeries1(:,1)~=timeAndSeries2(:,1))
    disp('time stamps of each series are not indentical- an interpolation needs to be done')
    return
end


thresholdPercentiles=[99,99,99]; %for the calculation of threshold value in each series
minPeakDistanceInDays=[3,3,3]; %minimum distance between monovariate peaks to be selected by findpeaks function
maxDistanceMultivariatePeaksInDays=[7,7,7]; %maximum distance between multivariates, first entry is the distance between first and second series,
                                             % second entry is the distance between second and third series and third entry is the distance between first and third series

timeWindow = 365.25*15; %to be used in tsEvaNonStationary
ciPercentile = 98;  %to be used in tsEvaNonStationary


[~,isort1]=sort(timeAndSeries1(:,1),'ascend');
[~,isort2]=sort(timeAndSeries2(:,1),'ascend');

try
    [~,isort3]=sort(timeAndSeries3(:,1),'ascend');
catch

end

timeAndSeries1=timeAndSeries1(isort1,:);
timeAndSeries2=timeAndSeries2(isort2,:);

try
    timeAndSeries3=timeAndSeries3(isort3,:);

catch

end

[nonStatEvaParams1, statTransfData1] = tsEvaNonStationary(timeAndSeries1, timeWindow, 'transfType', 'trendlinear',...
    'ciPercentile', ciPercentile,'potPercentiles',[97:0.5:99], 'minPeakDistanceInDays', 14);
plot(datetime(datevec(timeAndSeries1(:,1))),statTransfData1.stationarySeries,'--')

[nonStatEvaParams2, statTransfData2] = tsEvaNonStationary(timeAndSeries2, timeWindow, 'transfType', 'trendlinear',...
    'ciPercentile', ciPercentile,'potPercentiles',[97:0.5:99], 'minPeakDistanceInDays', 14);
hold on
plot(datetime(datevec(timeAndSeries1(:,1))),statTransfData2.stationarySeries,'--')

try
    [nonStatEvaParams3, statTransfData3] = tsEvaNonStationary(timeAndSeries3, timeWindow, 'transfType', 'trendlinear',...
        'ciPercentile', ciPercentile,'potPercentiles',[97:0.5:99], 'minPeakDistanceInDays', 14);
    plot(datetime(datevec(timeAndSeries3(:,1))),statTransfData3.stationarySeries,'--')
catch

end

try  %trivariate case
[jointExtremes,jointExtremes2, thresholds,tstotal,pkstotal] = tsCopulaSampleJointPeaks(timeAndSeries1(:,1), ...
    [statTransfData1.stationarySeries,statTransfData2.stationarySeries,statTransfData3.stationarySeries], ...
    thresholdPercentiles, ...
    minPeakDistanceInDays, ...
    maxDistanceMultivariatePeaksInDays);
catch   %bivariate case
[jointExtremes,jointExtremes2, thresholds,tstotal,pkstotal] = tsCopulaSampleJointPeaks(timeAndSeries1(:,1), ...
    [statTransfData1.stationarySeries,statTransfData2.stationarySeries], ...
    thresholdPercentiles, ...
    minPeakDistanceInDays, ...
    maxDistanceMultivariatePeaksInDays);
end

% plot(datetime(datevec(tstotal)),pkstotal,'k.')  % plot of all peaks returned by the findpeaks function

plot(datetime(datevec(jointExtremes(:,1,1))),jointExtremes(:,1,2),'ks')   %plot of peaks belonging to first series (joint peak events)

plot(datetime(datevec(jointExtremes(:,2,1))),jointExtremes(:,2,2),'r+')   %plot of peaks belonging to second series   (joint peak events)

try
    plot(datetime(datevec(jointExtremes(:,3,1))),jointExtremes(:,3,2),'b>')   %plot of peaks belonging to third series  (joint peak events)
catch

end


plot(datetime(datevec(jointExtremes2(:,1,1))),jointExtremes2(:,1,2),'gp')   %plot of peaks belonging to first series (joint non-peak events meaning that at least one peak (in one series) exceeded its respective threshold)

plot(datetime(datevec(jointExtremes2(:,2,1))),jointExtremes2(:,2,2),'mh')   %%plot of peaks belonging to second series (joint non-peak events meaning that at least one peak (in one series) exceeded its respective threshold)

try
    plot(datetime(datevec(jointExtremes2(:,3,1))),jointExtremes2(:,3,2),'cd') %%plot of peaks belonging to third series (joint non-peak events meaning that at least one peak (in one series) exceeded its respective threshold)
catch

end


 line(datetime(datevec(timeAndSeries1(:,1))),thresholds(1)*ones(1,length(timeAndSeries1(:,1)))) %plot of respective threshold of first series
 line(datetime(datevec(timeAndSeries2(:,1))),thresholds(2)*ones(1,length(timeAndSeries1(:,1))))  %plot of respective threshold of second series

 try
     line(datetime(datevec(timeAndSeries3(:,1))),thresholds(3)*ones(1,length(timeAndSeries1(:,1))))  %plot of respective threshold of third series

     legend('Timeseries1','TimeSeries2','TimeSeries3','jointpks1','jointpks2','jointpks3','nPks1','nPks2','nPks3','th1','th2','th3')

 catch
     legend('Timeseries1','TimeSeries2','jointpks1','jointpks2','nPks1','nPks2','th1','th2')

 end
