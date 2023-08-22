close all;
clc
clearvars;
addpath('../');

load testcase_multivariate %water level data of two stations in the Adriatic Sea plus another one from the Dutch coasts and another from Carcavelos
if any(timeAndSeries1(:,1)~=timeAndSeries2(:,1))
    disp('time stamps of each series must be indentical- interpolate beforehand')
    return
end
thresholdPercentiles=[99,99,99,99]; %threshold levels in each series
minPeakDistanceInDays=[3,3,3,3]; %minimum distance between monovariate peaks of each series
size(nchoosek([1:4],2),1); %this is the size of maxDistanceMultivariatePeakInDays
maxDistanceMultivariatePeaksInDays=[7,7,7,7,7,7]; %has to have a format and size matching size(nchoosek([1:numvar],2),1) where numvar is 2 in bivariate case, 3 in trivariate case, and so on

timeWindow = 365.25*15; %to be used in tsEvaNonStationary
ciPercentile = 98;  %to be used in tsEvaNonStationary, not to be confused with thresholdPercentiles


[~,isort1]=sort(timeAndSeries1(:,1),'ascend');
[~,isort2]=sort(timeAndSeries2(:,1),'ascend');
[~,isort3]=sort(timeAndSeries3(:,1),'ascend');
[~,isort4]=sort(timeAndSeries4(:,1),'ascend');

timeAndSeries1=timeAndSeries1(isort1,:);
timeAndSeries2=timeAndSeries2(isort2,:);
timeAndSeries3=timeAndSeries3(isort3,:);
timeAndSeries4=timeAndSeries4(isort4,:);


[nonStatEvaParams1, statTransfData1] = tsEvaNonStationary(timeAndSeries1, timeWindow, 'transfType', 'trendlinear',...
    'ciPercentile', ciPercentile,'potPercentiles',[97:0.5:99], 'minPeakDistanceInDays', 14);

plot(datetime(datevec(timeAndSeries1(:,1))),statTransfData1.stationarySeries,'--')
hold on

[nonStatEvaParams2, statTransfData2] = tsEvaNonStationary(timeAndSeries2, timeWindow, 'transfType', 'trendlinear',...
    'ciPercentile', ciPercentile,'potPercentiles',[97:0.5:99], 'minPeakDistanceInDays', 14);

plot(datetime(datevec(timeAndSeries1(:,1))),statTransfData2.stationarySeries,'--')

[nonStatEvaParams3, statTransfData3] = tsEvaNonStationary(timeAndSeries3, timeWindow, 'transfType', 'trendlinear',...
    'ciPercentile', ciPercentile,'potPercentiles',[97:0.5:99], 'minPeakDistanceInDays', 14);

plot(datetime(datevec(timeAndSeries1(:,1))),statTransfData3.stationarySeries,'--')

[nonStatEvaParams4, statTransfData4] = tsEvaNonStationary(timeAndSeries4, timeWindow, 'transfType', 'trendlinear',...
    'ciPercentile', ciPercentile,'potPercentiles',[97:0.5:99], 'minPeakDistanceInDays', 14);

plot(datetime(datevec(timeAndSeries1(:,1))),statTransfData4.stationarySeries,'--')

[jointExtremes,jointExtremes2,thresholds,tstotal,pkstotal] = tsCopulaSampleJointPeaksMultiVariatePruning(timeAndSeries1(:,1), ...
    [statTransfData1.stationarySeries,statTransfData2.stationarySeries,statTransfData3.stationarySeries,statTransfData4.stationarySeries], ...
    thresholdPercentiles, ...
    minPeakDistanceInDays, ...
    maxDistanceMultivariatePeaksInDays);

plot(datetime(datevec(jointExtremes(:,1,1))),jointExtremes(:,1,2),'ks')   %plot of joint peak events belonging to first series

plot(datetime(datevec(jointExtremes(:,2,1))),jointExtremes(:,2,2),'r+')   %plot of joint peak events belonging to second series
plot(datetime(datevec(jointExtremes(:,3,1))),jointExtremes(:,3,2),'b>')   %plot of joint peak events belonging to third series
plot(datetime(datevec(jointExtremes(:,4,1))),jointExtremes(:,4,2),'y<')   %plot of joint peak events belonging to fourth series

plot(datetime(datevec(jointExtremes2(:,1,1))),jointExtremes2(:,1,2),'gp')   %plot of joint non-peak events belonging to first series (a non-peak event is defined as a concurrent event where at least one series exceeds its corresponding threshold and not all series exceed their thresholds

plot(datetime(datevec(jointExtremes2(:,2,1))),jointExtremes2(:,2,2),'mh')   %plot of joint non-peak events belonging to second series
plot(datetime(datevec(jointExtremes2(:,3,1))),jointExtremes2(:,3,2),'cd')   %plot of joint non-peak events belonging to third series
plot(datetime(datevec(jointExtremes2(:,4,1))),jointExtremes2(:,4,2),'k|')   %plot of joint non-peak events belonging to fourth series

line(datetime(datevec(timeAndSeries1(:,1))),thresholds(1)*ones(1,length(timeAndSeries1(:,1)))) %line plot of threshold of first series
line(datetime(datevec(timeAndSeries2(:,1))),thresholds(2)*ones(1,length(timeAndSeries1(:,1))))  %line plot of threshold of second series
line(datetime(datevec(timeAndSeries3(:,1))),thresholds(3)*ones(1,length(timeAndSeries1(:,1))))  %line plot of threshold of third series
line(datetime(datevec(timeAndSeries4(:,1))),thresholds(4)*ones(1,length(timeAndSeries1(:,1))))  %line plot of threshold of fourth series

plot(datetime(datevec(tstotal)),pkstotal,'.k')
hl=legend('Series1','Series2','Series3','Series4','JointPk1','JointPk2','JointPk3','JointPk4','JointNpk1','JointNpk2','JointNpk3','JointNpk4','Thr1','Thr2','Thr3','Thr4','Peaks');


for ij=1:size(jointExtremes2,1)
    x=datetime(datevec((min(jointExtremes2(ij,:,1),[],2))));
    y=min(jointExtremes2(ij,:,2),[],2);
    w=datetime(datevec((min(jointExtremes2(ij,:,1),[],2))+max(jointExtremes2(ij,:,1),[],2)-min(jointExtremes2(ij,:,1),[],2)));
    x=[x,w];
    h=max(jointExtremes2(ij,:,2),[],2)-min(jointExtremes2(ij,:,2),[],2);
    y=[y,y+h];
    h=patch(x([1,2,2,1]),y([1,1,2,2]),'r');
    h.EdgeColor='k';
    h.FaceAlpha=0;
end

for ij=1:size(jointExtremes,1)
    x=datetime(datevec((min(jointExtremes(ij,:,1),[],2))));
    y=min(jointExtremes(ij,:,2),[],2);
    w=datetime(datevec((min(jointExtremes(ij,:,1),[],2))+max(jointExtremes(ij,:,1),[],2)-min(jointExtremes(ij,:,1),[],2)));
    x=[x,w];
    h=max(jointExtremes(ij,:,2),[],2)-min(jointExtremes(ij,:,2),[],2);
    y=[y,y+h];
    h=patch(x([1,2,2,1]),y([1,1,2,2]),'r');
    h.EdgeColor='r';
    h.FaceAlpha=0;
end

Str_leg=hl.String;
hl.String=Str_leg(1:17);
set(gcf,'position',1.0e+03 *[ 0.0010    0.0410    1.5360    0.7488])
set(gca,'position',[0.0724    0.1100    0.8958    0.8654])
ylabel('TWL (m)')
xlabel('Date (time)')
