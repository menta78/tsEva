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

[resampleLevel, resampleProb, resampleRetPer] = tsCopulaCompoundGPDMontecarlo(copulaAnalysis,...
    nResample,'timeIndex',1);


peakIndices=copulaAnalysis.jointExtremeIndices;
yMax=copulaAnalysis.jointExtremes; %non-stationary peaks
tMax=copulaAnalysis.jointExtremeTimeStamps;
thresholdC=copulaAnalysis.thresholdPotNS; %non-stationary threshold parameters
figure 

% the non-stationary series and threshold parameter of each 
% series is plotted to see if the peaks that were sampled from stationarized
% series also coincide with peaks in the non-stationary series

plot(datetime(datevec(timeAndSeries1(:,1))),timeAndSeries1(:,2))
hold on
plot(datetime(datevec(timeAndSeries2(:,1))),timeAndSeries2(:,2))
p1=copulaAnalysis.marginalAnalysis{1}{1}(2).objs.peakIndexes;
p2=copulaAnalysis.marginalAnalysis{2}{1}(2).objs.peakIndexes;

plot(datetime(datevec(timeAndSeries1(:,1))),thresholdC(:,1))
plot(datetime(datevec(timeAndSeries2(:,1))),thresholdC(:,2))
plot(datetime(datevec(tMax(:,1))),yMax(:,1),'.r')
plot(datetime(datevec(tMax(:,2))),yMax(:,2),'.k')
plot(datetime(datevec(timeAndSeries1(p1,1))),timeAndSeries1(p1,2),'sg')
plot(datetime(datevec(timeAndSeries2(p2,1))),timeAndSeries2(p2,2),'sb')
legend('Series1','Series2','Threshold-Param1','Threshold-param2','Peaks1','Peaks2','monoVarPeak-1','monoVarPeak-2')
title('Non-stationary series')




% plotting of jointExtremes and simulated return levels from the copula 
figHnd = tsCopulaPeakExtrPlotSctrBivar(resampleLevel, yMax, 'xlbl', 'Marshall-north', 'ylbl', 'Marshall-south');

%find the outliers
y1=yMax(:,1);
y2=yMax(:,2);
[xData, yData] = prepareCurveData( y1, y2 );

ft = fittype( 'poly1' );

[fitresult, gof] = fit( xData, yData, ft );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'y2 vs. y1', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'y1', 'Interpreter', 'none' );
ylabel( 'y2', 'Interpreter', 'none' );
grid on

xn=linspace(min(y1),max(y1),100);
yn=xn*(fitresult.p1)+fitresult.p2;
dx0=[];
for ij1=1:length(y1)
    dx=[];    
for ij=1:100
    dx=[dx,sqrt((y1(ij1)-xn(ij)).^2+(y2(ij1)-yn(ij)).^2)];
end
dx0=[dx0,min(dx)];
end
[~,idx0]=sort(dx0,'descend');
yMax(idx0(1:2),:)
datestr(tMax(idx0(1:2),:))
        
