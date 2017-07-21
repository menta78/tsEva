function sampleTASMaxSeries

%%
% A series of TAS (Temperature of Air Surface) yerly maxima is useful to understand how heat waves evolve (Alessandro Dosio, jrc).
% It is a series of yearly maxima, therefore it is fit for a GEV analysis, while a GPD analysis is meaningless.
%%

addpath('../');

dt = load('timeAndSeriesTASMax.mat');
timeAndSeries = dt.timeAndSeries;
clear dt;

timeWindow = 50*315.25;
minPeakDistanceInDays = 5*30.2;
returnPeriodsInYears = [20 50 100 300];

[nonStationaryEvaParams, stationaryTransformData, isValid] = tsEvaNonStationary( timeAndSeries, timeWindow, 'minPeakDistanceInDays', minPeakDistanceInDays, 'extremeLowThreshold', .1, 'evdType', 'GEV');

tsEvaPlotSeriesTrendStdDevFromAnalyisObj(nonStationaryEvaParams, stationaryTransformData, 'ylabel', 'TAS', 'legendLocation', 'northwest');
ylim([0 40])
text(datenum(2060, 1, 1), 35, 'Series and trends', 'fontsize', 25);

tsEvaPlotGEVImageScFromAnalysisObj((0:.001:40)', nonStationaryEvaParams, stationaryTransformData, 'ylabel', 'TAS');
text(datenum(1980, 1, 1), 35, 'Time varying GEV', 'fontsize', 30);

[returnLevels, returnLevelsErr] = tsEvaComputeReturnLevelsGEVFromAnalysisObj(nonStationaryEvaParams, returnPeriodsInYears);

timeIndex = 26;
rlRange = [0 14]
hndl = tsEvaPlotReturnLevelsGEVFromAnalysisObj(nonStationaryEvaParams, timeIndex, 'ylim', rlRange, 'ylabel', 'return levels (TAS)');
ax = gca;
ax.YTick = 0:2:16;
text(7, 12.5, 'Return level 1995', 'fontsize', 30);


timeIndex = size(timeAndSeries, 1) - 4;
rlRange = [0 70]
hndl = tsEvaPlotReturnLevelsGEVFromAnalysisObj(nonStationaryEvaParams, timeIndex, 'ylim', rlRange, 'ylabel', 'return levels (TAS)');
text(7, 65, 'Return level 2095', 'fontsize', 30);

