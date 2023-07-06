function sampleSPISeries_Gumbel

%%
% In a SPI (Standardized Precipitation Index) series peaks are distant at
% least 5 month one.
% It is possible to set the algorithm in order to fit a Gumbel instead of a
% full GEV
%%

addpath('../');

load('timeAndSeries_SPI_179_750E_-16.750N.mat');
timeAndSeries(:,2) = -timeAndSeries(:,2);

timeWindow = 50*315.25;
minPeakDistanceInDays = 5*30.2;
potEventsPerYear = 1;
returnPeriodsInYears = [10 20 50 100];

[nonStationaryEvaParams, stationaryTransformData, isValid] = tsEvaNonStationary( timeAndSeries, timeWindow, 'minPeakDistanceInDays', minPeakDistanceInDays, ...
  'transfType', 'trendCIPercentile', 'cipercentile', 80, 'gevType', 'gumbel', 'evdType', 'GEV');
tsEvaPlotSeriesTrendStdDevFromAnalyisObj(nonStationaryEvaParams, stationaryTransformData, 'plotpercentile', 95., 'ylabel', '-SPI', 'legendLocation', 'southwest');

[returnLevels, returnLevelsErr] = tsEvaComputeReturnLevelsGEVFromAnalysisObj(nonStationaryEvaParams, returnPeriodsInYears);
returnLevels = returnLevels*-1



