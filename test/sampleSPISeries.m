function sampleSPISeries

%%
% In a SPI (Standardized Precipitation Index) series peaks are distant at least 5 month one
% from the other. Therefore the concept of "annual maxima" is meaningless, 
% as it is the concept of "5 peaks over threshold per year".
% It is possible to set the algorithm to examine just one percentile for
% the POT, and doing just the GPD analysis and not the GEV one. This example shows how.
%%

addpath('../');

load('timeAndSeries_SPI_179_750E_-16.750N.mat');
timeAndSeries(:,2) = -timeAndSeries(:,2);

timeWindow = 50*315.25;
minPeakDistanceInDays = 5*30.2;
potEventsPerYear = 1;
returnPeriodsInYears = [20 50 100 300];

[nonStationaryEvaParams, stationaryTransformData, isValid] = tsEvaNonStationary( timeAndSeries, timeWindow, 'minPeakDistanceInDays', minPeakDistanceInDays, ...
  'transfType', 'trendCIPercentile', 'cipercentile', 80, 'potPercentiles', 80, 'evdType', 'GPD');
tsEvaPlotSeriesTrendStdDevFromAnalyisObj(nonStationaryEvaParams, stationaryTransformData, 'plotpercentile', 95., 'ylabel', '-SPI', 'legendLocation', 'southwest');

[returnLevels, returnLevelsErr] = tsEvaComputeReturnLevelsGPDFromAnalysisObj(nonStationaryEvaParams, returnPeriodsInYears);
returnLevels = returnLevels*-1;



