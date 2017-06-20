%%
% test/sampleEVAStationary.m: Sample application of the stationary
% approach.
%%

addpath('../');

load('timeAndSeriesHebrides.mat');
timeAndSeries = timeAndSeriesHebrides;

minPeakDistanceInDays = 3;

disp('stationary fit of extreme value distributions (GEV, GPD) to a time series');
statEvaParams = tsEvaStationary(timeAndSeries, 'minPeakDistanceInDays', minPeakDistanceInDays);

%computing and plotting the return levels for a given times
hndl = tsEvaPlotReturnLevelsGEVFromAnalysisObj(statEvaParams, 1, 'ylim', [.5 1.5]);
saveas(hndl{1}, 'GEV_ReturnLevels_STATIONARY.png', 'png');
hndl = tsEvaPlotReturnLevelsGPDFromAnalysisObj(statEvaParams, 1, 'ylim', [.5 1.5]);
saveas(hndl{1}, 'GPD_ReturnLevels_STATIONARY.png', 'png');
