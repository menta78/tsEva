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
[rlevGEV, rlevGEVErr] = tsEvaComputeReturnLevelsGEVFromAnalysisObj(statEvaParams, [10, 20, 50, 100]);
rlevGEV
hndl = tsEvaPlotReturnLevelsGEVFromAnalysisObj(statEvaParams, 1, 'ylim', [.5 1.5]);
title('GEV');
saveas(hndl{1}, 'GEV_ReturnLevels_STATIONARY.png', 'png');

[rlevGPD, rlevGPDErr] = tsEvaComputeReturnLevelsGPDFromAnalysisObj(statEvaParams, [10, 20, 50, 100]);
rlevGPD
hndl = tsEvaPlotReturnLevelsGPDFromAnalysisObj(statEvaParams, 1, 'ylim', [.5 1.5]);
title('GPD');
saveas(hndl{1}, 'GPD_ReturnLevels_STATIONARY.png', 'png');



disp('Same as before, but the POT is done with a fixed threshold');
potThreshold = prctile(timeAndSeries(:,2), 98);
statEvaParams = tsEvaStationary(timeAndSeries, 'minPeakDistanceInDays', minPeakDistanceInDays, 'doSampleData', false, 'potThreshold', potThreshold);

%computing and plotting the return levels for a given times
[rlevGPD, rlevGPDErr] = tsEvaComputeReturnLevelsGPDFromAnalysisObj(statEvaParams, [10, 20, 50, 100]);
rlevGPD
hndl = tsEvaPlotReturnLevelsGPDFromAnalysisObj(statEvaParams, 1, 'ylim', [.5 1.5]);
title('GPD');
saveas(hndl{1}, 'GPD_ReturnLevels_STATIONARY.png', 'png');



disp('same thing as before, but with a gumbel instead of a GEV');
statEvaParams = tsEvaStationary(timeAndSeries, 'minPeakDistanceInDays', minPeakDistanceInDays, 'gevtype', 'gumbel');

%computing and plotting the return levels for a given times
[rlevGEV, rlevGEVErr] = tsEvaComputeReturnLevelsGEVFromAnalysisObj(statEvaParams, [10, 20, 50, 100]);
rlevGEV
hndl = tsEvaPlotReturnLevelsGEVFromAnalysisObj(statEvaParams, 1, 'ylim', [.5 1.5]);
title('Gumbel');
saveas(hndl{1}, 'GEV_ReturnLevels_STATIONARY.png', 'png');
