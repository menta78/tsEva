%%
% this sample script illustrates how to execute the tsEva to estimate the
% long term variations of the extremes, using a moving percentile to
% estimate the amplitude of the series, instead of the moving standard
% deviation. This approach models better the variations of the extremes
% than the one based on the standard deviation, but is subject to stronger
% uncertainty.
%%

addpath('../');

load('timeAndSeriesHebrides.mat');
timeAndSeries = timeAndSeriesHebrides;
extremesRange = [.2 1.2];
rlRange = [.6 1.1];
seasonalExtrRange = [.1 1.1];
seriesDescr = 'Hebrides';

timeWindow = 365.25*6; % 6 years
minPeakDistanceInDays = 3;
ciPercentile = 98;

minTS = min(timeAndSeries(:,1));
maxTS = max(timeAndSeries(:,1));
axisFontSize = 20;
axisFontSize3d = 16;
labelFontSize = 28;
titleFontSize = 30;

% preparing xticks
years = (1980:2:2015)';
months = ones(size(years));
days = ones(size(years));
dtns = cat(2, years, months, days);
tickTmStmp = datenum(dtns);

wr = linspace(min(extremesRange), max(extremesRange), 1501);

disp('trend only statistics (transformation + eva + backtransformation)');
[nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, timeWindow, 'transfType', 'trendCiPercentile',... 
  'ciPercentile', ciPercentile, 'minPeakDistanceInDays', minPeakDistanceInDays);
disp('  plotting the series');
hndl = tsEvaPlotSeriesTrendStdDevFromAnalyisObj(nonStatEvaParams, statTransfData, 'legendLocation', 'northeast', ...
    'ylabel', 'Lvl (m)', 'title', seriesDescr, 'titleFontSize', titleFontSize, 'dateformat', 'yy', 'xtick', tickTmStmp);
disp('  saving the series plot');
saveas(hndl{1}, 'seriesTrendOnly_ciPercentile.png');    
% disp('  plotting and saving the 3D GEV graph');
% hndl = tsEvaPlotGEV3DFromAnalysisObj(wr, nonStatEvaParams, statTransfData, 'xlabel', 'Lvl (m)', 'axisfontsize', axisFontSize3d);
% title('GEV 3D', 'fontsize', titleFontSize);
% saveas(hndl{1}, 'GEV3DTrendOnly_ciPercentile.png', 'png');
disp('  plotting and saving the 2D GEV graph');
hndl = tsEvaPlotGEVImageScFromAnalysisObj(wr, nonStatEvaParams, statTransfData, 'ylabel', 'Lvl (m)', 'dateformat', 'yy', 'xtick', tickTmStmp);
title('GEV', 'fontsize', titleFontSize);
saveas(hndl{1}, 'GEV2DTrendOnly_ciPercentile.png', 'png');
disp('  plotting and saving the 2D GPD graph');
hndl = tsEvaPlotGPDImageScFromAnalysisObj(wr, nonStatEvaParams, statTransfData, 'ylabel', 'Lvl (m)', 'dateformat', 'yy', 'xtick', tickTmStmp);
title('GPD', 'fontsize', titleFontSize);
saveas(hndl{1}, 'GPD2DTrendOnly_ciPercentile.png', 'png');

%computing and plotting the return levels for a given times
timeIndex = 1000;
timeStamps = statTransfData.timeStamps;
dtvc = datevec(timeStamps(timeIndex));
tmstmpref = datenum(dtvc(1), dtvc(2), 1);
disp(['  plotting return levels for time ' datestr(timeStamps(timeIndex))]);
disp('  ... for GEV the sample is small and the confidence interval is broad');
[rlevGEV, rlevGEVErr] = tsEvaComputeReturnLevelsGEVFromAnalysisObj(nonStatEvaParams, [10, 20, 50, 100], 'timeindex', timeIndex);
rlevGEV
hndl = tsEvaPlotReturnLevelsGEVFromAnalysisObj(nonStatEvaParams, timeIndex, 'ylim', rlRange);
title(['GEV return levels for ' datestr(tmstmpref)]);
saveas(hndl{1}, 'GEV_ReturnLevels_ciPercentile.png', 'png');
[rlevGPD, rlevGPDErr] = tsEvaComputeReturnLevelsGPDFromAnalysisObj(nonStatEvaParams, [10, 20, 50, 100], 'timeindex', timeIndex);
rlevGPD
hndl = tsEvaPlotReturnLevelsGPDFromAnalysisObj(nonStatEvaParams, timeIndex, 'ylim', rlRange);
title(['GPD return levels for ' datestr(tmstmpref)]);
saveas(hndl{1}, 'GPD_ReturnLevels_ciPercentile.png', 'png');


disp('plotting and saving stationary series');
hndl = tsEvaPlotTransfToStatFromAnalysisObj(nonStatEvaParams, statTransfData, 'dateformat', 'yy', 'xtick', tickTmStmp);
saveas(hndl{1}, 'statSeriesTrendOnly_ciPercentile.png', 'png');


disp('seasonal statistics');
[nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, timeWindow, 'transfType', 'seasonalCiPercentile',... 
  'ciPercentile', ciPercentile, 'minPeakDistanceInDays', minPeakDistanceInDays);

wr = linspace(min(seasonalExtrRange), max(seasonalExtrRange), 1501);

disp('  plotting a slice of data ');
slice = { 1990 1995};
plotTitle = '1990-1995';
disp('    plotting the series');
hndl = tsEvaPlotSeriesTrendStdDevFromAnalyisObj(nonStatEvaParams, statTransfData,...
    'ylabel', 'Lvl (m)', 'title', plotTitle, 'minyear', slice{1}, 'maxyear', slice{2});
disp('    saving the series plot');
saveas(hndl{1}, 'seriesSeasonal.png');    
disp('plotting and saving stationary series');
hndl = tsEvaPlotTransfToStatFromAnalysisObj(nonStatEvaParams, statTransfData, 'dateformat', 'yy', 'minyear', slice{1}, 'maxyear', slice{2});
saveas(hndl{1}, 'statSeriesTrendOnly.png', 'png');
disp('    plotting and saving the 3D GEV graph');
hndl = tsEvaPlotGEV3DFromAnalysisObj(wr, nonStatEvaParams, statTransfData, 'xlabel', 'Lvl (m)', 'minyear', slice{1}, 'maxyear', slice{2}, 'axisfontsize', axisFontSize3d);
title(['GEV 3D, ' plotTitle], 'fontsize', titleFontSize);
saveas(hndl{1}, 'GEV3DSeasonal.png');
disp('    plotting and saving the 2D GEV graph');
hndl = tsEvaPlotGEVImageScFromAnalysisObj(wr, nonStatEvaParams, statTransfData, 'ylabel', 'Lvl (m)', 'minyear', slice{1}, 'maxyear', slice{2}, 'dateformat', 'yy');
title(['GEV ' plotTitle], 'fontsize', titleFontSize);
saveas(hndl{1}, 'GEV2DSeasonal.png');
disp('    plotting and saving the 2D GPD graph');
hndl = tsEvaPlotGPDImageScFromAnalysisObj(wr, nonStatEvaParams, statTransfData, 'ylabel', 'Lvl (m)', 'minyear', slice{1}, 'maxyear', slice{2}, 'dateformat', 'yy');
title(['GPD ' plotTitle], 'fontsize', titleFontSize);
saveas(hndl{1}, 'GPD2DSeasonal.png', 'png');

