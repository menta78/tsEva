addpath('../');

load('timeAndSeriesHebrides.mat');
timeAndSeries = timeAndSeriesHebrides;
extremesRange = [.2 1.2];
seasonalExtrRange = [.1 1.1];
seriesDescr = 'Hebrides';

timeWindow = 365.25*6; % 6 years
minPeakDistanceInDays = 3;

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
[nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, timeWindow, 'transfType', 'trend', 'minPeakDistanceInDays', minPeakDistanceInDays);
disp('  plotting the series');
hndl = tsEvaPlotSeriesTrendStdDevFromAnalyisObj(nonStatEvaParams, statTransfData,...
    'ylabel', 'Lvl (m)', 'title', seriesDescr, 'titleFontSize', titleFontSize, 'dateformat', 'yy', 'xtick', tickTmStmp);
disp('  saving the series plot');
saveas(hndl{1}, 'seriesTrendOnly.png');    
% disp('  plotting and saving the 3D GEV graph');
% hndl = tsEvaPlotGEV3DFromAnalysisObj(wr, nonStatEvaParams, statTransfData, 'xlabel', 'Lvl (m)', 'axisfontsize', axisFontSize3d);
% title('GEV 3D', 'fontsize', titleFontSize);
% saveas(hndl{1}, 'GEV3DTrendOnly.png', 'png');
disp('  plotting and saving the 2D GEV graph');
hndl = tsEvaPlotGEVImageScFromAnalysisObj(wr, nonStatEvaParams, statTransfData, 'ylabel', 'Lvl (m)', 'dateformat', 'yy', 'xtick', tickTmStmp);
title('GEV', 'fontsize', titleFontSize);
saveas(hndl{1}, 'GEV2DTrendOnly.png', 'png');
disp('  plotting and saving the 2D GPD graph');
hndl = tsEvaPlotGPDImageScFromAnalysisObj(wr, nonStatEvaParams, statTransfData, 'ylabel', 'Lvl (m)', 'dateformat', 'yy', 'xtick', tickTmStmp);
title('GPD', 'fontsize', titleFontSize);
saveas(hndl{1}, 'GPD2DTrendOnly.png', 'png');

%computing and plotting the return levels for a given times
timeIndex = 1000;
timeStamps = statTransfData.timeStamps;
disp(['  plotting return levels for time ' datestr(timeStamps(timeIndex))]);
disp('  ... for GEV the sample is small and the confidence interval is broad');
hndl = tsEvaPlotReturnLevelsGEVFromAnalysisObj(nonStatEvaParams, timeIndex, 'ylim', [.5 1.5]);
saveas(hndl{1}, 'GEVReturnLevels.png', 'png');
hndl = tsEvaPlotReturnLevelsGPDFromAnalysisObj(nonStatEvaParams, timeIndex, 'ylim', [.5 1.5]);
saveas(hndl{1}, 'GPDReturnLevels.png', 'png');


disp('plotting and saving stationary series');
hndl = tsEvaPlotTransfToStatFromAnalysisObj(nonStatEvaParams, statTransfData, 'dateformat', 'yy', 'xtick', tickTmStmp);
saveas(hndl{1}, 'statSeriesTrendOnly.png', 'png');

disp('seasonal statistics');
[nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, timeWindow, 'transfType', 'seasonal', 'minPeakDistanceInDays', minPeakDistanceInDays);

wr = linspace(min(seasonalExtrRange), max(seasonalExtrRange), 1501);

disp('  plotting a slice of data ');
slice = { 1988 1993};
plotTitle = '1988-1993';
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
