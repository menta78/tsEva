 %%
% this sample script illustrates how to execute the tsEva to estimate the
% long term variations of the extremes, using a moving percentile to
% estimate the amplitude of the series, instead of the moving standard
% deviation. This approach models better the variations of the extremes
% than the one based on the standard deviation, but is subject to stronger
% uncertainty.
%%

clearvars;
close all;
addpath('../');

dt = load("timeAndSeries_waves_SouthChina.mat");
timeAndSeries = dt.timeAndSeries;

extremesRange = [0 10]; 

timeWindow = 365.25*30; % 30 years
minPeakDistanceInDays = 3;
% minPeakDistanceInDays = 14;

axisFontSize = 20;
axisFontSize3d = 16;
labelFontSize = 28;
titleFontSize = 30;

% preparing xticks
% years = (1980:2:2015)';
years = (1950:5:2025)';
months = ones(size(years));
days = ones(size(years));
dtns = cat(2, years, months, days);
tickTmStmp = datenum(dtns);

wr = linspace(min(extremesRange), max(extremesRange), 1501);

ciPercentile = 99;

minTS = min(timeAndSeries(:,1));
maxTS = max(timeAndSeries(:,1));
seriesDescr = 'pt10';  % 得到 'pt0' 到 'pt9'

disp('trend only statistics (transformation + eva + backtransformation)');
[nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, timeWindow, 'transfType', 'trendlinear',... 
  'ciPercentile', ciPercentile, 'potPercentiles',[97:0.5:99], 'minPeakDistanceInDays', minPeakDistanceInDays);
disp('  plotting the series');


epsilon = nonStatEvaParams(2).parameters.epsilon;
pvalue = statTransfData.pValueChange;


aax=max(statTransfData.stdDevSeries);
bbx=max(statTransfData.nonStatSeries);
ul=bbx+aax;
ll=bbx-2*aax;
rlRange = [0 10];

tsEvaPlotSeriesTrendStdDevFromAnalyisObj(nonStatEvaParams, statTransfData, 'legendLocation', 'northeast', ...
    'ylabel', 'SWH(m)', 'title', seriesDescr, 'titleFontSize', titleFontSize, 'dateformat', 'yy', 'xtick', tickTmStmp,'Interpreter','none');
disp('  saving the series plot');
filename1 = 'seriesTrendLinear_pt10.png';
saveas(gcf, filename1);     

disp('  plotting the POT and some return levels');
tsPlotSeriesPotGPDRetLevFromAnalysisObj( nonStatEvaParams, statTransfData);
title('PotAndReturnLevelsLinearTrend pt10d', 'fontsize', 10);
filename2 = 'PotAndReturnLevelsLinearTrend_pt10.png';
saveas(gcf, filename2); 

disp('  plotting and saving the 2D GPD graph');
tsEvaPlotGPDImageScFromAnalysisObj(wr, nonStatEvaParams, statTransfData, 'ylabel', 'SWH(m)', 'dateformat', 'yy', 'xtick', tickTmStmp);
title('GPD pt10', 'fontsize', titleFontSize);
filename3 = 'GPD2DTrendLinear_pt10.png';
saveas(gcf, filename3);

%computing and plotting the return levels for a given times
timeIndex = 1; % 01/01/1950
timeStamps = statTransfData.timeStamps;
dtvc = datevec(timeStamps(timeIndex));
tmstmpref = datenum(dtvc(1), dtvc(2), 1);
disp(['  plotting return levels for time ' datestr(timeStamps(timeIndex))]);

[rlevGPD, rlevGPDErr] = tsEvaComputeReturnLevelsGPDFromAnalysisObj(nonStatEvaParams, [5, 10, 30, 100], 'timeindex', timeIndex);
tsEvaPlotReturnLevelsGPDFromAnalysisObj(nonStatEvaParams, timeIndex, 'ylim', rlRange,'maxReturnPeriodYears', 200);
yticks(0:1.5:10); 
    hold on;
legend(sprintf('negative shapeparemeter = %.3e', epsilon), 'Location', 'northwest');
title('pt10 GPD return levels-begin');
filename4 = 'GPD_ReturnLevels_ciPercentile_pt10.png';
saveas(gcf, filename4); 

timeIndex = 216199; %last date
timeStamps = statTransfData.timeStamps;
dtvc = datevec(timeStamps(timeIndex));
tmstmpref2 = datestr(timeStamps(timeIndex));
disp(['  plotting return levels for time ' datestr(timeStamps(timeIndex))]);

[rlevGPD, rlevGPDErr] = tsEvaComputeReturnLevelsGPDFromAnalysisObj(nonStatEvaParams, [5, 10, 30, 100], 'timeindex', timeIndex);
tsEvaPlotReturnLevelsGPDFromAnalysisObj(nonStatEvaParams, timeIndex, 'ylim', rlRange,'maxReturnPeriodYears', 200);
yticks(0:1.5:10); 
    hold on;
legend(sprintf('negative shapeparemeter = %.3e', epsilon), 'Location', 'northwest');
title('pt10 GPD return levels-end');
filename5 = 'GPD_ReturnLevels_ciPercentile_pt10_2023.png';
saveas(gcf,filename5); 
