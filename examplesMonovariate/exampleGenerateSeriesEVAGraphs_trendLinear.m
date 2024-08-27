%%
% this sample script illustrates how to execute the tsEva to estimate the
% long term variations of the extremes, using a moving percentile to
% estimate the amplitude of the series, instead of the moving standard
% deviation. This approach models better the variations of the extremes
% than the one based on the standard deviation, but is subject to stronger
% uncertainty.
%%

clearvars
addpath('../');

netcdf_filename = "jd_twl_ADRIATIC_NORTH.nc";
seriesDescr = "Adriatic TWL";

temp = ncread(netcdf_filename,'date');
ref = datenum(1950,1,1,0,0,0);
tm= ref + temp;
twl= ncread(netcdf_filename,'TWL');

timeAndSeries =[tm,twl];
extremesRange = [.5 2];
seasonalExtrRange = [.1 1.1];

timeWindow = 365.25*15; % 6 years
minPeakDistanceInDays = 14;
% minPeakDistanceInDays = 3;

minTS = min(timeAndSeries(:,1));
maxTS = max(timeAndSeries(:,1));
axisFontSize = 20;
axisFontSize3d = 16;
labelFontSize = 28;
titleFontSize = 30;

% preparing xticks
% years = (1980:2:2015)';
years = (1994:2:2022)';
months = ones(size(years));
days = ones(size(years));
dtns = cat(2, years, months, days);
tickTmStmp = datenum(dtns);

wr = linspace(min(extremesRange), max(extremesRange), 1501);

ciPercentile = 99;

disp('trend only statistics (transformation + eva + backtransformation)');
[nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, timeWindow, 'transfType', 'trendlinear',... 
  'ciPercentile', ciPercentile, 'potPercentiles',[97:0.5:99], 'minPeakDistanceInDays', minPeakDistanceInDays);
disp('  plotting the series');
aax=max(statTransfData.stdDevSeries);
bbx=max(statTransfData.nonStatSeries);
ul=bbx+aax;
ll=bbx-2*aax;
rlRange = [ll ul];
hndl = tsEvaPlotSeriesTrendStdDevFromAnalyisObj(nonStatEvaParams, statTransfData, 'legendLocation', 'northeast', ...
    'ylabel', 'TWL(m)', 'title', seriesDescr, 'titleFontSize', titleFontSize, 'dateformat', 'yy', 'xtick', tickTmStmp,'Interpreter','none');
disp('  saving the series plot');
saveas(hndl{1}, 'seriesTrendLinear.png');    

disp('  plotting the POT and some return levels');
hndl = tsPlotSeriesPotGPDRetLevFromAnalysisObj( nonStatEvaParams, statTransfData);
saveas(hndl{1}, 'PotAndReturnLevelsLinearTrend.png');    

disp('  plotting and saving the 2D GPD graph');
hndl = tsEvaPlotGPDImageScFromAnalysisObj(wr, nonStatEvaParams, statTransfData, 'ylabel', 'TWL(m)', 'dateformat', 'yy', 'xtick', tickTmStmp);
title('GPD', 'fontsize', titleFontSize);
% saveas(hndl{1}, [str{jx},'GPD2DTrendLinear.png'], 'png');

%computing and plotting the return levels for a given times
for lx=1:2
    if lx==1
        timeIndex = 1000;
        ttl='GPD return levels for beginning of the series';
    else
        timeIndex = 240000;
        ttl='GPD return levels for end of the series';
    end
    timeStamps = statTransfData.timeStamps;
    dtvc = datevec(timeStamps(timeIndex));
    tmstmpref = datenum(dtvc(1), dtvc(2), 1);
    [rlevGPD, rlevGPDErr] = tsEvaComputeReturnLevelsGPDFromAnalysisObj(nonStatEvaParams, [5, 10, 30, 100], 'timeindex', timeIndex);
    hndl = tsEvaPlotReturnLevelsGPDFromAnalysisObj(nonStatEvaParams, timeIndex, 'ylim', rlRange,'maxReturnPeriodYears', 200);
    % title(['GPD return levels for ' datestr(tmstmpref)]);
    title(ttl);
    if lx==1
        saveas(hndl{1}, 'GPD_ReturnLevels-end.png', 'png');
    else
        saveas(hndl{1}, 'GPD_ReturnLevels-beg.png', 'png');
    end
end
    
disp('plotting and saving stationary series');
hndl = tsEvaPlotTransfToStatFromAnalysisObj(nonStatEvaParams, statTransfData, 'dateformat', 'yy', 'xtick', tickTmStmp);    
saveas(hndl{1}, 'statSeries_trendLinear.png', 'png');
