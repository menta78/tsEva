# tsEVA 2.0 - Monovariate Examples Reference

This document contains all example scripts for monovariate extreme value analysis using the tsEVA 2.0 MATLAB toolbox.

---

## Example 1: Stationary EVA

**File:** `exampleEVAStationary.m`

**Purpose:** Demonstrates stationary extreme value analysis (GEV and GPD) on a time series.

**Key Features:**
- Stationary fit of GEV and GPD distributions
- Return level computation and plotting
- Fixed threshold POT analysis
- Gumbel distribution as alternative to full GEV

**Main Functions Used:**
- `tsEvaStationary()` - performs stationary EVA
- `tsEvaComputeReturnLevelsGEVFromAnalysisObj()` - computes GEV return levels
- `tsEvaPlotReturnLevelsGEVFromAnalysisObj()` - plots GEV return levels
- `tsEvaComputeReturnLevelsGPDFromAnalysisObj()` - computes GPD return levels
- `tsEvaPlotReturnLevelsGPDFromAnalysisObj()` - plots GPD return levels

**Key Parameters:**
- `minPeakDistanceInDays` - minimum distance between peaks
- `potThreshold` - fixed threshold for POT (optional)
- `gevtype` - can specify 'gumbel' instead of full GEV

**Code:**
```matlab
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
```

---

## Example 2: Non-Stationary EVA with Trend

**File:** `exampleGenerateSeriesEVAGraphs.m`

**Purpose:** Analyzes non-stationary time series with trend and seasonal components.

**Key Features:**
- Transformed-stationary approach for non-stationary data
- Trend-only transformation
- Seasonal transformation
- 2D and 3D visualization of time-varying distributions
- Return level computation at specific time indices

**Main Functions Used:**
- `tsEvaNonStationary()` - performs non-stationary EVA with transformation
- `tsEvaPlotSeriesTrendStdDevFromAnalyisObj()` - plots series with trend and std dev
- `tsEvaPlotGEVImageScFromAnalysisObj()` - 2D plot of time-varying GEV
- `tsEvaPlotGPDImageScFromAnalysisObj()` - 2D plot of time-varying GPD
- `tsEvaPlotGEV3DFromAnalysisObj()` - 3D plot of GEV
- `tsEvaPlotReturnLevelsGEVFromAnalysisObj()` - return level plot at specific time
- `tsEvaPlotTransfToStatFromAnalysisObj()` - plots transformed stationary series

**Key Parameters:**
- `timeWindow` - time window for detecting non-stationarity (e.g., 365.25*6 for 6 years)
- `transfType` - transformation type: 'trend' or 'seasonal'
- `minPeakDistanceInDays` - minimum distance between peaks

**Code:**
```matlab
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
saveas(hndl{1}, 'GEV_ReturnLevels.png', 'png');
hndl = tsEvaPlotReturnLevelsGPDFromAnalysisObj(nonStatEvaParams, timeIndex, 'ylim', [.5 1.5]);
saveas(hndl{1}, 'GPD_ReturnLevels.png', 'png');


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
```

---

## Example 3: Confidence Interval via Percentile

**File:** `exampleGenerateSeriesEVAGraphs_ciPercentile.m`

**Purpose:** Estimates long-term extreme variations using moving percentile instead of moving standard deviation.

**Key Features:**
- Uses moving percentile to estimate amplitude (more sensitive to extreme changes)
- Trend with CI percentile: `transfType = 'trendCiPercentile'`
- Seasonal with CI percentile: `transfType = 'seasonalCiPercentile'`
- Broader confidence intervals but better extreme modeling

**Additional Parameters:**
- `ciPercentile` - percentile for confidence interval estimation (e.g., 98)

**Code:**
```matlab
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
```

---

## Example 4: Negative Shape Parameter (GPD)

**File:** `exampleGenerateSeriesEVAGraphs_gpdNegShapeParam.m`

**Purpose:** Handles cases with negative GPD shape parameters (bounded upper tail).

**Key Features:**
- Linear trend transformation: `transfType = 'trendlinear'`
- Multiple POT percentiles: `potPercentiles = [97:0.5:99]`
- Handles negative shape parameter epsilon
- Wave height (SWH) analysis example

**Additional Parameters:**
- `potPercentiles` - array of percentiles to test for POT threshold

**Code:**
```matlab
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
```

---

## Example 5: Linear Trend Analysis

**File:** `exampleGenerateSeriesEVAGraphs_trendLinear.m`

**Purpose:** Analyzes total water level (TWL) with linear trend transformation.

**Key Features:**
- Linear trend for coastal flood applications
- Multiple POT percentiles
- Return level comparison at beginning and end of series
- Plotting POT and return levels together

**Main Functions Used:**
- `tsPlotSeriesPotGPDRetLevFromAnalysisObj` - plots the time series, the POT, and some return levels. Handy to diagnose
- `tsPlotSeriesYearMaxGEVRetLevFromAnalysisObj` - plots the time series, the annual maxima, and some return levels. Handy to diagnose

**Code:**
```matlab
load EOatSEE.mat
seriesDescr = "Adriatic TWL";

% temp = ncread(netcdf_filename,'date');
% ref = datenum(1950,1,1,0,0,0);
% tm= ref + temp;
% twl= ncread(netcdf_filename,'TWL');

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

disp('  plotting the year maxima and some GEV return levels');
hndl = tsPlotSeriesYearMaxGEVRetLevFromAnalysisObj( nonStatEvaParams, statTransfData);
saveas(hndl{1}, 'YearMaxGEVReturnLevelsLinearTrend.png');    

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
```

---

## Example 6: SPI Series (Sparse Peaks)

**File:** `exampleSPISeries.m`

**Purpose:** Analyzes Standardized Precipitation Index where peaks are widely separated.

**Key Features:**
- Large minimum peak distance (5 months)
- GPD-only analysis (no GEV due to sparse annual maxima)
- Single POT percentile
- Inverted series (analyzing droughts as negative SPI)

**Key Parameters:**
- `evdType = 'GPD'` - performs only GPD analysis, skips GEV
- `potPercentiles` - single value (e.g., 80) instead of array

**Code Structure:**
```matlab
timeWindow = 50*315.25;
minPeakDistanceInDays = 5*30.2;  % 5 months

[nonStationaryEvaParams, stationaryTransformData, isValid] = tsEvaNonStationary(
    timeAndSeries, timeWindow, 
    'minPeakDistanceInDays', minPeakDistanceInDays,
    'transfType', 'trendCIPercentile', 
    'cipercentile', 80, 
    'potPercentiles', 80, 
    'evdType', 'GPD');

[returnLevels, returnLevelsErr] = tsEvaComputeReturnLevelsGPDFromAnalysisObj(
    nonStationaryEvaParams, returnPeriodsInYears);
```

---

## Example 7: SPI with Gumbel

**File:** `exampleSPISeries_Gumbel.m`

**Purpose:** Fits Gumbel distribution (GEV with shape parameter = 0) to SPI data.

**Key Features:**
- Gumbel as special case of GEV
- GEV-only analysis
- Suitable when tail behavior suggests exponential decay

**Key Parameters:**
- `gevType = 'gumbel'` - fits Gumbel instead of full GEV
- `evdType = 'GEV'` - performs only GEV analysis

**Code Structure:**
```matlab
[nonStationaryEvaParams, stationaryTransformData, isValid] = tsEvaNonStationary(
    timeAndSeries, timeWindow,
    'minPeakDistanceInDays', minPeakDistanceInDays,
    'transfType', 'trendCIPercentile',
    'cipercentile', 80,
    'gevType', 'gumbel',
    'evdType', 'GEV');

[returnLevels, returnLevelsErr] = tsEvaComputeReturnLevelsGEVFromAnalysisObj(
    nonStationaryEvaParams, returnPeriodsInYears);
```

---

## Example 8: Temperature Series (Annual Maxima)

**File:** `exampleTASMaxSeries.m`

**Purpose:** Analyzes annual maximum temperature (heat wave evolution).

**Key Features:**
- Yearly maxima series - ideal for GEV, GPD meaningless
- GEV-only analysis
- Low threshold for extremes: `extremeLowThreshold = 0.1`
- Multiple time point comparisons (e.g., 1995 vs 2095)

**Key Parameters:**
- `evdType = 'GEV'` - only GEV, no GPD
- `extremeLowThreshold` - threshold below which values are not considered

**Code Structure:**
```matlab
[nonStationaryEvaParams, stationaryTransformData, isValid] = tsEvaNonStationary(
    timeAndSeries, timeWindow,
    'minPeakDistanceInDays', minPeakDistanceInDays,
    'extremeLowThreshold', .1,
    'evdType', 'GEV');

% Compare return levels at different times
timeIndex = 26;  % 1995
hndl = tsEvaPlotReturnLevelsGEVFromAnalysisObj(nonStationaryEvaParams, timeIndex);

timeIndex = size(timeAndSeries, 1) - 4;  % 2095
hndl = tsEvaPlotReturnLevelsGEVFromAnalysisObj(nonStationaryEvaParams, timeIndex);
```

---

## Common Workflow Patterns

### Basic Non-Stationary Analysis
```matlab
% 1. Load data
timeAndSeries = [timestamps, values];

% 2. Set parameters
timeWindow = 365.25 * years;
minPeakDistanceInDays = days;

% 3. Perform analysis
[nonStatEvaParams, statTransfData] = tsEvaNonStationary(
    timeAndSeries, timeWindow, 
    'transfType', 'trend',  % or 'seasonal', 'trendlinear', 'trendCiPercentile'
    'minPeakDistanceInDays', minPeakDistanceInDays);

% 4. Visualize
tsEvaPlotSeriesTrendStdDevFromAnalyisObj(nonStatEvaParams, statTransfData);
tsEvaPlotGEVImageScFromAnalysisObj(wr, nonStatEvaParams, statTransfData);
tsEvaPlotGPDImageScFromAnalysisObj(wr, nonStatEvaParams, statTransfData);

% 5. Compute return levels
[rlevGPD, rlevGPDErr] = tsEvaComputeReturnLevelsGPDFromAnalysisObj(
    nonStatEvaParams, [10, 20, 50, 100], 'timeindex', timeIndex);
```

### Transformation Types
- `'trend'` - trend + moving std dev
- `'seasonal'` - seasonal + moving std dev
- `'trendCiPercentile'` - trend + moving percentile
- `'seasonalCiPercentile'` - seasonal + moving percentile
- `'trendlinear'` - linear trend

### Distribution Types
- `'evdType', 'GPD'` - GPD only
- `'evdType', 'GEV'` - GEV only
- Default: both GEV and GPD

### Special Options
- `'gevType', 'gumbel'` - fit Gumbel instead of full GEV
- `'potPercentiles', [97:0.5:99]` - multiple thresholds to test
- `'potThreshold', value` - fixed POT threshold
- `'ciPercentile', 98` - percentile for CI estimation
- `'extremeLowThreshold', 0.1` - minimum value threshold

---



## References

Mentaschi, L., et al. (2016). The transformed-stationary approach: a generic and simplified methodology for non-stationary extreme value analysis. Hydrol. Earth Syst. Sci., 20, 3527-3547.

Bahmanpour, F., et al. (2025). Transformed-Stationary EVA 2.0: A Generalized Framework for Non-Stationary Joint Extremes Analysis. Hydrol. Earth Syst. Sci. (under review).
