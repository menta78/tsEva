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

**Code Structure:**
```matlab
% Load data
load('timeAndSeriesHebrides.mat');
timeAndSeries = timeAndSeriesHebrides;

% Perform stationary analysis
statEvaParams = tsEvaStationary(timeAndSeries, 'minPeakDistanceInDays', minPeakDistanceInDays);

% Compute and plot return levels
[rlevGEV, rlevGEVErr] = tsEvaComputeReturnLevelsGEVFromAnalysisObj(statEvaParams, [10, 20, 50, 100]);
hndl = tsEvaPlotReturnLevelsGEVFromAnalysisObj(statEvaParams, 1, 'ylim', [.5 1.5]);
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

**Code Structure:**
```matlab
% Trend-only analysis
[nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, timeWindow, 
    'transfType', 'trend', 'minPeakDistanceInDays', minPeakDistanceInDays);

% Plot series with trend
hndl = tsEvaPlotSeriesTrendStdDevFromAnalyisObj(nonStatEvaParams, statTransfData,
    'ylabel', 'Lvl (m)', 'title', seriesDescr);

% Plot 2D time-varying distributions
hndl = tsEvaPlotGEVImageScFromAnalysisObj(wr, nonStatEvaParams, statTransfData);
hndl = tsEvaPlotGPDImageScFromAnalysisObj(wr, nonStatEvaParams, statTransfData);

% Return levels at specific time
timeIndex = 1000;
hndl = tsEvaPlotReturnLevelsGEVFromAnalysisObj(nonStatEvaParams, timeIndex);
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

**Code Structure:**
```matlab
ciPercentile = 98;

[nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, timeWindow, 
    'transfType', 'trendCiPercentile', 'ciPercentile', ciPercentile, 
    'minPeakDistanceInDays', minPeakDistanceInDays);

% Compute return levels at specific time
timeIndex = 1000;
[rlevGPD, rlevGPDErr] = tsEvaComputeReturnLevelsGPDFromAnalysisObj(
    nonStatEvaParams, [10, 20, 50, 100], 'timeindex', timeIndex);
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

**Code Structure:**
```matlab
[nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, timeWindow, 
    'transfType', 'trendlinear', 'ciPercentile', ciPercentile, 
    'potPercentiles', [97:0.5:99], 'minPeakDistanceInDays', minPeakDistanceInDays);

% Extract shape parameter
epsilon = nonStatEvaParams(2).parameters.epsilon;

% Plot POT and return levels
tsPlotSeriesPotGPDRetLevFromAnalysisObj(nonStatEvaParams, statTransfData);

% Return levels at beginning and end
timeIndex = 1; % beginning
[rlevGPD, rlevGPDErr] = tsEvaComputeReturnLevelsGPDFromAnalysisObj(
    nonStatEvaParams, [5, 10, 30, 100], 'timeindex', timeIndex);
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
- `tsPlotSeriesPotGPDRetLevFromAnalysisObj()` - combined POT and return level plot

**Code Structure:**
```matlab
[nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, timeWindow, 
    'transfType', 'trendlinear', 'ciPercentile', ciPercentile, 
    'potPercentiles', [97:0.5:99], 'minPeakDistanceInDays', minPeakDistanceInDays);

% Compare return levels at different times
for lx=1:2
    if lx==1
        timeIndex = 1000;  % beginning
    else
        timeIndex = 240000;  % end
    end
    [rlevGPD, rlevGPDErr] = tsEvaComputeReturnLevelsGPDFromAnalysisObj(
        nonStatEvaParams, [5, 10, 30, 100], 'timeindex', timeIndex);
    hndl = tsEvaPlotReturnLevelsGPDFromAnalysisObj(nonStatEvaParams, timeIndex, 
        'ylim', rlRange, 'maxReturnPeriodYears', 200);
end
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
