# tsEVA 2.0 Monovariate Examples

This document provides comprehensive walkthroughs of monovariate Extreme Value Analysis (EVA) examples in tsEVA. Each example demonstrates practical applications with real-world data, explaining both the implementation steps and the scientific reasoning behind parameter choices.

## Overview of Monovariate EVA in tsEVA

Monovariate Extreme Value Analysis focuses on understanding the statistical behavior of extreme values in a single variable time series. The tsEVA toolbox provides two main approaches:

1. **GEV (Generalized Extreme Value)**: For block maxima (typically annual maxima)
2. **GPD (Generalized Pareto Distribution)**: For Peaks Over Threshold (POT) analysis

Both approaches can handle:
- **Stationary data**: Statistics remain constant over time
- **Non-stationary data**: Statistics vary with time (trend and/or seasonality)

### Key Concepts

- **Time Window**: The period over which moving statistics are computed (typically 5-30 years)
- **Transformation Types**:
  - `'trend'`: Accounts for long-term trends using moving standard deviation
  - `'trendCiPercentile'`: Uses moving percentiles (better for extremes, more uncertain)
  - `'trendLinear'`: Assumes linear trend in the data
  - `'seasonal'`: Includes seasonal cycles in addition to trends
  - `'seasonalCiPercentile'`: Seasonal analysis with percentile-based transformation
- **Return Levels**: Values expected to be exceeded with a given probability
- **Return Periods**: Average time between exceedances of a given level

---

## Example 1: Basic Non-Stationary GEV Analysis
### `exampleGenerateSeriesEVAGraphs.m`

**Purpose**: Comprehensive analysis of residual water levels with time-varying statistics and seasonality.

**Use Case**: Ocean water levels, wave heights, or any environmental variable with both long-term trends and seasonal patterns.

### Data Description

- **Dataset**: `timeAndSeriesHebrides.mat` - 30 years of residual water levels at the Hebrides islands
- **Format**: `timeAndSeries` = [timestamps, water levels (m)]
- **Characteristics**: Non-stationary with trends and seasonal variations

### Key Parameters

```matlab
timeWindow = 365.25*6;          % 6 years: minimum window for non-stationarity detection
minPeakDistanceInDays = 3;      % Ensure peak independence (3 days for water levels)
extremesRange = [.2 1.2];       % Range for plotting extremes (meters)
seasonalExtrRange = [.1 1.1];   % Range for seasonal extremes
```

**Why these values?**
- **6-year window**: Below this, series is considered quasi-stationary; above this, detects long-term changes
- **3-day separation**: Water level storms typically last 1-2 days; 3 days ensures independence
- **Ranges**: Set based on physical bounds of the phenomenon

### Step-by-Step Walkthrough

#### Part 1: Trend-Only Analysis

```matlab
% Load data
addpath('../');
load('timeAndSeriesHebrides.mat');
timeAndSeries = timeAndSeriesHebrides;

% Define parameters
timeWindow = 365.25*6;
minPeakDistanceInDays = 3;

% Perform non-stationary EVA with trend transformation
[nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, ...
    timeWindow, ...
    'transfType', 'trend', ...
    'minPeakDistanceInDays', minPeakDistanceInDays);
```

**What happens here?**
1. Time series is divided into overlapping windows
2. For each window: compute mean and standard deviation
3. Transform to stationarity by removing trend
4. Fit GEV and GPD to transformed extremes
5. Back-transform to original scale

#### Visualization: Series with Trend

```matlab
% Plot original series with detected trend and variance
hndl = tsEvaPlotSeriesTrendStdDevFromAnalyisObj(nonStatEvaParams, ...
    statTransfData, ...
    'ylabel', 'Lvl (m)', ...
    'title', 'Hebrides', ...
    'dateformat', 'yy', ...
    'xtick', tickTmStmp);
saveas(hndl{1}, 'seriesTrendOnly.png');
```

**Key Insight**: This plot shows the original data overlaid with:
- Moving mean (captures long-term trends)
- Moving standard deviation bands (captures changing variability)

#### Visualization: Time-Varying GEV and GPD

```matlab
% Define range for extreme values
wr = linspace(min(extremesRange), max(extremesRange), 1501);

% Plot 2D GEV (time on x-axis, level on y-axis, probability as color)
hndl = tsEvaPlotGEVImageScFromAnalysisObj(wr, nonStatEvaParams, ...
    statTransfData, ...
    'ylabel', 'Lvl (m)', ...
    'dateformat', 'yy', ...
    'xtick', tickTmStmp);
title('GEV');
saveas(hndl{1}, 'GEV2DTrendOnly.png');

% Plot 2D GPD
hndl = tsEvaPlotGPDImageScFromAnalysisObj(wr, nonStatEvaParams, ...
    statTransfData, ...
    'ylabel', 'Lvl (m)', ...
    'dateformat', 'yy', ...
    'xtick', tickTmStmp);
title('GPD');
saveas(hndl{1}, 'GPD2DTrendOnly.png');
```

**Interpretation**: 
- **GEV**: Shows probability distribution of annual maxima changing over time
- **GPD**: Shows probability distribution of threshold exceedances changing over time
- Color indicates probability density
- Horizontal slices = distribution at specific time
- Vertical slices = how a specific level's probability changes over time

#### Computing Return Levels

```matlab
% Choose a time instant (index 1000 in the time series)
timeIndex = 1000;
timeStamps = statTransfData.timeStamps;

% Plot GEV return levels at this time instant
hndl = tsEvaPlotReturnLevelsGEVFromAnalysisObj(nonStatEvaParams, ...
    timeIndex, ...
    'ylim', [.5 1.5]);
saveas(hndl{1}, 'GEV_ReturnLevels.png');

% Plot GPD return levels
hndl = tsEvaPlotReturnLevelsGPDFromAnalysisObj(nonStatEvaParams, ...
    timeIndex, ...
    'ylim', [.5 1.5]);
saveas(hndl{1}, 'GPD_ReturnLevels.png');
```

**Understanding Return Levels**:
- X-axis: Return period in years (e.g., 10, 50, 100)
- Y-axis: Water level (m)
- Shaded area: Confidence interval
- A 100-year return level of 1.2m means: this level is exceeded on average once every 100 years

**Note**: "For GEV the sample is small and the confidence interval is broad" - annual maxima give fewer data points than POT approach, leading to wider uncertainty.

#### Transformation Diagnostic

```matlab
% Plot the transformed stationary series
hndl = tsEvaPlotTransfToStatFromAnalysisObj(nonStatEvaParams, ...
    statTransfData, ...
    'dateformat', 'yy', ...
    'xtick', tickTmStmp);
saveas(hndl{1}, 'statSeriesTrendOnly.png');
```

**Purpose**: Verify transformation quality - transformed series should appear stationary (constant mean and variance).

#### Part 2: Seasonal Analysis

```matlab
% Perform non-stationary EVA with seasonal transformation
[nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, ...
    timeWindow, ...
    'transfType', 'seasonal', ...
    'minPeakDistanceInDays', minPeakDistanceInDays);

% Define new range for seasonal extremes
wr = linspace(min(seasonalExtrRange), max(seasonalExtrRange), 1501);
```

**Why seasonal analysis?**
- Water levels often show seasonal patterns (winter storms vs. summer calm)
- Seasonal transformation removes both trend AND annual cycle
- Provides more accurate extreme value estimates

```matlab
% Plot a time slice (1988-1993)
slice = {1988, 1993};

% Plot series for this period
hndl = tsEvaPlotSeriesTrendStdDevFromAnalyisObj(nonStatEvaParams, ...
    statTransfData, ...
    'ylabel', 'Lvl (m)', ...
    'title', '1988-1993', ...
    'minyear', slice{1}, ...
    'maxyear', slice{2});
saveas(hndl{1}, 'seriesSeasonal.png');

% Plot 3D GEV for this period
hndl = tsEvaPlotGEV3DFromAnalysisObj(wr, nonStatEvaParams, ...
    statTransfData, ...
    'xlabel', 'Lvl (m)', ...
    'minyear', slice{1}, ...
    'maxyear', slice{2});
title('GEV 3D, 1988-1993');
saveas(hndl{1}, 'GEV3DSeasonal.png');
```

**3D Plot Interpretation**:
- X-axis: Water level
- Y-axis: Time
- Z-axis: Probability density
- Shows how distribution evolves through time and seasons

### Expected Outputs

1. **seriesTrendOnly.png**: Original data with trend line and variability bands
2. **GEV2DTrendOnly.png**: Time-varying GEV distribution
3. **GPD2DTrendOnly.png**: Time-varying GPD distribution  
4. **GEV_ReturnLevels.png**: Return level plot for GEV at specific time
5. **GPD_ReturnLevels.png**: Return level plot for GPD at specific time
6. **statSeriesTrendOnly.png**: Transformed stationary series
7. **seriesSeasonal.png**: Series with seasonal component for 1988-1993
8. **GEV3DSeasonal.png**: 3D visualization of seasonal GEV
9. **GEV2DSeasonal.png**: 2D GEV with seasonality
10. **GPD2DSeasonal.png**: 2D GPD with seasonality

### Key Takeaways

1. **Two-stage approach**: First analyze trend-only, then add seasonality
2. **Always plot diagnostics**: Check that transformation produces stationary series
3. **GEV vs GPD tradeoff**: GEV has wider confidence intervals but is theoretically rigorous for maxima; GPD uses more data but requires threshold selection
4. **Time windows matter**: Too small = noisy estimates; too large = miss real changes
5. **Return levels are time-dependent**: In non-stationary analysis, a 100-year event today may be different from a 100-year event 20 years from now

---

## Example 2: Percentile-Based Confidence Intervals
### `exampleGenerateSeriesEVAGraphs_ciPercentile.m`

**Purpose**: Use moving percentiles instead of moving standard deviation to capture extreme value variability.

**Use Case**: When extremes show different variability patterns than the bulk of the data, or when data has heavy tails.

### Why Percentiles?

Standard deviation captures variability of the entire distribution, but extremes may behave differently. Using a high percentile (e.g., 98th) directly captures extreme value behavior.

**Trade-off**:
- ✅ Better models extreme variability
- ⚠️ Higher uncertainty (fewer data points at high percentiles)

### Key Parameters

```matlab
ciPercentile = 98;  % Use 98th percentile to capture extreme variability
```

**Choice of percentile**:
- Too low (e.g., 90%): Not capturing true extremes
- Too high (e.g., 99.5%): Too few points, high uncertainty
- Typical range: 95-99%

### Code Walkthrough

```matlab
% Trend-only analysis with percentile-based CI
[nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, ...
    timeWindow, ...
    'transfType', 'trendCiPercentile', ...
    'ciPercentile', ciPercentile, ...
    'minPeakDistanceInDays', minPeakDistanceInDays);
```

**What's different**: 
- `'transfType', 'trendCiPercentile'` instead of `'trend'`
- Must specify `'ciPercentile'` parameter
- Transformation uses moving 98th percentile instead of moving std dev

### Computing and Displaying Return Levels

```matlab
timeIndex = 1000;
timeStamps = statTransfData.timeStamps;

% Compute return levels explicitly
[rlevGEV, rlevGEVErr] = tsEvaComputeReturnLevelsGEVFromAnalysisObj(...
    nonStatEvaParams, [10, 20, 50, 100], ...
    'timeindex', timeIndex);
rlevGEV

% Plot with custom range
hndl = tsEvaPlotReturnLevelsGEVFromAnalysisObj(nonStatEvaParams, ...
    timeIndex, ...
    'ylim', rlRange);
title(['GEV return levels for ' datestr(tmstmpref)]);
saveas(hndl{1}, 'GEV_ReturnLevels_ciPercentile.png');
```

**Note**: This example also computes return levels explicitly before plotting, allowing you to inspect the numerical values.

### Seasonal Analysis with Percentiles

```matlab
% Seasonal analysis with percentile-based transformation
[nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, ...
    timeWindow, ...
    'transfType', 'seasonalCiPercentile', ...
    'ciPercentile', ciPercentile, ...
    'minPeakDistanceInDays', minPeakDistanceInDays);
```

Combines both seasonality AND percentile-based variability estimation.

### Expected Outputs

1. **seriesTrendOnly_ciPercentile.png**: Series with 98th percentile bands (instead of std dev)
2. **GEV2DTrendOnly_ciPercentile.png**: GEV with percentile transformation
3. **GEV_ReturnLevels_ciPercentile.png**: Return levels with percentile-based uncertainty
4. Similar plots for GPD and seasonal analysis

### Key Takeaways

1. **Choose transformation based on data**: If extremes show different variability patterns, use percentile approach
2. **Expect wider uncertainty**: Percentile-based methods have larger confidence intervals
3. **Compare approaches**: Run both standard deviation and percentile methods to understand sensitivity
4. **Percentile selection**: Test different percentiles (95, 98, 99) to assess robustness

---

## Example 3: Comparing Different CI Approaches  
### `exampleCompareDifferentCI.m`

**Purpose**: Systematically compare moving standard deviation vs. multiple percentile-based transformations.

**Use Case**: Determining which transformation approach is most appropriate for your data.

### Data Description

- **Dataset**: `timeAndSeries_waves_015_220E_055_509N.mat` - Wave height time series
- Can also use alternative wave datasets (commented in code)

### Key Parameters

```matlab
timeWindow = 30*365.25;  % 30 years (longer for wave data)
minPeakDistanceInDays = 3;
```

**Why 30 years?** Wave climate changes occur over longer timescales than water level variations.

### Step-by-Step Comparison

#### Method 1: Standard Deviation

```matlab
[nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, ...
    timeWindow, ...
    'minPeakDistanceInDays', minPeakDistanceInDays);

% Define plotting range based on data
minext = (max(statTransfData.nonStatSeries) + 3*min(statTransfData.nonStatSeries))/4;
maxext = max(statTransfData.nonStatSeries)*1.2;
xext = minext:.01:maxext;

% Plot GEV
tsEvaPlotGEVImageScFromAnalysisObj(xext, nonStatEvaParams, ...
    statTransfData, 'ylabel', 'Hs (m)');
title('standard deviation');

% Transformation diagnostic
tsEvaPlotTransfToStatFromAnalysisObj(nonStatEvaParams, statTransfData);
title('standard deviation');
```

#### Method 2: 98th Percentile

```matlab
[nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, ...
    timeWindow, ...
    'transfType', 'trendCiPercentile', ...
    'ciPercentile', 98, ...
    'minPeakDistanceInDays', minPeakDistanceInDays);

tsEvaPlotGEVImageScFromAnalysisObj(xext, nonStatEvaParams, ...
    statTransfData, 'ylabel', 'Hs (m)');
title('98th percentile');

tsEvaPlotTransfToStatFromAnalysisObj(nonStatEvaParams, statTransfData);
title('98th percentile');
```

#### Method 3: 98.5th Percentile

```matlab
[nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, ...
    timeWindow, ...
    'transfType', 'trendCiPercentile', ...
    'ciPercentile', 98.5, ...
    'minPeakDistanceInDays', minPeakDistanceInDays);

tsEvaPlotGEVImageScFromAnalysisObj(xext, nonStatEvaParams, ...
    statTransfData, 'ylabel', 'Hs (m)');
title('98.5th percentile');

tsEvaPlotTransfToStatFromAnalysisObj(nonStatEvaParams, statTransfData);
title('98.5th percentile');
```

#### Method 4: 99th Percentile

```matlab
[nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, ...
    timeWindow, ...
    'transfType', 'trendCiPercentile', ...
    'ciPercentile', 99, ...
    'minPeakDistanceInDays', minPeakDistanceInDays);

tsEvaPlotGEVImageScFromAnalysisObj(xext, nonStatEvaParams, ...
    statTransfData, 'ylabel', 'Hs (m)');
title('99th percentile');

tsEvaPlotTransfToStatFromAnalysisObj(nonStatEvaParams, statTransfData);
title('99th percentile');
```

### How to Compare Results

For each method, examine:

1. **GEV plots**: Are the patterns similar or different?
   - Similar patterns → results are robust
   - Very different patterns → data sensitive to transformation choice

2. **Transformation diagnostic plots**: Is the transformed series truly stationary?
   - Should show constant mean and variance
   - If not, transformation is inadequate

3. **Width of extreme value bands**: How do uncertainty estimates compare?
   - Percentile methods typically show wider bands
   - Very narrow bands may underestimate uncertainty

### Expected Outputs

- 4 pairs of plots (GEV + transformation diagnostic) for each method
- Visual comparison of how transformation choice affects results

### Key Takeaways

1. **No universal best method**: Choice depends on data characteristics
2. **Robustness check**: If results dramatically change with percentile choice, be cautious
3. **Document your choice**: Explain why you selected a particular transformation
4. **Consider physics**: Choose based on understanding of the phenomenon, not just statistics
5. **Standard deviation is safer**: If unsure, start with standard deviation (more stable)
6. **Percentiles for heavy tails**: Use percentile approach when data has very heavy tails or extreme outliers

---

## Example 4: GPD/POT for Drought Analysis
### `exampleSPISeries.m`

**Purpose**: Analyze Standardized Precipitation Index (SPI) using GPD only (no GEV), with long peak separation.

**Use Case**: Phenomena where:
- Events are rare and widely separated
- Annual maxima approach is inappropriate
- GPD/POT is more suitable than GEV

### Data Description

- **Dataset**: `timeAndSeries_SPI_179_750E_-16.750N.mat` - SPI time series
- **Variable**: SPI (Standardized Precipitation Index, inverted to analyze droughts as peaks)
- **Special characteristic**: Drought events are separated by months, not days

### Why SPI is Different

```matlab
% Invert series: analyze droughts (negative SPI) as positive peaks
timeAndSeries(:,2) = -timeAndSeries(:,2);
```

SPI is already normalized and standardized. Negative values indicate drought. By inverting, we analyze drought severity as if they were peaks.

### Key Parameters

```matlab
timeWindow = 50*315.25;         % 50 years: very long for climate analysis
minPeakDistanceInDays = 5*30.2; % ~5 months: droughts are long-lasting events
returnPeriodsInYears = [10 20 50 100];
```

**Critical choices**:
- **5-month separation**: Drought events typically last months; ensures independence
- **50-year window**: Climate trends emerge over multi-decadal timescales
- **GPD only**: With such long peak separation, "peaks per year" is very small; POT/GPD approach is more appropriate than annual maxima

### Code Walkthrough

```matlab
[nonStationaryEvaParams, stationaryTransformData, isValid] = tsEvaNonStationary(...
    timeAndSeries, ...
    timeWindow, ...
    'minPeakDistanceInDays', minPeakDistanceInDays, ...
    'transfType', 'trendCIPercentile', ...
    'cipercentile', 80, ...
    'potPercentiles', 80, ...  % Use 80th percentile as threshold
    'evdType', 'GPD');          % ONLY fit GPD, not GEV
```

**Key parameters explained**:
- `'evdType', 'GPD'`: Skip GEV fitting entirely
- `'potPercentiles', 80`: Use only one threshold (80th percentile) for POT
- `'cipercentile', 80`: Match transformation percentile to POT threshold

**Why 80th percentile?**
- Lower than typical (95-98%) because SPI is already normalized
- Captures sufficient drought events for GPD fitting
- Higher values would give too few events

### Visualization

```matlab
tsEvaPlotSeriesTrendStdDevFromAnalyisObj(nonStationaryEvaParams, ...
    stationaryTransformData, ...
    'plotpercentile', 95., ...
    'ylabel', '-SPI', ...
    'legendLocation', 'southwest');
```

**Note**: `'plotpercentile', 95.` adds 95th percentile line to plot for reference, independent of the 80th percentile used in analysis.

### Computing Return Levels

```matlab
[returnLevels, returnLevelsErr] = tsEvaComputeReturnLevelsGPDFromAnalysisObj(...
    nonStationaryEvaParams, returnPeriodsInYears);

% Invert back to original SPI scale
returnLevels = returnLevels * -1
```

**Interpretation**: A 100-year return level of -2.5 means a drought of SPI = -2.5 is expected once per century.

### Expected Outputs

- Series plot showing inverted SPI with trend and percentile bands
- Return level values for 10, 20, 50, 100-year drought events

### Key Takeaways

1. **Adapt method to data**: Long-separated events need different parameters
2. **GPD without GEV**: Valid when annual maxima concept doesn't apply
3. **Single threshold**: Can use just one percentile for POT (not the typical multi-percentile exploration)
4. **Transform for interpretation**: Remember to invert results back to original scale
5. **Climate timescales**: Use very long time windows (50+ years) for climate data

---

## Example 5: GEV-Only Analysis for Temperature Extremes
### `exampleTASMaxSeries.m`

**Purpose**: Analyze yearly maximum temperatures (heat waves) using GEV exclusively.

**Use Case**: When data consists of one value per year (annual maxima) - GPD/POT is meaningless.

### Data Description

- **Dataset**: `timeAndSeriesTASMax.mat` - Annual maximum air surface temperature
- **Format**: One value per year (1850-2100 or similar)
- **Application**: Climate change impact on heat waves

### Why GEV Only?

```matlab
'evdType', 'GEV'  % Only fit GEV, skip GPD
```

When you have exactly one value per year (annual maxima):
- GEV is theoretically appropriate (models maxima directly)
- GPD requires threshold exceedances within years → not applicable here
- No need to select peaks or thresholds

### Key Parameters

```matlab
timeWindow = 50*315.25;         % 50 years
minPeakDistanceInDays = 5*30.2; % Not really used (already annual maxima)
returnPeriodsInYears = [20 50 100 300];
```

### Code Walkthrough

```matlab
[nonStationaryEvaParams, stationaryTransformData, isValid] = tsEvaNonStationary(...
    timeAndSeries, ...
    timeWindow, ...
    'minPeakDistanceInDays', minPeakDistanceInDays, ...
    'extremeLowThreshold', .1, ...  % Exclude very low values (data quality)
    'evdType', 'GEV');              % GEV only
```

**New parameter**: `'extremeLowThreshold', .1` - excludes the lowest 10% of values, useful for filtering erroneous or missing data.

### Visualization: Time Series

```matlab
tsEvaPlotSeriesTrendStdDevFromAnalyisObj(nonStationaryEvaParams, ...
    stationaryTransformData, ...
    'ylabel', 'TAS', ...
    'legendLocation', 'northwest');
ylim([0 40]);
text(datenum(2060, 1, 1), 35, 'Series and trends', 'fontsize', 25);
```

Shows temperature data with trends - clear for climate change visualization.

### Visualization: Time-Varying GEV

```matlab
tsEvaPlotGEVImageScFromAnalysisObj((0:.001:40)', ...
    nonStationaryEvaParams, ...
    stationaryTransformData, ...
    'ylabel', 'TAS');
text(datenum(1980, 1, 1), 35, 'Time varying GEV', 'fontsize', 30);
```

Shows how the probability of extreme heat changes over time (climate change impact).

### Return Levels at Different Time Points

```matlab
% Early period (1995)
timeIndex = 26;
rlRange = [0 14];
hndl = tsEvaPlotReturnLevelsGEVFromAnalysisObj(nonStationaryEvaParams, ...
    timeIndex, ...
    'ylim', rlRange, ...
    'ylabel', 'return levels (TAS)');
ax = gca;
ax.YTick = 0:2:16;
text(7, 12.5, 'Return level 1995', 'fontsize', 30);

% Late period (2095)  
timeIndex = size(timeAndSeries, 1) - 4;
rlRange = [0 70];
hndl = tsEvaPlotReturnLevelsGEVFromAnalysisObj(nonStationaryEvaParams, ...
    timeIndex, ...
    'ylim', rlRange, ...
    'ylabel', 'return levels (TAS)');
text(7, 65, 'Return level 2095', 'fontsize', 30);
```

**Key insight**: Note the dramatic difference in y-axis ranges:
- 1995: 0-14°C
- 2095: 0-70°C  

This illustrates the increasing severity of heat wave extremes under climate change.

### Computing Return Levels

```matlab
[returnLevels, returnLevelsErr] = tsEvaComputeReturnLevelsGEVFromAnalysisObj(...
    nonStationaryEvaParams, returnPeriodsInYears);
```

Provides numerical values for 20, 50, 100, 300-year heat wave temperatures.

### Expected Outputs

1. Series plot showing temperature trends
2. Time-varying GEV distribution image
3. Return level plots for 1995 and 2095 showing dramatic increase
4. Numerical return level values

### Key Takeaways

1. **Annual maxima = GEV**: When data is already annual maxima, use GEV exclusively
2. **Skip POT**: No need for threshold selection or peak finding
3. **Climate change visualization**: Time-varying return levels clearly show changing risk
4. **Long time windows**: Essential for capturing climate trends (50+ years)
5. **Projections**: Can analyze future projections if data includes model outputs

---

## Example 6: Stationary Analysis  
### `exampleEVAStationary.m`

**Purpose**: Apply traditional stationary EVA when time-varying statistics are not needed.

**Use Case**: 
- Short time series where non-stationarity can't be reliably estimated
- Exploratory analysis before attempting non-stationary approach
- Data with no apparent trends

### Key Difference

```matlab
statEvaParams = tsEvaStationary(timeAndSeries, ...
    'minPeakDistanceInDays', minPeakDistanceInDays);
```

Uses `tsEvaStationary` instead of `tsEvaNonStationary` - fits a single GEV and GPD to entire series.

### Code Walkthrough

#### Basic Stationary Analysis

```matlab
% Fit stationary GEV and GPD
statEvaParams = tsEvaStationary(timeAndSeries, ...
    'minPeakDistanceInDays', minPeakDistanceInDays);

% Compute return levels (single set of values, not time-varying)
[rlevGEV, rlevGEVErr] = tsEvaComputeReturnLevelsGEVFromAnalysisObj(...
    statEvaParams, [10, 20, 50, 100]);
rlevGEV

% Plot return levels (timeIndex=1 because there's only one time point)
hndl = tsEvaPlotReturnLevelsGEVFromAnalysisObj(statEvaParams, 1, ...
    'ylim', [.5 1.5]);
title('GEV');
saveas(hndl{1}, 'GEV_ReturnLevels_STATIONARY.png');
```

#### Fixed Threshold POT

```matlab
% Manually specify POT threshold
potThreshold = prctile(timeAndSeries(:,2), 98);

statEvaParams = tsEvaStationary(timeAndSeries, ...
    'minPeakDistanceInDays', minPeakDistanceInDays, ...
    'doSampleData', false, ...      % Skip automatic threshold selection
    'potThreshold', potThreshold);  % Use fixed threshold

[rlevGPD, rlevGPDErr] = tsEvaComputeReturnLevelsGPDFromAnalysisObj(...
    statEvaParams, [10, 20, 50, 100]);
rlevGPD
```

**Why fixed threshold?**
- Gives full control over threshold selection
- Useful for sensitivity analysis
- Can align with physical thresholds (e.g., flood stage)

#### Gumbel Distribution

```matlab
% Fit Gumbel (special case of GEV with shape parameter = 0)
statEvaParams = tsEvaStationary(timeAndSeries, ...
    'minPeakDistanceInDays', minPeakDistanceInDays, ...
    'gevtype', 'gumbel');

[rlevGEV, rlevGEVErr] = tsEvaComputeReturnLevelsGEVFromAnalysisObj(...
    statEvaParams, [10, 20, 50, 100]);
rlevGEV

hndl = tsEvaPlotReturnLevelsGEVFromAnalysisObj(statEvaParams, 1, ...
    'ylim', [.5 1.5]);
title('Gumbel');
```

**When to use Gumbel?**
- When shape parameter is close to zero
- For simplicity (fewer parameters to estimate)
- When data suggests exponential tail behavior

### Expected Outputs

- Return level plots showing single curve with confidence intervals
- Numerical return level values for GEV, GPD, and Gumbel

### Key Takeaways

1. **Simpler but limited**: Stationary analysis is easier but assumes no temporal changes
2. **Useful for validation**: Compare stationary vs non-stationary to assess if complexity is justified
3. **Fixed thresholds**: Give more control in POT analysis
4. **Gumbel option**: Reduces parameters when appropriate
5. **Return levels constant**: Unlike non-stationary case, return levels don't vary with time

---

## Example 7: Linear Trend Analysis
### `exampleGenerateSeriesEVAGraphs_trendLinear.m`

**Purpose**: Analyze extremes when a linear trend is appropriate (common in climate change studies).

**Use Case**: Data showing steady increase/decrease over time without complex variations.

### Key Parameters

```matlab
timeWindow = 365.25*30;  % 30 years
ciPercentile = 99;       % 99th percentile for transformation
```

### Code Walkthrough

```matlab
[nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, ...
    timeWindow, ...
    'transfType', 'trendlinear', ...
    'ciPercentile', ciPercentile, ...
    'potPercentiles', [97:0.5:99], ...  % Multiple thresholds for robustness
    'minPeakDistanceInDays', minPeakDistanceInDays);
```

**What's different**:
- `'transfType', 'trendlinear'`: Assumes linear change in mean
- `'potPercentiles', [97:0.5:99]`: Tests multiple thresholds (97%, 97.5%, 98%, ..., 99%)

### Extracting Statistics

```matlab
epsilon = nonStatEvaParams(2).parameters.epsilon;  % GPD shape parameter
pvalue = statTransfData.pValueChange;              % Significance of trend
```

Useful for assessing:
- **epsilon**: Tail behavior (negative = bounded, zero = exponential, positive = heavy)
- **pvalue**: Is the detected trend statistically significant?

### Multiple Time Point Analysis

```matlab
% Beginning of series
timeIndex = 1;
[rlevGPD, rlevGPDErr] = tsEvaComputeReturnLevelsGPDFromAnalysisObj(...
    nonStatEvaParams, [5, 10, 30, 100], 'timeindex', timeIndex);
tsEvaPlotReturnLevelsGPDFromAnalysisObj(nonStatEvaParams, timeIndex, ...
    'ylim', rlRange, 'maxReturnPeriodYears', 200);
title('GPD return levels-begin');

% End of series  
timeIndex = size(timeAndSeries, 1);
[rlevGPD, rlevGPDErr] = tsEvaComputeReturnLevelsGPDFromAnalysisObj(...
    nonStatEvaParams, [5, 10, 30, 100], 'timeindex', timeIndex);
tsEvaPlotReturnLevelsGPDFromAnalysisObj(nonStatEvaParams, timeIndex, ...
    'ylim', rlRange, 'maxReturnPeriodYears', 200);
title('GPD return levels-end');
```

Directly compares how return levels have changed from beginning to end of time series.

### Visualization of POT

```matlab
tsPlotSeriesPotGPDRetLevFromAnalysisObj(nonStatEvaParams, statTransfData);
title('PotAndReturnLevelsLinearTrend');
```

Shows original series with:
- POT threshold(s)
- Selected peaks
- Return level curves

### Expected Outputs

1. Series plot with linear trend
2. POT visualization with threshold and peaks
3. 2D GPD distribution over time
4. Return level plots for beginning and end of series
5. Numerical epsilon and p-value for trend assessment

### Key Takeaways

1. **Linear trend assumption**: Simpler than general trend but may not capture complex variations
2. **Multiple thresholds**: Using a range of percentiles improves robustness
3. **Shape parameter**: Check epsilon to understand tail behavior
4. **Statistical significance**: Use p-value to validate that trend is real
5. **Temporal comparison**: Plot return levels at multiple times to quantify change

---

## Example 8: GPD with Negative Shape Parameter
### `exampleGenerateSeriesEVAGraphs_gpdNegShapeParam.m`

**Purpose**: Handle cases where GPD shape parameter is negative (bounded upper tail).

**Use Case**: Variables with physical upper bounds or when fitting suggests negative shape.

### Key Parameters

```matlab
timeWindow = 365.25*15;         % 15 years
minPeakDistanceInDays = 14;     % 2 weeks (longer than typical)
ciPercentile = 99;
potPercentiles = [97:0.5:99];   % Multiple thresholds
```

### Code Walkthrough

```matlab
[nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, ...
    timeWindow, ...
    'transfType', 'trendlinear', ...
    'ciPercentile', ciPercentile, ...
    'potPercentiles', [97:0.5:99], ...
    'minPeakDistanceInDays', minPeakDistanceInDays);

% Extract shape parameter
epsilon = nonStatEvaParams(2).parameters.epsilon;
```

**Negative epsilon interpretation**:
- Distribution has an upper bound
- Upper bound = location + scale/abs(epsilon)
- Return levels eventually plateau
- Extrapolation beyond observed range is unreliable

### Visualization with Shape Parameter

```matlab
tsEvaPlotReturnLevelsGPDFromAnalysisObj(nonStatEvaParams, timeIndex, ...
    'ylim', rlRange, 'maxReturnPeriodYears', 200);
hold on;
legend(sprintf('negative shapeparemeter = %.3e', epsilon), ...
    'Location', 'northwest');
title('GPD return levels-begin');
```

Annotates plot with actual shape parameter value for transparency.

### Comparing Beginning and End

```matlab
% Beginning (timeIndex = 1)
epsilon_begin = nonStatEvaParams(2).parameters.epsilon;
[rlevGPD_begin, rlevGPDErr_begin] = tsEvaComputeReturnLevelsGPDFromAnalysisObj(...
    nonStatEvaParams, [5, 10, 30, 100], 'timeindex', 1);

% End (timeIndex = last)  
epsilon_end = nonStatEvaParams(2).parameters.epsilon;
[rlevGPD_end, rlevGPDErr_end] = tsEvaComputeReturnLevelsGPDFromAnalysisObj(...
    nonStatEvaParams, [5, 10, 30, 100], 'timeindex', timeIndex);
```

If negative shape parameter, check if:
- Upper bound is changing over time
- Bound is consistent with physical constraints

### Expected Outputs

1. Series plot showing bounded behavior
2. POT and return level visualization
3. 2D GPD showing upper bound constraint
4. Return level plots with shape parameter annotated
5. Comparison of bounds at beginning and end

### Key Takeaways

1. **Negative shape is valid**: Not an error, indicates bounded distribution
2. **Physical interpretation**: Makes sense for variables with natural upper limits
3. **Extrapolation caution**: Don't extrapolate far beyond observed data with negative shape
4. **Upper bound formula**: Can compute exact upper bound from parameters
5. **Time-varying bounds**: In non-stationary case, upper bound may change over time

---

## Example 9: Gumbel for SPI
### `exampleSPISeries_Gumbel.m`

**Purpose**: Use Gumbel distribution (GEV with shape=0) for SPI drought analysis.

**Use Case**: When shape parameter is close to zero or when simpler model is preferred.

### Key Difference from Example 4

```matlab
[nonStationaryEvaParams, stationaryTransformData, isValid] = tsEvaNonStationary(...
    timeAndSeries, ...
    timeWindow, ...
    'minPeakDistanceInDays', minPeakDistanceInDays, ...
    'transfType', 'trendCIPercentile', ...
    'cipercentile', 80, ...
    'gevType', 'gumbel', ...  % Use Gumbel instead of full GEV
    'evdType', 'GEV');         % Fit GEV (Gumbel variant), not GPD
```

**Gumbel vs. full GEV**:
- Gumbel: 2 parameters (location, scale)
- GEV: 3 parameters (location, scale, shape)
- Gumbel: Simpler, more stable with limited data
- GEV: More flexible, can capture heavier or lighter tails

### When to Choose Gumbel

1. **Limited data**: Fewer parameters = more stable estimates
2. **Shape ≈ 0**: If full GEV gives shape parameter near zero, Gumbel is appropriate
3. **Theoretical grounds**: Some phenomena have exponential tail behavior
4. **Simplicity**: When interpretability is important

### Code Walkthrough

```matlab
% Invert SPI for drought analysis
timeAndSeries(:,2) = -timeAndSeries(:,2);

% Fit Gumbel with percentile transformation
[nonStationaryEvaParams, stationaryTransformData, isValid] = tsEvaNonStationary(...
    timeAndSeries, timeWindow, ...
    'minPeakDistanceInDays', 5*30.2, ...
    'transfType', 'trendCIPercentile', ...
    'cipercentile', 80, ...
    'gevType', 'gumbel', ...
    'evdType', 'GEV');

% Plot series
tsEvaPlotSeriesTrendStdDevFromAnalyisObj(nonStationaryEvaParams, ...
    stationaryTransformData, ...
    'plotpercentile', 95., ...
    'ylabel', '-SPI', ...
    'legendLocation', 'southwest');

% Compute return levels (Gumbel-based)
[returnLevels, returnLevelsErr] = tsEvaComputeReturnLevelsGEVFromAnalysisObj(...
    nonStationaryEvaParams, returnPeriodsInYears);

% Invert back to SPI scale
returnLevels = returnLevels * -1
```

### Expected Outputs

- Series plot with Gumbel-based transformation
- Return level values for drought events (inverted SPI)

### Key Takeaways

1. **Gumbel = simplified GEV**: Appropriate when shape parameter is zero
2. **Stability**: More stable parameter estimation with limited data
3. **Compare with full GEV**: Check if simplification is justified
4. **SPI application**: Particularly relevant for normalized indices like SPI
5. **Use GEV functions**: Even though it's Gumbel, use GEV computation/plotting functions

---

## Common Patterns Across Examples

### 1. Data Loading and Setup

```matlab
addpath('../');                          % Add tsEVA functions to path
load('datafile.mat');                    % Load time series data
timeAndSeries = loadedVariable;          % Ensure format: [time, values]
```

### 2. Parameter Selection Logic

| Parameter | Typical Range | Selection Criteria |
|-----------|---------------|-------------------|
| `timeWindow` | 5-50 years | Longer for climate data, shorter for engineering |
| `minPeakDistanceInDays` | 1-150 days | Based on event duration and independence |
| `ciPercentile` | 95-99% | Higher for better extreme capture, but more uncertainty |
| `potPercentiles` | 95-99% or range | Single value or range for robustness testing |
| `evdType` | 'GEV', 'GPD', or both | GEV for annual maxima, GPD for POT, both for comparison |

### 3. Analysis Flow

```matlab
% 1. Non-stationary analysis
[nonStatEvaParams, statTransfData] = tsEvaNonStationary(...);

% 2. Plot series with trends
tsEvaPlotSeriesTrendStdDevFromAnalyisObj(...);

% 3. Visualize distributions  
tsEvaPlotGEVImageScFromAnalysisObj(...);  % and/or GPD

% 4. Compute return levels
[returnLevels, returnLevelsErr] = tsEvaComputeReturnLevelsGEVFromAnalysisObj(...);

% 5. Plot return levels
tsEvaPlotReturnLevelsGEVFromAnalysisObj(...);

% 6. Check transformation
tsEvaPlotTransfToStatFromAnalysisObj(...);
```

### 4. Scientific Reasoning

**Time Window Selection**:
- Too small → Noisy estimates, can't detect long-term trends
- Too large → Over-smoothing, miss real variations
- Rule of thumb: 5-10 years for engineering, 30-50 years for climate

**Peak Distance Selection**:
- Must ensure independence of extreme events
- Typical event duration + safety margin
- Examples: 3 days for water levels, 14 days for large storms, 150 days for droughts

**Transformation Type**:
- Start with `'trend'` (standard deviation based)
- Use `'trendCiPercentile'` if extremes show different variability pattern
- Add seasonality (`'seasonal'`) if clear annual cycle
- Use `'trendLinear'` for simple monotonic trends

**GEV vs GPD**:
- GEV: Theoretically sound for annual maxima, but fewer data points
- GPD: Uses more data (all threshold exceedances), but requires threshold selection
- Both: Run both and compare for robustness

### 5. Diagnostic Checks

Always examine:

1. **Transformation plot**: Transformed series should be stationary
2. **Shape parameter**: Check if reasonable (negative → bounded, positive → heavy tails)
3. **P-values**: Assess statistical significance of trends
4. **Confidence intervals**: Wide intervals indicate high uncertainty
5. **Sensitivity**: Test different parameters to assess robustness

### 6. Common Pitfalls

❌ **Avoid**:
- Time window smaller than parameter of interest (e.g., 5-year window for 100-year return level)
- Ignoring seasonality in seasonal data
- Over-interpreting results with wide confidence intervals
- Extrapolating far beyond observed data range
- Using GEV for threshold exceedance or GPD for annual maxima

✅ **Best Practices**:
- Plot diagnostics before believing results
- Compare multiple approaches (GEV vs GPD, std dev vs percentile)
- Document all parameter choices and reasoning
- Report confidence intervals, not just point estimates
- Validate with physical understanding of the phenomenon

---

## Summary Table of Examples

| Example | Transformation | Distribution | Data Type | Key Learning |
|---------|---------------|--------------|-----------|--------------|
| 1. Basic Non-Stationary | Trend + Seasonal | GEV + GPD | Water levels | Complete workflow, both approaches |
| 2. Percentile CI | Trend/Seasonal CI Percentile | GEV + GPD | Water levels | Alternative variability estimation |
| 3. Compare CI | Multiple methods | GEV + GPD | Wave heights | Systematic comparison of methods |
| 4. SPI/GPD | Trend CI Percentile | GPD only | Drought index | Long peaks, GPD-only approach |
| 5. TAS/GEV | Trend | GEV only | Temperature | Annual maxima, GEV-only approach |
| 6. Stationary | None (stationary) | GEV + GPD + Gumbel | Water levels | Traditional EVA, no time-variation |
| 7. Linear Trend | Trend Linear | GEV + GPD | Various | Linear trend assumption, multiple thresholds |
| 8. Negative Shape | Trend Linear | GPD with negative ε | Various | Bounded distributions |
| 9. Gumbel SPI | Trend CI Percentile | Gumbel (GEV) | Drought index | Simplified GEV model |

---

## Next Steps

After working through these examples:

1. **Apply to your data**: Adapt the examples to your specific dataset
2. **Explore parameter sensitivity**: Test different time windows, peak distances, percentiles
3. **Compare methods**: Run multiple approaches and understand differences
4. **Validate results**: Use physical reasoning and domain knowledge
5. **Read copula examples**: For multivariate extremes, see copula documentation
6. **Consult function reference**: See `2_Function_Reference.md` for detailed parameter descriptions

---

## References

For theoretical background on these methods, see:
- `1_Core_Methodology.md` - Mathematical foundations
- `3_Workflow_Patterns.md` - Standard analysis patterns
- Coles, S. (2001). *An Introduction to Statistical Modeling of Extreme Values*. Springer.
- Extreme value theory literature for GEV and GPD distributions

---

*Document version: 1.0*  
*Last updated: 2024*  
*Part of tsEVA 2.0 documentation suite*
