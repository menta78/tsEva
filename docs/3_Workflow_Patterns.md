# tsEVA 2.0 Workflow Patterns

This document describes standard analysis workflows and common patterns for using tsEVA 2.0. All code examples use only documented functions.

## Standard Workflows

### Workflow 1: Basic Non-Stationary GEV Analysis

**Use Case**: Analyzing annual maxima with time-varying statistics (trend only, no seasonality)

```matlab
% 1. Load data
load('timeAndSeriesData.mat');  % Assumes variable timeAndSeries [time, values]

% 2. Define parameters
timeWindow = 365.25 * 6;  % 6 years in days
minPeakDistanceInDays = 3;  % Minimum separation between peaks

% 3. Perform non-stationary EVA
[nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, ...
    timeWindow, ...
    'transfType', 'trend', ...
    'minPeakDistanceInDays', minPeakDistanceInDays);

% 4. Compute return levels
returnPeriods = [10, 50, 100, 500];  % Years
timeInstants = linspace(min(timeAndSeries(:,1)), max(timeAndSeries(:,1)), 100);
[rl] = tsEvaComputeReturnLevelsGEVFromAnalysisObj(returnPeriods, ...
    timeInstants, nonStatEvaParams, statTransfData);

% 5. Visualize results
tsEvaPlotGEVImageScFromAnalysisObj(wr, nonStatEvaParams, statTransfData);
tsEvaPlotReturnLevelsGEVFromAnalysisObj(returnPeriods, nonStatEvaParams, statTransfData);
```

**Key Pattern**: `NonStationary → ComputeReturnLevels → Plot`

---

### Workflow 2: Seasonal Non-Stationary Analysis  

**Use Case**: Data with clear seasonal cycle (e.g., sea levels, temperature)

```matlab
% 1. Load data
load('timeAndSeriesData.mat');

% 2. Define parameters
timeWindow = 365.25 * 6;  % 6 years
minPeakDistanceInDays = 3;

% 3. Non-stationary EVA with seasonality
[nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, ...
    timeWindow, ...
    'transfType', 'seasonal', ...  % KEY: Include seasonality
    'minPeakDistanceInDays', minPeakDistanceInDays);

% 4. Plot seasonal component
tsEvaPlotSeasonalityGevFromAnalysisObj(nonStatEvaParams, statTransfData);

% 5. Compute and plot return levels
returnPeriods = [10, 50, 100];
tsEvaPlotReturnLevelsGEVFromAnalysisObj(returnPeriods, nonStatEvaParams, statTransfData);
```

**Key Pattern**: Use `'transfType', 'seasonal'` for seasonal data

---

### Workflow 3: GPD (Peaks Over Threshold)

**Use Case**: When events are not limited to annual maxima (e.g., all storms above threshold)

```matlab
% 1. Load data
load('timeAndSeriesData.mat');

% 2. Define parameters  
timeWindow = 365.25 * 50;  % 50 years for long-term trends
minPeakDistanceInDays = 60;  % Longer separation for independence
thresholdQuantile = 0.95;  % Use 95th percentile as threshold

% 3. Non-stationary EVA (GPD approach)
[nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, ...
    timeWindow, ...
    'transfType', 'trend', ...
    'minPeakDistanceInDays', minPeakDistanceInDays, ...
    'thresholdQuantile', thresholdQuantile);

% 4. Compute GPD return levels
returnPeriods = [10, 50, 100, 500];
timeInstants = linspace(min(timeAndSeries(:,1)), max(timeAndSeries(:,1)), 50);
[rl] = tsEvaComputeReturnLevelsGPDFromAnalysisObj(returnPeriods, ...
    timeInstants, nonStatEvaParams, statTransfData);

% 5. Visualize GPD results
tsEvaPlotGPDImageScFromAnalysisObj(wr, nonStatEvaParams, statTransfData);
tsEvaPlotReturnLevelsGPDFromAnalysisObj(returnPeriods, nonStatEvaParams, statTransfData);
```

**Key Pattern**: Set `thresholdQuantile` for POT approach, use GPD-specific plotting functions

---

### Workflow 4: Percentile-Based Confidence Interval

**Use Case**: When extremes change more than the overall distribution (more sensitive to tail behavior)

```matlab
% 1. Load data  
load('timeAndSeriesData.mat');

% 2. Define parameters
timeWindow = 365.25 * 6;
minPeakDistanceInDays = 3;
ciPercentile = 98.5;  % Use 98.5th percentile for CI estimation

% 3. Transformation using percentile-based CI
[stSeries] = tsEvaTransformSeriesToStat_ciPercentile(timeAndSeries, ...
    timeWindow, ciPercentile);

% 4. Perform stationary EVA on transformed series
[evaParams] = tsEvaStationary(stSeries, 'minPeakDistanceInDays', minPeakDistanceInDays);

% 5. Back-transform and compute return levels
% (Manual back-transformation required - consult exampleGenerateSeriesEVAGraphs_ciPercentile.m)
```

**Key Pattern**: Use `tsEvaTransformSeriesToStat_ciPercentile` when extremes are changing faster than mean

---

### Workflow 5: Stationary EVA (Baseline)

**Use Case**: Data without significant trends (or for comparison with non-stationary approach)

```matlab
% 1. Load data
load('timeAndSeriesData.mat');

% 2. Stationary EVA  
[evaParams] = tsEvaStationary(timeAndSeries, ...
    'minPeakDistanceInDays', 3);

% 3. Compute return levels (stationary)
returnPeriods = [10, 50, 100, 500];
[rl] = tsEvaComputeReturnLevelsGEV(returnPeriods, evaParams.gevParams);

% 4. Plot (stationary version)
tsEvaPlotReturnLevelsGEVStationary(returnPeriods, evaParams);
```

**Key Pattern**: Use `tsEvaStationary` instead of `tsEvaNonStationary` for baseline analysis

---

### Workflow 6: Bivariate Copula Analysis (Non-Stationary)

**Use Case**: Joint analysis of two related variables (e.g., wave height + storm surge)

```matlab
% 1. Load data
load('bivariateData.mat');  % timeAndSeries1, timeAndSeries2

% 2. Define parameters
timeWindow = 365.25 * 40;  % 40 years
minPeakDistanceInDays = 30;
thresholdQuantile1 = 0.95;  % 95th percentile for variable 1  
thresholdQuantile2 = 0.99;  % 99th percentile for variable 2
maxTimeLagInDays = 45;  % Maximum lag for joint occurrence

% 3. Sample joint peaks
[peaks] = tsCopulaSampleJointPeaksMultiVariatePruning(...
    {timeAndSeries1, timeAndSeries2}, ...
    minPeakDistanceInDays, ...
    [thresholdQuantile1, thresholdQuantile2], ...
    maxTimeLagInDays);

% 4. Fit copula extremes
copulaFamily = 'Gumbel';  % Or 'Gaussian' or 'Frank' (bivariate only)
[copulaParams, marginalParams] = tsCopulaExtremes(peaks, timeWindow, ...
    'copulaFamily', copulaFamily, ...
    'minPeakDistanceInDays', minPeakDistanceInDays);

% 5. Compute bivariate return periods
returnPeriods = [10, 50, 100];
[rpMatrix] = tsCopulaComputeBivarRP(returnPeriods, copulaParams, marginalParams);

% 6. Visualize
tsCopulaPlotBivariate(copulaParams, marginalParams);
tsPlotBivarReturnPeriod(rpMatrix, marginalParams);
```

**Key Pattern**: 
1. Sample joint peaks with `tsCopulaSampleJointPeaksMultiVariatePruning`
2. Fit with `tsCopulaExtremes`
3. Compute return periods with `tsCopulaComputeBivarRP`
4. Visualize with copula-specific plotting functions

---

### Workflow 7: Trivariate Copula Analysis

**Use Case**: Joint analysis of three variables (e.g., significant wave height, peak period, direction)

```matlab
% 1. Load data  
load('trivariateData.mat');  % timeAndSeries1, timeAndSeries2, timeAndSeries3

% 2. Define parameters
timeWindow = 365.25 * 40;
minPeakDistanceInDays = 30;
thresholds = [0.95, 0.95, 0.90];  % Quantile thresholds for each variable
maxTimeLagInDays = 45;

% 3. Sample joint peaks (trivariate)
[peaks] = tsCopulaSampleJointPeaksMultiVariatePruning(...
    {timeAndSeries1, timeAndSeries2, timeAndSeries3}, ...
    minPeakDistanceInDays, ...
    thresholds, ...
    maxTimeLagInDays);

% 4. Fit copula extremes (Gaussian or Gumbel for multivariate)
copulaFamily = 'Gumbel';  % MUST be Gaussian or Gumbel (not Frank)
[copulaParams, marginalParams] = tsCopulaExtremes(peaks, timeWindow, ...
    'copulaFamily', copulaFamily, ...
    'minPeakDistanceInDays', minPeakDistanceInDays);

% 5. Visualize trivariate results
tsCopulaPlotTrivariate(copulaParams, marginalParams);
```

**Key Pattern**: 
- Use Gaussian or Gumbel copula (NOT Frank) for trivariate analysis
- Adjust thresholds to ensure sufficient joint events (50+ recommended)

---

### Workflow 8: Year Extremes Copula (Annual Maxima Joint Analysis)

**Use Case**: Analyzing joint distribution of annual maxima across multiple variables

```matlab
% 1. Load data
load('bivariateData.mat');

% 2. Compute annual maxima for each variable
annualMaxima1 = tsEvaComputeAnnualMaxima(timeAndSeries1);
annualMaxima2 = tsEvaComputeAnnualMaxima(timeAndSeries2);

% 3. Fit year extremes copula  
copulaFamily = 'Gumbel';
[copulaFit] = tsCopulaYearExtrFit({annualMaxima1, annualMaxima2}, copulaFamily);

% 4. Generate random samples from fitted copula
nSamples = 1000;
[samples] = tsCopulaYearExtrRnd(copulaFit, nSamples);

% 5. Visualize
tsCopulaYearExtrPlotSctrBivar(copulaFit);
```

**Key Pattern**: Use `YearExtr` functions for annual maxima joint analysis

---

## Common Patterns

### Pattern A: Data Preparation

```matlab
% Load and inspect
load('rawData.mat');
timeAndSeries = [time(:), values(:)];  % Ensure column vectors

% Remove NaNs or fill gaps
timeAndSeries = timeAndSeries(~isnan(timeAndSeries(:,2)), :);
% OR
[filledSeries] = tsEvaFillSeries(timeAndSeries);

% Check time step
dt = tsEvaGetTimeStep(timeAndSeries);
fprintf('Time step: %.2f days\n', dt);

% Compute annual maxima if needed
annualMaxima = tsEvaComputeAnnualMaxima(timeAndSeries);
```

---

### Pattern B: Parameter Selection Guidelines

```matlab
% Time window selection
dataSpan = max(timeAndSeries(:,1)) - min(timeAndSeries(:,1));
timeWindow = dataSpan / 5;  % Start with 1/5 of data span
% Typical: 5-10 years for short series, 30-50 years for climate projections

% Peak distance selection  
dt = tsEvaGetTimeStep(timeAndSeries);
minPeakDistanceInDays = max(3, 3*dt);  % At least 3 time steps
% For monthly data: 30-60 days
% For daily data: 2-5 days  
% For sub-daily: 1-3 days

% Threshold quantile
thresholdQuantile = 0.95;  % Gives ~18 events/year for daily data
% Lower (0.90): More events, less bias, more variance
% Higher (0.98): Fewer events, more bias, less variance
```

---

### Pattern C: Ensemble Analysis

**Use Case**: Quantify uncertainty by varying parameters

```matlab
% Define parameter ranges
timeWindows = [365.25*5, 365.25*6, 365.25*8];  % Multiple time windows
thresholdQuantiles = [0.93, 0.95, 0.97];  % Multiple thresholds

% Run ensemble
[ensembleParams] = tsEnsembleEvaParams(timeAndSeries, ...
    timeWindows, thresholdQuantiles);

% Analyze spread
meanReturnLevel = mean(ensembleParams.returnLevels, 3);
stdReturnLevel = std(ensembleParams.returnLevels, 0, 3);
```

---

### Pattern D: Model Comparison (Stationary vs Non-Stationary)

```matlab
% Fit both models
[statParams] = tsEvaStationary(timeAndSeries, 'minPeakDistanceInDays', 3);
[nonStatParams, statTransfData] = tsEvaNonStationary(timeAndSeries, timeWindow, ...
    'transfType', 'trend', 'minPeakDistanceInDays', 3);

% Compare return levels
returnPeriods = [10, 50, 100];
timeInstant = max(timeAndSeries(:,1));  % Present day

rlStat = tsEvaComputeReturnLevelsGEV(returnPeriods, statParams.gevParams);
rlNonStat = tsEvaComputeReturnLevelsGEVFromAnalysisObj(returnPeriods, ...
    timeInstant, nonStatParams, statTransfData);

fprintf('Return Level Comparison (Present):\n');
fprintf('RP\tStationary\tNon-Stationary\tDifference\n');
for i = 1:length(returnPeriods)
    fprintf('%d\t%.3f\t\t%.3f\t\t%.3f\n', ...
        returnPeriods(i), rlStat(i), rlNonStat(i), rlNonStat(i) - rlStat(i));
end
```

---

### Pattern E: Goodness-of-Fit Testing

```matlab
% For copulas
[gofResults] = tsCopulaGOFNonStat(copulaParams, marginalParams, peaks);

% Visual inspection
% - Check transformation to stationarity
tsEvaPlotTransfToStatFromAnalysisObj(nonStatEvaParams, statTransfData);

% - Check trend and variance
tsEvaPlotSeriesTrendStdDevFromAnalyisObj(nonStatEvaParams, statTransfData);
```

---

## Function Call Sequences

### Monovariate EVA
1. `tsEvaNonStationary` OR `tsEvaStationary`
2. `tsEvaComputeReturnLevels[GEV|GPD]FromAnalysisObj` OR direct GEV/GPD functions
3. `tsEvaPlot*FromAnalysisObj` OR direct plotting functions

### Copula EVA
1. `tsCopulaSampleJointPeaksMultiVariatePruning`
2. `tsCopulaExtremes` OR `tsCopulaYearExtrFit`
3. `tsCopulaComputeBivarRP` OR `tsCopulaYearExtrGetMltvrtRetPeriod`
4. `tsCopulaPlot*` OR `tsPlotBivarReturnPeriod`

### Data Transformation (Advanced)
1. `tsEvaTransformSeriesToStationaryTrend*` (if manual transformation needed)
2. `tsEvaStationary` on transformed data
3. Manual back-transformation using transformation parameters

---

## Best Practices

### DO
- Start with visual inspection of time series
- Test for stationarity before assuming non-stationarity
- Use appropriate time windows relative to data span
- Ensure peak independence through adequate separation
- Report confidence intervals and uncertainty
- Validate results against physical understanding

### DON'T  
- Use Frank copula for >2 variables (multivariate)
- Set time window longer than 1/3 of data span
- Ignore seasonal patterns in seasonal data
- Over-interpret small samples (< 30 years for long-term trends)
- Extrapolate beyond reasonable return periods (max ~3x data span)

---

## Troubleshooting

### Issue: "Not enough peaks"
**Cause**: Threshold too high or peak separation too long  
**Solution**: Lower `thresholdQuantile` or reduce `minPeakDistanceInDays`

### Issue: "Transformation doesn't look stationary"
**Cause**: Time window too short or data has multiple scales of variability  
**Solution**: Increase time window or use seasonal transformation

### Issue: "Copula fit fails"
**Cause**: Insufficient joint events or incompatible copula family  
**Solution**: Check joint peak count (need 50+), verify copula family for dimensionality

### Issue: "Return levels seem unrealistic"
**Cause**: Extrapolation too far, poor fit, or data quality issues  
**Solution**: Check data quality, assess goodness-of-fit, limit return period range

---

## Example Script References

For complete working examples, see:
- **Monovariate**: `examplesMonovariate/exampleGenerateSeriesEVAGraphs.m`
- **Seasonal**: `examplesMonovariate/exampleGenerateSeriesEVAGraphs_ciPercentile.m`
- **GPD**: `examplesMonovariate/exampleSPISeries.m`
- **Bivariate Copula**: `examplesCopula/exampleCopulaBivariateNonStationary.m`
- **Case Studies**: `examplesCopula/caseStudy01.m`, `caseStudy02.m`, `caseStudy03.m`

Always adapt these examples to your specific data and research questions.
