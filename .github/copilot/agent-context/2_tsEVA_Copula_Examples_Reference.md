# tsEVA 2.0 – Copula / Multivariate Examples Reference

This document describes the copula-based multivariate extreme value analysis workflows
implemented in the official tsEVA 2.0 example scripts:

- `caseStudy01.m` – bivariate compound flooding (GPD margins)
- `caseStudy02.m` – trivariate spatial dependence (GPD margins)
- `caseStudy03.m` – bivariate compound climate extremes (GEV margins)

All examples implement the Transformed-Stationary (TS) approach:

1. Transform non-stationary marginals to stationary
2. Apply stationary EVA (GPD or GEV)
3. Model dependence with a copula
4. Generate joint extremes via Monte Carlo simulation
5. Assess goodness of fit and compute return periods
6. Visualize joint distributions and diagnostics

---

## Quick Function Index

### Core Analysis Functions
- `tsCopulaExtremes()` - Main copula analysis (margins + dependence + event pairing)
- `tsCopulaMontecarlo()` - Monte Carlo simulation of joint extremes
- `tsCopulaGOFNonStat()` - Goodness-of-fit testing for copula model
- `tsCopulaComputeBivarRP()` - Compute bivariate return periods (OR/AND)
- `tsCopulaComputeandPlotBivarRP()` - Compute and plot bivariate return periods

### Plotting Functions
- `tsCopulaPlotBivariate()` - Comprehensive bivariate diagnostic plots (scatter, marginals, GOF, return period)
- `tsCopulaPlotTrivariate()` - Comprehensive trivariate diagnostic plots (scatter matrix, marginals, GOF)
- `tsCopulaPlotTrivariateWithMap()` - Trivariate plots with spatial map overlay
- `tsCopulaPlotJointReturnPeriod()` - Joint return period contour plot
- `tsCopulaPeakExtrPlotSctrBivar()` - Peak extremes scatter plot (bivariate)

### Advanced/Specialized Functions
- `tsCopulaFit()` - Low-level copula fitting (use `tsCopulaExtremes` instead in most cases)
- `tsCopulaRnd()` - Generate random samples from fitted copula
- `tsCopulaCdfFromSamples()` - Empirical CDF from samples
- `tsCopulaSampleJointPeaksMultiVariatePruning()` - Advanced event pairing with pruning

### Year-Extremes Functions (GEV margins)
- `tsCopulaYearExtrFit()` - Fit copula to annual extremes
- `tsCopulaYearExtrRnd()` - Generate samples from year-extremes copula
- `tsCopulaYearExtrDistribution()` - Year-extremes joint distribution
- `tsCopulaYearExtrGetMltvrtRetPeriod()` - Multivariate return period for annual extremes
- `tsCopulaYearExtrPlotSctrBivar()` - Bivariate scatter (annual extremes)
- `tsCopulaYearExtrPlotSctrTrivar()` - Trivariate scatter (annual extremes)
- `tsCopulaYearExtrPlotJdistTrivar()` - Joint distribution plot (annual extremes, trivariate)

### Helper Functions
- `tsCopulaGetFamilyId()` - Get numeric ID from copula family name
- `tsCopulaGetFamilyFromId()` - Get copula family name from numeric ID

---

## Common Copula Workflow

### 1. Copula analysis object

```matlab
copulaAnalysis = tsCopulaExtremes(time, dataMatrix, ...
    'minPeakDistanceInDaysMonovarSampling', minDeltaUnivarSampli, ...
    'maxPeakDistanceInDaysMultivarSampling', maxDeltaMultivarSampli, ...
    'copulaFamily', copulaFamily, ...
    'transfType', transfType, ...
    'timeWindow', timeWindow, ...
    'ciPercentile', ciPercentile, ...
    'potPercentiles', potPercentiles, ...
    ... % optional: samplingOrder, peakType, marginalDistributions, smoothInd
);
```

- `time` must be numeric (MATLAB datenum)
- each column of `dataMatrix` represents one variable
- non-stationarity is handled through the TS transformation

---

### 2. Monte Carlo simulation

```matlab
monteCarloAnalysis = tsCopulaMontecarlo(copulaAnalysis, ...
    'nResample', N, ...
    'timeIndex', 'middle', ...
    ... % optional: nonStationarity
);
```

---

### 3. Goodness-of-fit

```matlab
gofStatistics = tsCopulaGOFNonStat(copulaAnalysis, monteCarloAnalysis);
```

Optional smoothing:

```matlab
gofStatistics = tsCopulaGOFNonStat(copulaAnalysis, monteCarloAnalysis, ...
                                  'smoothInd', 10);
```

---

### 4. Return periods (bivariate only)

```matlab
retPerAnalysis = tsCopulaComputeBivarRP(copulaAnalysis, monteCarloAnalysis);
```

---

### 5. Plotting

Bivariate:

```matlab
tsCopulaPlotBivariate(copulaAnalysis, monteCarloAnalysis, ...
    'gofStatistics', gofStatistics, ...
    'retPerAnalysis', retPerAnalysis, ...
    'ylbl', {'Var1','Var2'});
```

Trivariate:

```matlab
tsCopulaPlotTrivariate(copulaAnalysis, monteCarloAnalysis, ...
    'gofStatistics', gofStatistics, ...
    'varLabels', {'Var1','Var2','Var3'});
```

---

## Case Study 01 – Compound Flooding (Bivariate, GPD)

File: `caseStudy01.m`
Variables: river discharge and significant wave height
Copula: Gumbel
Marginal distributions: GPD
Transformation: `trendlinear`

### Key parameters

```matlab
ciPercentile = [99,99];
potPercentiles = [{95},{99}];

timeWindowJointDist = 365.25*40;

minDeltaUnivarSampli = [30,30];
maxDeltaMultivarSampli = 45;

copulaFamily = 'gumbel';
samplingOrder = [2,1];
marginalDistributions = 'gpd';
```

### Monte Carlo

```matlab
monteCarloAnalysis = tsCopulaMontecarlo(copulaAnalysis, ...
                                       'nResample',1000, ...
                                       'timeIndex','middle');
```

### Full code:
```matlab
load caseStudy01_data

% find the overlapping part of both data sources; define a 3-hourly time
% frame
timeCommon=(max(timeSWH(1),timeRiverDisch(1))):3/24:(min(timeSWH(end),timeRiverDisch(end)));

indexGoodData=find(~isnan(riverineDischarge));
indexGoodDataW=find(~isnan(SWH));

%interpolate river and SWH data using common time frame
riverineDischarge_=interp1(timeRiverDisch(indexGoodData),riverineDischarge(indexGoodData),timeCommon);
timeSeriesRiver=[timeCommon',riverineDischarge_'];

SWH_=interp1(timeSWH(indexGoodDataW),SWH(indexGoodDataW),timeCommon);
timeSeriesSWH=[timeCommon',SWH_'];

%percentile levels of univariate series (used for transformation)
ciPercentile = [99,99];  

% peak-over-threshold levels used for sampling of univariate series; it
% should be a 1-d cell array where each cell can have different size
potPercentiles=[{95},{99}];   %95 99                           

%Non-stationary time window (in days) used for time-varying joint
%distribution
timeWindowJointDist = 365.25*40;   
%minimum distance (in days) between univariate peaks 
minDeltaUnivarSampli=[30,30]; %30, 30

 %maximum distance (in days) between multivariate peaks; can either take
 %one value or has to have a format and size matching
 %size(nchoosek([1:numvar],2),1) where numvar is 2 in bivariate case, 3 in
 %trivariate case, and so on
maxDeltaMultivarSampli=45; %45   

%copula family; Gumbel, gaussian and Frank are possible choices
copulaFamily='gumbel';  

%methodology to perform univariate transformation from non-stationary to
%stationary
transfType='trendlinear';

marginalDistributions='gpd';
samplingOrder=[2,1];
%

[copulaAnalysis] = tsCopulaExtremes(timeSeriesRiver(:,1), ...
    [timeSeriesRiver(:,2),timeSeriesSWH(:,2)], ...
    'minPeakDistanceInDaysMonovarSampling',minDeltaUnivarSampli, ...
    'maxPeakDistanceInDaysMultivarSampling',maxDeltaMultivarSampli, ...
    'copulaFamily',copulaFamily,...
    'transfType',transfType,'timeWindow',timeWindowJointDist,...
    'ciPercentile',ciPercentile,'potPercentiles',potPercentiles,...
    'marginalDistributions',marginalDistributions,'samplingOrder',samplingOrder);

[monteCarloAnalysis] = tsCopulaMontecarlo(copulaAnalysis,...
    'nResample',1000,'timeIndex','middle');

gofStatistics = tsCopulaGOFNonStat(copulaAnalysis, monteCarloAnalysis);

retPerAnalysis = tsCopulaComputeBivarRP(copulaAnalysis, monteCarloAnalysis);
    
axxArray = tsCopulaPlotBivariate(copulaAnalysis, monteCarloAnalysis, ...
    'gofStatistics', gofStatistics, ...
    'retPerAnalysis', retPerAnalysis, ...
    'ylbl', {'River discharge (m^3s^{-1})','SWH (m)'});
```

---

## Case Study 02 – Spatial Wave Extremes (Trivariate, GPD)

File: `caseStudy02.m`
Variables: significant wave height at three locations
Copula: Gumbel
Marginal distributions: GPD
Transformation: `trendlinear`

### Key parameters

```matlab
ciPercentile = [99,99,99];
potPercentiles = [{99},{99},{99}];

timeWindowNonStat = 365*40;

minDeltaUnivarSampli = [0.5,0.5,0.5];
maxDeltaMultivarSampli = 0.5;

copulaFamily = {'Gumbel'};
peakType = 'allExceedThreshold';
```

### Two-stage Monte Carlo

```matlab
mcStats = tsCopulaMontecarlo(copulaAnalysis,'nResample',10000);
mcPlot  = tsCopulaMontecarlo(copulaAnalysis,'nResample',300);
```

### Full code:
```matlab
%load data
load caseStudy02_data

%percentile levels of univariate series (used for transformation)

ciPercentile = [99,99,99];     

% peak-over-threshold levels used for sampling of univariate series; it
% should be a 1-d cell array where each cell can have different size
potPercentiles=[{99},{99},{99}];  

%Non-stationary time window (in days) used for time-varying joint
%distribution
timeWindowNonStat=365*40;
%minimum distance (in days) between univariate peaks 
minDeltaUnivarSampli=[0.5,0.5,0.5];

 %maximum distance (in days) between multivariate peaks; can either take
 %one value or has to have a format and size matching
 %size(nchoosek([1:numvar],2),1) where numvar is 2 in bivariate case, 3 in
 %trivariate case, and so on
maxDeltaMultivarSampli=0.5; 

%copula family; can be gaussian or gumbel
%copulaFamily = 'Gaussian';
copulaFamily = {'Gumbel'};  

%methodology to perform univariate transformation from non-stationary to
%stationary
transfType='trendlinear';
peakType='allExceedThreshold';
[copulaAnalysis] = tsCopulaExtremes(timeAndSeries1(:,1), ...
    [timeAndSeries1(:,2),timeAndSeries2(:,2),timeAndSeries3(:,2)], ...
    'minPeakDistanceInDaysMonovarSampling',minDeltaUnivarSampli, ...
    'maxPeakDistanceInDaysMultivarSampling',maxDeltaMultivarSampli, ...
    'copulaFamily',copulaFamily,...
    'transfType',transfType,'timeWindow',timeWindowNonStat,...
    'ciPercentile',ciPercentile,'potPercentiles',potPercentiles,...
    'peakType',peakType);

[monteCarloAnalysis1] = tsCopulaMontecarlo(copulaAnalysis,...
    'nResample',10000,'timeIndex','middle'); % large montecarlo good for statistics
[monteCarloAnalysis2] = tsCopulaMontecarlo(copulaAnalysis,...
    'nResample',300,'timeIndex','middle'); % smaller montecarlo good for plotting

[gofStatistics] = tsCopulaGOFNonStat(copulaAnalysis, monteCarloAnalysis1, 'smoothInd',10);

axxArray = tsCopulaPlotTrivariate(copulaAnalysis, monteCarloAnalysis2, ...
    'gofStatistics', gofStatistics, ...
    'varLabels', {'{Loc 1}_{SWH (m)}','{Loc 2}_{SWH (m)}','{Loc 3}_{SWH (m)}'});
```

---

## Case Study 03 – Temperature and Drought (Bivariate, GEV)

File: `caseStudy03.m`
Variables: surface temperature and SPEI
Copula: Gumbel
Marginal distributions: GEV
Transformation: `trendlinear`

### Key parameters

```matlab
marginalDistributions = 'gev';
peakType = 'allExceedThreshold';

ciPercentile = [99,99];
potPercentiles = [{75},{97}];

timeWindowNonStat = 365*35;

minDeltaUnivarSampli = [30,30];
maxDeltaMultivarSampli = 12*30;
```

### Monte Carlo with margins-only non-stationarity

```matlab
monteCarloAnalysis = tsCopulaMontecarlo(copulaAnalysis, ...
    'nResample',10000, ...
    'timeIndex','middle', ...
    'nonStationarity','margins');
```

### Full code:
```matlab
load caseStudy03_data

%percentile levels of univariate series (used for transformation)

ciPercentile = [99,99]; 

% peak-over-threshold levels used for sampling of univariate series; it
% should be a 1-d cell array where each cell can have different size
potPercentiles=[{75},{97}]; 

%Non-stationary time window (in days) used for time-varying joint
%distribution
timeWindowNonStat=365*35;

%minimum distance (in days) between univariate peaks 
minDeltaUnivarSampli=[30,30];  
%maximum distance (in days) between multivariate peaks; can either take
%one value or has to have a format and size matching
%size(nchoosek([1:numvar],2),1) where numvar is 2 in bivariate case, 3 in
%trivariate case, and so on
maxDeltaMultivarSampli=12*30; %12

%copula family; Gumbel, gaussian and Frank are possible choices
copulaFamily={'gumbel'};  

%methodology to perform univariate transformation from non-stationary to
%stationary
transfType='trendlinear';
peakType='allExceedThreshold';
marginalDistributions='gev';
[copulaAnalysis] = tsCopulaExtremes(timeAndSeries1(:,1), ...
    [timeAndSeries2(:,2),timeAndSeries1(:,2)], ...
    'minPeakDistanceInDaysMonovarSampling',minDeltaUnivarSampli, ...
    'maxPeakDistanceInDaysMultivarSampling',maxDeltaMultivarSampli, ...
    'copulaFamily',copulaFamily,...
    'transfType',transfType,'timeWindow',timeWindowNonStat,...
    'ciPercentile',ciPercentile,'potPercentiles',potPercentiles,...
    'peakType',peakType, ...
    'marginalDistributions',marginalDistributions,'smoothInd',10);

% large montecarlo good for statistics computation
[monteCarloAnalysis1] = tsCopulaMontecarlo(copulaAnalysis,...
    'nResample',10000,'timeIndex','middle','nonStationarity','margins');
% small montecarlo good for plotting
[monteCarloAnalysis2] = tsCopulaMontecarlo(copulaAnalysis,...
    'nResample',1000,'timeIndex','middle','nonStationarity','margins');

[gofStatistics] = tsCopulaGOFNonStat(copulaAnalysis, monteCarloAnalysis1, 'smoothInd',10);

[retPerAnalysis] = tsCopulaComputeBivarRP(copulaAnalysis, monteCarloAnalysis1);

axxArray = tsCopulaPlotBivariate(copulaAnalysis, monteCarloAnalysis2, ...
    'gofStatistics', gofStatistics, ...
    'retPerAnalysis', retPerAnalysis, ...
    'ylbl', {'- SPEI','Temp. °K'},'smoothInd',10,'rpPlot');
```


---

## Notes and Conventions

- `copulaFamily` may be passed as a string (e.g. `'gumbel'`)
  or as a cell array (e.g. `{'Gumbel'}`), as shown in the scripts.
- Two Monte Carlo runs (large for statistics, small for plotting) are recommended.
- GEV margins are fully supported in copula workflows when block-maxima logic is required.
- Non-stationarity can be applied to margins only via
  `'nonStationarity','margins'` in the Monte Carlo step.
