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

---

## Notes and Conventions

- `copulaFamily` may be passed as a string (e.g. `'gumbel'`)
  or as a cell array (e.g. `{'Gumbel'}`), as shown in the scripts.
- Two Monte Carlo runs (large for statistics, small for plotting) are recommended.
- GEV margins are fully supported in copula workflows when block-maxima logic is required.
- Non-stationarity can be applied to margins only via
  `'nonStationarity','margins'` in the Monte Carlo step.
