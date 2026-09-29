# tsEVA 2.0 Function Reference
This document lists ALL documented functions in tsEVA 2.0. **Only functions listed here may be referenced by the AI assistant.**
**Total Functions:** 98
---
## Table of Contents
1. [Core Analysis Functions](#core-analysis-functions)
   - [Monovariate EVA](#monovariate-eva) (25 functions)
   - [Copula Analysis](#copula-analysis) (23 functions)
2. [Plotting & Visualization](#plotting--visualization) (23 functions)
3. [Data Preparation & Utilities](#data-preparation--utilities) (27 functions)
4. [Function Details](#function-details)

---
## Quick Summary
| Category | Count | Description |
|----------|-------|-------------|
| Monovariate EVA | 25 | GEV, GPD, stationary, non-stationary analysis |
| Copula Analysis | 23 | Copula fitting, simulation, multivariate analysis |
| Plotting & Visualization | 23 | All plotting and graphical output functions |
| Data Preparation & Utilities | 27 | Data transformation, helper functions, utilities |

---
## Core Analysis Functions
### Monovariate EVA
Functions for univariate extreme value analysis using GEV and GPD distributions.

| Function | Description |
|----------|-------------|
| `Modified_MannKendall_test` | % FUNCTION INPUTS AND OUTPUTS |
| `tsApproxP` | No description available |
| `tsEstimateAverageSeasonality` | estimating the first 2 fourier components |
| `tsEvaComputeAnnualMaxima` | No description available |
| `tsEvaComputeAnnualMaximaMtx` | No description available |
| `tsEvaComputeMonthlyMaxima` | No description available |
| `tsEvaComputeReturnLevelsGEV` | tsEvaComputeReturnLevelsGEV: returns the return levels given the gev parameters  |
| `tsEvaComputeReturnLevelsGEVFromAnalysisObj` | No description available |
| `tsEvaComputeReturnLevelsGPD` | reference: Stuart Coles 2001, pag 81. sampleTimeHorizon and returnPeriods must b |
| `tsEvaComputeReturnLevelsGPDFromAnalysisObj` | percentile = nonStationaryEvaParams(2).parameters.percentile; dtSample = nonStat |
| `tsEvaDetrendTimeSeries` | No description available |
| `tsEvaNanRunningMean` | at both extremeties of the series, half windowSize is used which gradually incre |
| `tsEvaNanRunningPercentile` | tsEvaNanRunningPercentile: computes a runnig percentile for a given series, usin |
| `tsEvaNanRunningStatistics` | No description available |
| `tsEvaNanRunningVariance` | !!! series must be 0 averaged!! |
| `tsEvaNonStationary` | tsEvaNonStationary: performs the TS EVA analysis as described by Mentaschi et al |
| `tsEvaReduceOutputObjSize` | redStatTransData.stationarySeries = redStatTransData.stationarySeries(tsIndxs);  |
| `tsEvaRunningMeanTrend` | No description available |
| `tsEvaStationary` | tsEvaStationary: executes a regular stationary EVA on timeAndSeries. stationaryE |
| `tsEvaTransformSeriesToStatSeasonal_ciPercentile` | this function decomposes the series into a season-dependent trend and a season-d |
| `tsEvaTransformSeriesToStationaryMultiplicativeSeasonality` | this function decomposes the series into a season-dependent trend and a season-d |
| `tsEvaTransformSeriesToStationaryTrendLinear` | this function first calculates linear trend of the series and also return linear |
| `tsEvaTransformSeriesToStationaryTrendOnly` | further smoothing |
| `tsEvaTransformSeriesToStationaryTrendOnly_ciPercentile` | normalizing to standard deviation (just to be able to make acceptable graphs wit |
| `tsGpdNegShapeFit` | TSGPDNEGSHAPEFIT Fit GPD with negative shape constraint (k < 0) [paramEsts, para |

### Copula Analysis
Functions for multivariate dependence modeling using copulas.

| Function | Description |
|----------|-------------|
| `tsCopulaCdfFromSamples` | tsCopulaCdfFromSamples  Empirical copula CDF from sample points u       : (q x d |
| `tsCopulaComputeBivarRP` | tsCopulaComputeBivarRP computing of bivariate return period of type "AND" |
| `tsCopulaComputeandPlotBivarRP` | tsCopulaComputeandPlotBivarRP computing and plotting of bivariate return period  |
| `tsCopulaExtremes` | tsCopulaExtremes joint distribution of non-stationary compound extremes |
| `tsCopulaFit` | Replace negatives in the lower triangle |
| `tsCopulaGOFNonStat` | tsCopulaGOFNonStat estimation of copula goodness-of-fit and other battery of sta |
| `tsCopulaGetFamilyFromId` | No description available |
| `tsCopulaGetFamilyId` | No description available |
| `tsCopulaMontecarlo` | tsCopulaCompoundGPDMontecarlo pefrom Monte-Carlo simulation (resampling) from a  |
| `tsCopulaPeakExtrPlotSctrBivar` | No description available |
| `tsCopulaPlotBivariate` | tsCopulaPlotBivariate  plotting of joint peaks fitted by a copula |
| `tsCopulaPlotJointReturnPeriod` | plots the multivariate return period according to AND/OR Scenarios see https://d |
| `tsCopulaPlotTrivariate` | tsCopulaPlotTrivariate  plotting joint peaks and Monte-Carlo resampled values |
| `tsCopulaPlotTrivariateWithMap` | tsCopulaPlotTrivariate  plotting joint peaks and Monte-Carlo resampled values |
| `tsCopulaRnd` | generates the random vector from a copula distribution. if multivariate t or gau |
| `tsCopulaSampleJointPeaksMultiVariatePruning` | tsCopulaSampleJointPeaksMultiVariatePruning       multivariate peak-over-thresho |
| `tsCopulaYearExtrDistribution` | normal copula |
| `tsCopulaYearExtrFit` | this is a stationary set of return levels |
| `tsCopulaYearExtrGetMltvrtRetPeriod` | computes the multivariate return period according to Slavadori and De Michele 20 |
| `tsCopulaYearExtrPlotJdistTrivar` | No description available |
| `tsCopulaYearExtrPlotSctrBivar` | No description available |
| `tsCopulaYearExtrPlotSctrTrivar` | No description available |
| `tsCopulaYearExtrRnd` | No description available |

## Plotting & Visualization
Functions for creating plots and visualizations.

| Function | Description |
|----------|-------------|
| `plotGEV3D` | No description available |
| `plotGPD3D` | No description available |
| `plotGPD3DFromAnalysisObj` | No description available |
| `plotReturnLevelsGEVStationary` | nonStationaryEvaParams, stationaryTransformData are as these returned by functio |
| `plotReturnLevelsGPDStationary` | nonStationaryEvaParams, stationaryTransformData are as these returned by functio |
| `tsEvaPlotGEV3DFromAnalysisObj` | No description available |
| `tsEvaPlotGEVImageSc` | npdf = (year(maxTS) - year(minTS) + 1)*args.nPlottedTimesByYear; |
| `tsEvaPlotGEVImageScFromAnalysisObj` | No description available |
| `tsEvaPlotGPDImageSc` | npdf = (year(maxTS) - year(minTS) + 1)*args.nPlottedTimesByYear; |
| `tsEvaPlotGPDImageScFromAnalysisObj` | No description available |
| `tsEvaPlotReturnLevelsGEV` | No description available |
| `tsEvaPlotReturnLevelsGEVFromAnalysisObj` | nonStationaryEvaParams, stationaryTransformData are as these returned by functio |
| `tsEvaPlotReturnLevelsGPD` | No description available |
| `tsEvaPlotReturnLevelsGPDFromAnalysisObj` | nonStationaryEvaParams, stationaryTransformData are as these returned by functio |
| `tsEvaPlotSeasonalityGev` | tsEvaPlotSeasonalityGev plots a single year of data adding the series of monthly |
| `tsEvaPlotSeasonalityGevFromAnalysisObj` | No description available |
| `tsEvaPlotSeriesTrendStdDev` | No description available |
| `tsEvaPlotSeriesTrendStdDevFromAnalyisObj` | No description available |
| `tsEvaPlotTransfToStat` | stdThirdMom = third root of the third statistical momentum stdFouthMom = fourth  |
| `tsEvaPlotTransfToStatFromAnalysisObj` | std3mom = nthroot(stationaryTransformData.statSer3Mom, 3.); std4mom = nthroot(st |
| `tsLcSubplotManager` | No description available |
| `tsPlotBivarReturnPeriod` | tsCopulaPlotJointReturnPeriod plotting of multivariate return periods |
| `tsPlotSeriesPotGPDRetLevFromAnalysisObj` | computing the return levels |

## Data Preparation & Utilities
Helper functions for data manipulation and transformation.

| Function | Description |
|----------|-------------|
| `cvineOrder` | Root-first C-vine order from a Gumbel-θ matrix (alpha). Heuristic: pick node wit |
| `tsEVstatistics` | Evangelos Voukouvalas, Michalis Vousdoukas 2015 gevMaxima can be annual or month |
| `tsEasyParseNamedArgs` | No description available |
| `tsEmpirical` | tsEmpirical Empirical copula [ C ] = tsEmpirical( U ) retuns a variable C contai |
| `tsEnsemble` | calls tsEnsembleEvaParams and tsEnsembleStatTransfData. ASSUMING THAT ALL THE SE |
| `tsEnsembleEvaParams` | From a cell array of non nonStatEvaParams computes the average nonStatEvaParamsA |
| `tsEnsembleStatTransfData` | From a cell array of non stationaryTransformData computes the average stationary |
| `tsEstimateConfidenceIntervalOfRL` | % tsEstimateConfidenceIntervalOfRL: estimates the confidence interval of the ret |
| `tsEvaFillSeries` | ensuring monotonic time vector and no dubplicates. In case of dubplicate time st |
| `tsEvaGetReturnPeriodOfLevelGEV` | GEV |
| `tsEvaGetReturnPeriodOfLevelGPD` | No description available |
| `tsEvaGetTimeStep` | No description available |
| `tsEvaSampleData` | args.pcts = [50 70 85:2:99 99.1:0.1:99.9 99.91:0.01:99.99]; |
| `tsGet` | No description available |
| `tsGetNumberPerYear` | function nperYear=tsGetNumberPerYear(ms,locs) Gives number of events per year fr |
| `tsGetPOT` | function POTdata=tsGetPOT(ms,pcts,desiredEventsPerYear) Gets POT using an automa |
| `tsGetReturnPeriodOfLevel` | given a list of retPeriod and corresponding retLevel with error, estimates the r |
| `tsInterp1Extrap` | function nvals=interp1Extrap(Tr,values,RP,logExtrap) This function interpolates  |
| `tsLinearExtrapolation` | Linearly extrapolates time series x,y to points extrax, according to the slope o |
| `tsLoglogExtrapolation` | Extrapolates exponentiallytime series x,y to points extrax, according to the slo |
| `tsPseudoObservations` | No description available |
| `tsRankmax` | tsRankmax Returns vector of one-based ranks for each element |
| `tsRemoveConstantSubseries` | No description available |
| `tsRoundSDate` | function [sdround,dvec,sdunique,dvunique]=tsRoundSDate(sd,sd_precision) ROunds s |
| `tsSameValuesSegmentation` | function [inds,rinds]=tsSameValuesSegmentation(iii,val) separates segments of sa |
| `tsTimeSeriesToPointData` | tsTimeSeriesToPointData: given a ms produces a structure pointData like the one  |
| `tsYear` | No description available |

---
## Function Details
Comprehensive documentation for each function.

### Monovariate EVA
#### Modified_MannKendall_test
- **Signature**: `[tau, z, p, H] = Modified_MannKendall_test(t, X, alpha, alpha_ac)`
- **Description**: % FUNCTION INPUTS AND OUTPUTS
- **Category**: Monovariate EVA
- **Inputs**:
  - `t`
  - `X`
  - `alpha`
  - `alpha_ac`
- **Outputs**:
  - `tau`
  - `z`
  - `p`
  - `H`

#### tsApproxP
- **Signature**: `[Pval]=tsApproxP(N,copulaFamily,rho,nu,snSample,s2Sample)`
- **Category**: Monovariate EVA
- **Inputs**:
  - `N`
  - `copulaFamily`
  - `rho`
  - `nu`
  - `snSample`
  - `s2Sample`
- **Outputs**:
  - `Pval`

#### tsEstimateAverageSeasonality
- **Signature**: `averageSeasonalitySeries = tsEstimateAverageSeasonality( timeStamps, seasonalitySeries )`
- **Description**: estimating the first 2 fourier components
- **Category**: Monovariate EVA
- **Inputs**:
  - `timeStamps`
  - `seasonalitySeries`
- **Outputs**:
  - `averageSeasonalitySeries`

#### tsEvaComputeAnnualMaxima
- **Signature**: `[annualMax, annualMaxDate, annualMaxIndx] = tsEvaComputeAnnualMaxima(timeAndSeries)`
- **Category**: Monovariate EVA
- **Inputs**:
  - `timeAndSeries`
- **Outputs**:
  - `annualMax`
  - `annualMaxDate`
  - `annualMaxIndx`

#### tsEvaComputeAnnualMaximaMtx
- **Signature**: `[annualMax, annualMaxDate, annualMaxIndx] = tsEvaComputeAnnualMaximaMtx(timeStamps, srs)`
- **Category**: Monovariate EVA
- **Inputs**:
  - `timeStamps`
  - `srs`
- **Outputs**:
  - `annualMax`
  - `annualMaxDate`
  - `annualMaxIndx`

#### tsEvaComputeMonthlyMaxima
- **Signature**: `[monthlyMax, monthlyMaxDate, monthlyMaxIndx] = tsEvaComputeMonthlyMaxima(timeAndSeries)`
- **Category**: Monovariate EVA
- **Inputs**:
  - `timeAndSeries`
- **Outputs**:
  - `monthlyMax`
  - `monthlyMaxDate`
  - `monthlyMaxIndx`

#### tsEvaComputeReturnLevelsGEV
- **Signature**: `[returnLevels, returnLevelsErr] = tsEvaComputeReturnLevelsGEV( epsilon, sigma, mu, epsilonStdErr, sigmaStdErr, muStdErr, returnPeriodsInDts, varargin )`
- **Description**: tsEvaComputeReturnLevelsGEV: returns the return levels given the gev parameters and their standard error. The parameter returnPeriodsInDts contains the return period expressed in a time unit that corresponds to the size of the time segments where we
- **Category**: Monovariate EVA
- **Inputs**:
  - `epsilon`
  - `sigma`
  - `mu`
  - `epsilonStdErr`
  - `sigmaStdErr`
  - `muStdErr`
  - `returnPeriodsInDts`
  - `varargin`
- **Outputs**:
  - `returnLevels`
  - `returnLevelsErr`

#### tsEvaComputeReturnLevelsGEVFromAnalysisObj
- **Signature**: `[returnLevels, returnLevelsErr, returnLevelsErrFit, returnLevelsErrTransf] = tsEvaComputeReturnLevelsGEVFromAnalysisObj(nonStationaryEvaParams, returnPeriodsInYears, varargin)`
- **Category**: Monovariate EVA
- **Inputs**:
  - `nonStationaryEvaParams`
  - `returnPeriodsInYears`
  - `varargin`
- **Outputs**:
  - `returnLevels`
  - `returnLevelsErr`
  - `returnLevelsErrFit`
  - `returnLevelsErrTransf`

#### tsEvaComputeReturnLevelsGPD
- **Signature**: `[returnLevels, returnLevelsErr] = tsEvaComputeReturnLevelsGPD( epsilon, sigma, threshold, epsilonStdErr, sigmaStdErr, thresholdStdErr, nPeaks, sampleTimeHorizon, returnPeriods )`
- **Description**: reference: Stuart Coles 2001, pag 81. sampleTimeHorizon and returnPeriods must be in the same units, e.g. years
- **Category**: Monovariate EVA
- **Inputs**:
  - `epsilon`
  - `sigma`
  - `threshold`
  - `epsilonStdErr`
  - `sigmaStdErr`
  - `thresholdStdErr`
  - `nPeaks`
  - `sampleTimeHorizon`
  - `returnPeriods`
- **Outputs**:
  - `returnLevels`
  - `returnLevelsErr`

#### tsEvaComputeReturnLevelsGPDFromAnalysisObj
- **Signature**: `[returnLevels, returnLevelsErr, returnLevelsErrFit, returnLevelsErrTransf] = tsEvaComputeReturnLevelsGPDFromAnalysisObj(nonStationaryEvaParams, returnPeriodsInYears, varargin )`
- **Description**: percentile = nonStationaryEvaParams(2).parameters.percentile; dtSample = nonStationaryEvaParams(2).parameters.timeDeltaYears;
- **Category**: Monovariate EVA
- **Inputs**:
  - `nonStationaryEvaParams`
  - `returnPeriodsInYears`
  - `varargin`
- **Outputs**:
  - `returnLevels`
  - `returnLevelsErr`
  - `returnLevelsErrFit`
  - `returnLevelsErrTransf`

#### tsEvaDetrendTimeSeries
- **Signature**: `[ detrendSeries, trendSeries, filledTimeStamps, filledSeries, nRunMn ] = tsEvaDetrendTimeSeries( timeStamps, series, timeWindow, varargin )`
- **Category**: Monovariate EVA
- **Inputs**:
  - `timeStamps`
  - `series`
  - `timeWindow`
  - `varargin`
- **Outputs**:
  - `detrendSeries`
  - `trendSeries`
  - `filledTimeStamps`
  - `filledSeries`
  - `nRunMn`

#### tsEvaNanRunningMean
- **Signature**: `[ rnmn ] = tsEvaNanRunningMean( series, windowSize )`
- **Description**: at both extremeties of the series, half windowSize is used which gradually increases to reach windowSize; once windowSize is reached, windowSize is rolled throghout the series
- **Category**: Monovariate EVA
- **Inputs**:
  - `series`
  - `windowSize`
- **Outputs**:
  - `rnmn`

#### tsEvaNanRunningPercentile
- **Signature**: `[ rnprcnt, stdError ] = tsEvaNanRunningPercentile( series, windowSize, percent, varargin )`
- **Description**: tsEvaNanRunningPercentile: computes a runnig percentile for a given series, using a window with a size given by windowSize.
- **Category**: Monovariate EVA
- **Inputs**:
  - `windowSize`: size of the window for the running percentile. Cannot be < 1000
  - `percent`: percent level to which the percentile is compute.
  - `percentDelta`: delta for the computation of a percentile interval
  - `parameter`: around the requested percentage. If for example
  - `parameter`: percent==90 and percentDelta==1, then the 89th, 90th and
  - `parameter`: 91st percentiles are computed. Default value: 1 if
  - `parameter`: windowSize > 2000, 2 if 2000 > windowsize > 1000.
  - `nLowLimit`: minimum number of non nan elements for a window for
  - `parameter`: percentile computation.
  - `rnprcnt`: approximated running percentile.
  - `parameter`: How it works:
  - `parameter`: let's suppose that percent == 90.
  - `parameter`: For the first window we compute the right percentile using matlab
  - `parameter`: function prctile, for percentages 89, 90, 91.
  - `parameter`: Then for each step, we update these percentages on the basis
  - `parameter`: of the quitting values and incoming values,
  - `parameter`: and interpolate an approximated percentile for the requested percentage.
- **Outputs**:
  - `rnprcnt`
  - `stdError`

#### tsEvaNanRunningStatistics
- **Signature**: `[ rnmn, rnvar, rn3mom, rn4mom ] = tsEvaNanRunningStatistics( series, windowSize )`
- **Category**: Monovariate EVA
- **Inputs**:
  - `series`
  - `windowSize`
- **Outputs**:
  - `output`: to the forth.
  - `rnmn`: running mean
  - `rnvar`: running variance
  - `rn3mom`: running third statistical momentum
  - `rn4mom`: running fourth statistical momentum

#### tsEvaNanRunningVariance
- **Signature**: `[ rnmn ] = tsEvaNanRunningVariance( series, windowSize )`
- **Description**: !!! series must be 0 averaged!!
- **Category**: Monovariate EVA
- **Inputs**:
  - `series`
  - `windowSize`
- **Outputs**:
  - `rnmn`

#### tsEvaNonStationary
- **Signature**: `[nonStationaryEvaParams, stationaryTransformData, isValid] = tsEvaNonStationary( timeAndSeries, timeWindow, varargin )`
- **Description**: tsEvaNonStationary: performs the TS EVA analysis as described by Mentaschi et al 2016.
- **Category**: Monovariate EVA
- **Inputs**:
  - `timeAndSeries`: array with shape nx2, with the time stamps in the first
  - `parameter`: column and the values in the second.
  - `timeWindow`: time window for the transformation expressed in days.
  - `transfType`: can assume values
  - `parameter`: 1) 'trend': long term variability. The trend is computed
  - `parameter`: with a running mean, the ci with the running standard deviation.
  - `parameter`: 2) 'seasonal': long term + seasonal variability. The trend is computed
  - `parameter`: with a running mean, the ci with the running
  - `parameter`: standard deviation.
  - `parameter`: 3) 'trendCIPercentile': long term variability. The trend is computed
  - `parameter`: with a running mean, the ci with the running xx percentile.
  - `parameter`: Using this option the argument ciPercentile is
  - `parameter`: mandatory.
  - `parameter`: 4)  'trendlinear':long term variability. the trend is
  - `parameter`: computed with a linear fit, the ci with a linear fit
  - `parameter`: of percentile set by the user
  - `parameter`: % sample calls
  - `parameter`: nonStatEvaParams = tsEvaNonStationary(ms, timeWindow, 'potPercentiles',[95], 'minPeakDistanceInDays', 3)
  - `parameter`: samples POT data using a fixed 95 percentile threshold, with peaks at
  - `parameter`: a minimum distance of 3 days, looking for a threshold so that we have an average
  - `parameter`: of 5 events every year.
  - `parameter`: nonStatEvaParams = tsEvaNonStationary(ms, timeWindow, 'minPeakDistanceInDays', 3, 'desiredeventsperyear', 6)
  - `parameter`: samples POT data looking for a threshold so that we have an average
  - `parameter`: of 6 events every year.
  - `parameter`: nonStatEvaParams = tsEvaNonStationary(ms, timeWindow, 'minPeakDistanceInDays', 3, 'trasftype', 'trendCIPercentile', 'ciPercentile', 99)
  - `parameter`: for the transformation uses instead of the moving standard deviation,
  - `parameter`: the moving 99th percentile.
  - `parameter`: % %%%%%%%%%%%%%
- **Outputs**:
  - `nonStationaryEvaParams`
  - `stationaryTransformData`
  - `isValid`

#### tsEvaReduceOutputObjSize
- **Signature**: `[ redNonStatEvaParams, redStatTransData ] = tsEvaReduceOutputObjSize( nonStationaryEvaParams, stationaryTransformData, newTimeStamps, varargin )`
- **Description**: redStatTransData.stationarySeries = redStatTransData.stationarySeries(tsIndxs); redStatTransData.nonStatSeries = redStatTransData.nonStatSeries(tsIndxs);
- **Category**: Monovariate EVA
- **Inputs**:
  - `nonStationaryEvaParams`
  - `stationaryTransformData`
  - `newTimeStamps`
  - `varargin`
- **Outputs**:
  - `redNonStatEvaParams`
  - `redStatTransData`

#### tsEvaRunningMeanTrend
- **Signature**: `[ trendSeries, filledTimeStamps, filledSeries, nRunMn ] = tsEvaRunningMeanTrend( timeStamps, series, timeWindow)`
- **Category**: Monovariate EVA
- **Inputs**:
  - `timeStamps`
  - `series`
  - `timeWindow`
- **Outputs**:
  - `trendSeries`
  - `filledTimeStamps`
  - `filledSeries`
  - `nRunMn`

#### tsEvaStationary
- **Signature**: `[stationaryEvaParams, isValid] = tsEvaStationary( timeAndSeries, varargin )`
- **Description**: tsEvaStationary: executes a regular stationary EVA on timeAndSeries. stationaryEvaParams includes the parameters estimated for GEV and GPD % %%%%%%%%%%%%%
- **Category**: Monovariate EVA
- **Inputs**:
  - `timeAndSeries`
  - `varargin`
- **Outputs**:
  - `stationaryEvaParams`
  - `isValid`

#### tsEvaTransformSeriesToStatSeasonal_ciPercentile
- **Signature**: `[trasfData] = tsEvaTransformSeriesToStatSeasonal_ciPercentile( timeStamps, series, timeWindow, percentile, varargin  )`
- **Description**: this function decomposes the series into a season-dependent trend and a season-dependent standard deviation. The season-dependent standard deviation is given by a seasonal factor multiplied by a slowly varying standard deviation.
- **Category**: Monovariate EVA
- **Inputs**:
  - `timeStamps`
  - `series`
  - `timeWindow`
  - `percentile`
  - `varargin`
- **Outputs**:
  - `trasfData`

#### tsEvaTransformSeriesToStationaryMultiplicativeSeasonality
- **Signature**: `[trasfData] = tsEvaTransformSeriesToStationaryMultiplicativeSeasonality( timeStamps, series, timeWindow, varargin )`
- **Description**: this function decomposes the series into a season-dependent trend and a season-dependent standard deviation. The season-dependent standard deviation is given by a seasonal factor multiplied by a slowly varying standard deviation.
- **Category**: Monovariate EVA
- **Inputs**:
  - `timeStamps`
  - `series`
  - `timeWindow`
  - `varargin`
- **Outputs**:
  - `trasfData`

#### tsEvaTransformSeriesToStationaryTrendLinear
- **Signature**: `[trasfData] = tsEvaTransformSeriesToStationaryTrendLinear( timeStamps, series, timeWindow, percentile, varargin )`
- **Description**: this function first calculates linear trend of the series and also return linear trend of percetile series (at any particular percentile level set by the user) through a fitting algorithm
- **Category**: Monovariate EVA
- **Inputs**:
  - `timeStamps`
  - `series`
  - `timeWindow`
  - `percentile`
  - `varargin`
- **Outputs**:
  - `trasfData`

#### tsEvaTransformSeriesToStationaryTrendOnly
- **Signature**: `[trasfData] = tsEvaTransformSeriesToStationaryTrendOnly( timeStamps, series, timeWindow, varargin )`
- **Description**: further smoothing
- **Category**: Monovariate EVA
- **Inputs**:
  - `timeStamps`
  - `series`
  - `timeWindow`
  - `varargin`
- **Outputs**:
  - `trasfData`

#### tsEvaTransformSeriesToStationaryTrendOnly_ciPercentile
- **Signature**: `[trasfData] = tsEvaTransformSeriesToStationaryTrendOnly_ciPercentile( timeStamps, series, timeWindow, percentile, varargin )`
- **Description**: normalizing to standard deviation (just to be able to make acceptable graphs with the scripts of this library)
- **Category**: Monovariate EVA
- **Inputs**:
  - `timeStamps`
  - `series`
  - `timeWindow`
  - `percentile`
  - `varargin`
- **Outputs**:
  - `trasfData`

#### tsGpdNegShapeFit
- **Signature**: `[paramEsts, paramCIs] = tsGpdNegShapeFit(data, alphaCI)`
- **Description**: TSGPDNEGSHAPEFIT Fit GPD with negative shape constraint (k < 0) [paramEsts, paramCIs] = tsGpdNegShapeFit(data, alphaCI) mimics gpfit but constrains the shape parameter to be negative.
- **Category**: Monovariate EVA
- **Inputs**:
  - `data`
  - `alphaCI`
- **Outputs**:
  - `paramEsts`
  - `paramCIs`

### Copula Analysis
#### tsCopulaCdfFromSamples
- **Signature**: `[C, se] = tsCopulaCdfFromSamples(u, Usample)`
- **Description**: tsCopulaCdfFromSamples  Empirical copula CDF from sample points u       : (q x d) query points in [0,1]^d Usample : (M x d) sample points from the fitted copula C       : (q x 1) empirical CDF estimates se      : (q x 1) standard errors sqrt(C(1-C)/M) [optional]
- **Category**: Copula Analysis
- **Inputs**:
  - `u`
  - `Usample`
- **Outputs**:
  - `C`
  - `se`

#### tsCopulaComputeBivarRP
- **Signature**: ` [rpAnalysis] = tsCopulaComputeBivarRP(copulaAnalysis, monteCarloAnalysis, varargin)`
- **Description**: tsCopulaComputeBivarRP computing of bivariate return period of type "AND"
- **Category**: Copula Analysis
- **Inputs**:
  - `copulaAnalysis`
  - `monteCarloAnalysis`
  - `varargin`
- **Outputs**:
  - `rpAnalysis`

#### tsCopulaComputeandPlotBivarRP
- **Signature**: ` [rpAnalysis,hFig]=tsCopulaComputeandPlotBivarRP(copulaAnalysis,varargin)`
- **Description**: tsCopulaComputeandPlotBivarRP computing and plotting of bivariate return period of type "AND"
- **Category**: Copula Analysis
- **Inputs**:
  - `copulaAnalysis`
  - `varargin`
- **Outputs**:
  - `rpAnalysis`
  - `hFig`

#### tsCopulaExtremes
- **Signature**: `[CopulaAnalysis] = tsCopulaExtremes(inputtimestamps,inputtimeseries, varargin)`
- **Description**: tsCopulaExtremes joint distribution of non-stationary compound extremes
- **Category**: Copula Analysis
- **Inputs**:
  - `inputtimestamps`
  - `inputtimeseries`
  - `varargin`
- **Outputs**:
  - `CopulaAnalysis`

#### tsCopulaFit
- **Signature**: `copulaParam = tsCopulaFit(copulaFamily, uProb)`
- **Description**: Replace negatives in the lower triangle
- **Category**: Copula Analysis
- **Inputs**:
  - `copulaFamily`
  - `uProb`
- **Outputs**:
  - `copulaParam`

#### tsCopulaGOFNonStat
- **Signature**: `[gofStatistics] = tsCopulaGOFNonStat(copulaAnalysis, monteCarloAnalysis, varargin)`
- **Description**: tsCopulaGOFNonStat estimation of copula goodness-of-fit and other battery of statistics
- **Category**: Copula Analysis
- **Inputs**:
  - `copulaAnalysis`
  - `monteCarloAnalysis`
  - `varargin`
- **Outputs**:
  - `gofStatistics`

#### tsCopulaGetFamilyFromId
- **Signature**: `copulaFamily = tsCopulaGetFamilyFromId(familyId)`
- **Category**: Copula Analysis
- **Inputs**:
  - `familyId`
- **Outputs**:
  - `copulaFamily`

#### tsCopulaGetFamilyId
- **Signature**: `familyId = tsCopulaGetFamilyId(copulaFamily)`
- **Category**: Copula Analysis
- **Inputs**:
  - `copulaFamily`
- **Outputs**:
  - `familyId`

#### tsCopulaMontecarlo
- **Signature**: `[monteCarloAnalysis] = tsCopulaMontecarlo(copulaAnalysis, varargin)`
- **Description**: tsCopulaCompoundGPDMontecarlo pefrom Monte-Carlo simulation (resampling) from a pre-determined copula function
- **Category**: Copula Analysis
- **Inputs**:
  - `copulaAnalysis`
  - `varargin`
- **Outputs**:
  - `monteCarloAnalysis`

#### tsCopulaPeakExtrPlotSctrBivar
- **Signature**: `[handles] = tsCopulaPeakExtrPlotSctrBivar(monteCarloRsmpl, yMaxLevel, varargin)`
- **Category**: Copula Analysis
- **Inputs**:
  - `monteCarloRsmpl`
  - `yMaxLevel`
  - `varargin`
- **Outputs**:
  - `handles`

#### tsCopulaPlotBivariate
- **Signature**: `[axxArray] = tsCopulaPlotBivariate(copulaAnalysis, monteCarloAnalysis, varargin)`
- **Description**: tsCopulaPlotBivariate  plotting of joint peaks fitted by a copula
- **Category**: Copula Analysis
- **Inputs**:
  - `copulaAnalysis`
  - `monteCarloAnalysis`
  - `varargin`
- **Outputs**:
  - `axxArray`

#### tsCopulaPlotJointReturnPeriod
- **Signature**: ` tsCopulaPlotJointReturnPeriod(copulaAnalysis,varargin)`
- **Description**: plots the multivariate return period according to AND/OR Scenarios see https://doi.org/10.1002/2015WR017225 Parts of the code were reworked from MvCAT toolbox https://doi.org/10.1002/2016WR020242 Bahmanpour, M.H., 2023
- **Category**: Copula Analysis
- **Inputs**:
  - `copulaAnalysis`
  - `varargin`

#### tsCopulaPlotTrivariate
- **Signature**: `[axxArray] = tsCopulaPlotTrivariate(copulaAnalysis, monteCarloAnalysis, varargin)`
- **Description**: tsCopulaPlotTrivariate  plotting joint peaks and Monte-Carlo resampled values
- **Category**: Copula Analysis
- **Inputs**:
  - `copulaAnalysis`
  - `monteCarloAnalysis`
  - `varargin`
- **Outputs**:
  - `axxArray`

#### tsCopulaPlotTrivariateWithMap
- **Signature**: `[axxArray] = tsCopulaPlotTrivariateWithMap(copulaAnalysis,gofStatistics,varargin)`
- **Description**: tsCopulaPlotTrivariate  plotting joint peaks and Monte-Carlo resampled values
- **Category**: Copula Analysis
- **Inputs**:
  - `copulaAnalysis`
  - `gofStatistics`
  - `varargin`
- **Outputs**:
  - `axxArray`

#### tsCopulaRnd
- **Signature**: `u = tsCopulaRnd(family, copulaPar, N, uProb)`
- **Description**: generates the random vector from a copula distribution. if multivariate t or gaussian or bivariate archimedean, it uses copularnd else (this would be the case of multivariate archimedean) it implements a simple C-Vine copula
- **Category**: Copula Analysis
- **Inputs**:
  - `family`: copula family
  - `copulaPar`: rho for gaussian, alpha for archimedean
  - `N`: sample size
  - `uProb`: pseudo-observation data (that is, transformed in uniform probability space)
  - `parameter`: from which the copula was built. Necessary in case the
  - `c`: vine copula must be constructed
- **Outputs**:
  - `u`

#### tsCopulaSampleJointPeaksMultiVariatePruning
- **Signature**: `[samplingAnalysis] = tsCopulaSampleJointPeaksMultiVariatePruning(inputtimestamps,inputtimeseries,varargin)`
- **Description**: tsCopulaSampleJointPeaksMultiVariatePruning       multivariate peak-over-threshold sampling of compound events
- **Category**: Copula Analysis
- **Inputs**:
  - `inputtimestamps`
  - `inputtimeseries`
  - `varargin`
- **Outputs**:
  - `samplingAnalysis`

#### tsCopulaYearExtrDistribution
- **Signature**: `[jpdf, jcdf] = tsCopulaYearExtrDistribution(retPeriod, copulaParam, varargin)`
- **Description**: normal copula
- **Category**: Copula Analysis
- **Inputs**:
  - `retPeriod`
  - `copulaParam`
  - `varargin`
- **Outputs**:
  - `jpdf`
  - `jcdf`

#### tsCopulaYearExtrFit
- **Signature**: `[retLev, copulaParam, yRetPer, yProb] = tsCopulaYearExtrFit(retPeriod, retLev, yMax, varargin)`
- **Description**: this is a stationary set of return levels
- **Category**: Copula Analysis
- **Inputs**:
  - `retPeriod`
  - `retLev`
  - `yMax`
  - `varargin`
- **Outputs**:
  - `retLev`
  - `copulaParam`
  - `yRetPer`
  - `yProb`

#### tsCopulaYearExtrGetMltvrtRetPeriod
- **Signature**: `[returnPeriod, prob] = tsCopulaYearExtrGetMltvrtRetPeriod(randomSample, level)`
- **Description**: computes the multivariate return period according to Slavadori and De Michele 2004, Salvadori et al. 2011, used by Zscheischler et al. 2017
- **Category**: Copula Analysis
- **Inputs**:
  - `randomSample`
  - `level`
- **Outputs**:
  - `returnPeriod`
  - `prob`

#### tsCopulaYearExtrPlotJdistTrivar
- **Signature**: `handles = tsCopulaYearExtrPlotJdistTrivar( retLev, jdist, varargin )`
- **Category**: Copula Analysis
- **Inputs**:
  - `retLev`
  - `jdist`
  - `varargin`
- **Outputs**:
  - `handles`

#### tsCopulaYearExtrPlotSctrBivar
- **Signature**: `[handles] = tsCopulaYearExtrPlotSctrBivar(monteCarloRsmpl, yMaxLevel, varargin)`
- **Category**: Copula Analysis
- **Inputs**:
  - `monteCarloRsmpl`
  - `yMaxLevel`
  - `varargin`
- **Outputs**:
  - `handles`

#### tsCopulaYearExtrPlotSctrTrivar
- **Signature**: `handles = tsCopulaYearExtrPlotSctrTrivar(monteCarloRsmpl, yMaxLevel, varargin)`
- **Category**: Copula Analysis
- **Inputs**:
  - `monteCarloRsmpl`
  - `yMaxLevel`
  - `varargin`
- **Outputs**:
  - `handles`

#### tsCopulaYearExtrRnd
- **Signature**: `[monteCarloRsmpl, resampleProb, resampleRetPer] = tsCopulaYearExtrRnd(retPeriod, retLev, copulaParam, nResample, varargin)`
- **Category**: Copula Analysis
- **Inputs**:
  - `retPeriod`
  - `retLev`
  - `copulaParam`
  - `nResample`
  - `varargin`
- **Outputs**:
  - `monteCarloRsmpl`
  - `resampleProb`
  - `resampleRetPer`

### Plotting & Visualization
#### plotGEV3D
- **Signature**: `phandles = plotGEV3D( X, timeStamps, epsilon, sigma, mu, varargin )`
- **Category**: Plotting & Visualization
- **Inputs**:
  - `X`
  - `timeStamps`
  - `epsilon`
  - `sigma`
  - `mu`
  - `varargin`
- **Outputs**:
  - `phandles`

#### plotGPD3D
- **Signature**: `phandles = plotGPD3D( X, timeStamps, epsilon, sigma, threshold, varargin )`
- **Category**: Plotting & Visualization
- **Inputs**:
  - `X`
  - `timeStamps`
  - `epsilon`
  - `sigma`
  - `threshold`
  - `varargin`
- **Outputs**:
  - `phandles`

#### plotGPD3DFromAnalysisObj
- **Signature**: `phandles = plotGPD3DFromAnalysisObj( X, nonStationaryEvaParams, stationaryTransformData )`
- **Category**: Plotting & Visualization
- **Inputs**:
  - `X`
  - `nonStationaryEvaParams`
  - `stationaryTransformData`
- **Outputs**:
  - `phandles`

#### plotReturnLevelsGEVStationary
- **Signature**: `phandles = plotReturnLevelsGEVStationary( nonStationaryEvaParams, varargin )`
- **Description**: nonStationaryEvaParams, stationaryTransformData are as these returned by function nonStationaryEvaJRCApproach timeIndex is the index at which the time varying analysis should be estimated.
- **Category**: Plotting & Visualization
- **Inputs**:
  - `nonStationaryEvaParams`
  - `varargin`
- **Outputs**:
  - `phandles`

#### plotReturnLevelsGPDStationary
- **Signature**: `phandles = plotReturnLevelsGPDStationary( nonStationaryEvaParams, varargin )`
- **Description**: nonStationaryEvaParams, stationaryTransformData are as these returned by function nonStationaryEvaJRCApproach timeIndex is the index at which the time varying analysis should be estimated.
- **Category**: Plotting & Visualization
- **Inputs**:
  - `nonStationaryEvaParams`
  - `varargin`
- **Outputs**:
  - `phandles`

#### tsEvaPlotGEV3DFromAnalysisObj
- **Signature**: `phandles = tsEvaPlotGEV3DFromAnalysisObj( X, nonStationaryEvaParams, stationaryTransformData, varargin )`
- **Category**: Plotting & Visualization
- **Inputs**:
  - `X`
  - `nonStationaryEvaParams`
  - `stationaryTransformData`
  - `varargin`
- **Outputs**:
  - `phandles`

#### tsEvaPlotGEVImageSc
- **Signature**: `phandles = tsEvaPlotGEVImageSc( Y, timeStamps, epsilon, sigma, mu, varargin )`
- **Description**: npdf = (year(maxTS) - year(minTS) + 1)*args.nPlottedTimesByYear;
- **Category**: Plotting & Visualization
- **Inputs**:
  - `Y`
  - `timeStamps`
  - `epsilon`
  - `sigma`
  - `mu`
  - `varargin`
- **Outputs**:
  - `phandles`

#### tsEvaPlotGEVImageScFromAnalysisObj
- **Signature**: `phandles = tsEvaPlotGEVImageScFromAnalysisObj( X, nonStationaryEvaParams, stationaryTransformData, varargin )`
- **Category**: Plotting & Visualization
- **Inputs**:
  - `X`
  - `nonStationaryEvaParams`
  - `stationaryTransformData`
  - `varargin`
- **Outputs**:
  - `phandles`

#### tsEvaPlotGPDImageSc
- **Signature**: `phandles = tsEvaPlotGPDImageSc( Y, timeStamps, epsilon, sigma, threshold, varargin )`
- **Description**: npdf = (year(maxTS) - year(minTS) + 1)*args.nPlottedTimesByYear;
- **Category**: Plotting & Visualization
- **Inputs**:
  - `Y`
  - `timeStamps`
  - `epsilon`
  - `sigma`
  - `threshold`
  - `varargin`
- **Outputs**:
  - `phandles`

#### tsEvaPlotGPDImageScFromAnalysisObj
- **Signature**: `phandles = tsEvaPlotGPDImageScFromAnalysisObj( Y, nonStationaryEvaParams, stationaryTransformData, varargin )`
- **Category**: Plotting & Visualization
- **Inputs**:
  - `Y`
  - `nonStationaryEvaParams`
  - `stationaryTransformData`
  - `varargin`
- **Outputs**:
  - `phandles`

#### tsEvaPlotReturnLevelsGEV
- **Signature**: `phandles = tsEvaPlotReturnLevelsGEV( epsilon, sigma, mu, epsilonStdErr, sigmaStdErr, muStdErr, varargin  )`
- **Category**: Plotting & Visualization
- **Inputs**:
  - `epsilon`
  - `sigma`
  - `mu`
  - `epsilonStdErr`
  - `sigmaStdErr`
  - `muStdErr`
  - `varargin`
- **Outputs**:
  - `phandles`

#### tsEvaPlotReturnLevelsGEVFromAnalysisObj
- **Signature**: `phandles = tsEvaPlotReturnLevelsGEVFromAnalysisObj( nonStationaryEvaParams, timeIndex, varargin )`
- **Description**: nonStationaryEvaParams, stationaryTransformData are as these returned by function nonStationaryEvaJRCApproach timeIndex is the index at which the time varying analysis should be estimated.
- **Category**: Plotting & Visualization
- **Inputs**:
  - `nonStationaryEvaParams`
  - `timeIndex`
  - `varargin`
- **Outputs**:
  - `phandles`

#### tsEvaPlotReturnLevelsGPD
- **Signature**: `[phandles, returnPeriods, returnLevels, retrunLevelsErrs] = tsEvaPlotReturnLevelsGPD( epsilon, sigma, threshold, epsilonStdErr, sigmaStdErr, thresholdStdErr, nPeaks, timeHorizonInYears, varargin  )`
- **Category**: Plotting & Visualization
- **Inputs**:
  - `epsilon`
  - `sigma`
  - `threshold`
  - `epsilonStdErr`
  - `sigmaStdErr`
  - `thresholdStdErr`
  - `nPeaks`
  - `timeHorizonInYears`
  - `varargin`
- **Outputs**:
  - `phandles`
  - `returnPeriods`
  - `returnLevels`
  - `retrunLevelsErrs`

#### tsEvaPlotReturnLevelsGPDFromAnalysisObj
- **Signature**: `[phandles, returnPeriods, returnLevels, retrunLevelsErrs] = tsEvaPlotReturnLevelsGPDFromAnalysisObj( nonStationaryEvaParams, timeIndex, varargin )`
- **Description**: nonStationaryEvaParams, stationaryTransformData are as these returned by function nonStationaryEvaJRCApproach timeIndex is the index at which the time varying analysis should be estimated.
- **Category**: Plotting & Visualization
- **Inputs**:
  - `nonStationaryEvaParams`
  - `timeIndex`
  - `varargin`
- **Outputs**:
  - `phandles`
  - `returnPeriods`
  - `returnLevels`
  - `retrunLevelsErrs`

#### tsEvaPlotSeasonalityGev
- **Signature**: `phandles = tsEvaPlotSeasonalityGev(extremesRange, referenceYear, timeStamps, epsilon, sigma, mu, monthlyMaxIndexes, series, trend, stddev, varargin)`
- **Description**: tsEvaPlotSeasonalityGev plots a single year of data adding the series of monthly maxima renormalized to the considered year and superimposed to the series. I furtherly plots the time varying location and scale parameters.
- **Category**: Plotting & Visualization
- **Inputs**:
  - `extremesRange`
  - `referenceYear`
  - `timeStamps`
  - `epsilon`
  - `sigma`
  - `mu`
  - `monthlyMaxIndexes`
  - `series`
  - `trend`
  - `stddev`
  - `varargin`
- **Outputs**:
  - `phandles`

#### tsEvaPlotSeasonalityGevFromAnalysisObj
- **Signature**: `phandles = tsEvaPlotSeasonalityGevFromAnalysisObj( extremesRange, referenceYear, nonStationaryEvaParams, stationaryTransformData, varargin )`
- **Category**: Plotting & Visualization
- **Inputs**:
  - `extremesRange`
  - `referenceYear`
  - `nonStationaryEvaParams`
  - `stationaryTransformData`
  - `varargin`
- **Outputs**:
  - `phandles`

#### tsEvaPlotSeriesTrendStdDev
- **Signature**: `phandles = tsEvaPlotSeriesTrendStdDev( timeStamps, series, trend, stdDev, varargin )`
- **Category**: Plotting & Visualization
- **Inputs**:
  - `timeStamps`
  - `series`
  - `trend`
  - `stdDev`
  - `varargin`
- **Outputs**:
  - `phandles`

#### tsEvaPlotSeriesTrendStdDevFromAnalyisObj
- **Signature**: `phandles = tsEvaPlotSeriesTrendStdDevFromAnalyisObj( nonStationaryEvaParams, stationaryTransformData, varargin )`
- **Category**: Plotting & Visualization
- **Inputs**:
  - `nonStationaryEvaParams`
  - `stationaryTransformData`
  - `varargin`
- **Outputs**:
  - `phandles`

#### tsEvaPlotTransfToStat
- **Signature**: `phandles = tsEvaPlotTransfToStat( timeStamps, statSeries, srsmean, stdDev, thirdMom, fourthMom, varargin )`
- **Description**: stdThirdMom = third root of the third statistical momentum stdFouthMom = fourth root of the fourth statistical momentum
- **Category**: Plotting & Visualization
- **Inputs**:
  - `timeStamps`
  - `statSeries`
  - `srsmean`
  - `stdDev`
  - `thirdMom`
  - `fourthMom`
  - `varargin`
- **Outputs**:
  - `phandles`

#### tsEvaPlotTransfToStatFromAnalysisObj
- **Signature**: `phandles = tsEvaPlotTransfToStatFromAnalysisObj( nonStationaryEvaParams, stationaryTransformData, varargin )`
- **Description**: std3mom = nthroot(stationaryTransformData.statSer3Mom, 3.); std4mom = nthroot(stationaryTransformData.statSer4Mom, 4.);
- **Category**: Plotting & Visualization
- **Inputs**:
  - `nonStationaryEvaParams`
  - `stationaryTransformData`
  - `varargin`
- **Outputs**:
  - `phandles`

#### tsLcSubplotManager
- **Signature**: `obj = tsLcSubplotManager(N, M, varargin)`
- **Category**: Plotting & Visualization
- **Inputs**:
  - `N`
  - `M`
  - `varargin`
- **Outputs**:
  - `obj`

#### tsPlotBivarReturnPeriod
- **Signature**: ` [copulaAnalysis,hFig1] = tsPlotBivarReturnPeriod(copulaAnalysis,axxArray,varargin)`
- **Description**: tsCopulaPlotJointReturnPeriod plotting of multivariate return periods
- **Category**: Plotting & Visualization
- **Inputs**:
  - `copulaAnalysis`
  - `axxArray`
  - `varargin`
- **Outputs**:
  - `copulaAnalysis`
  - `hFig1`

#### tsPlotSeriesPotGPDRetLevFromAnalysisObj
- **Signature**: `phandles = tsPlotSeriesPotGPDRetLevFromAnalysisObj( nonStationaryEvaParams, stationaryTransformData, varargin )`
- **Description**: computing the return levels
- **Category**: Plotting & Visualization
- **Inputs**:
  - `nonStationaryEvaParams`
  - `stationaryTransformData`
  - `varargin`
- **Outputs**:
  - `phandles`

### Data Preparation & Utilities
#### cvineOrder
- **Signature**: `order = cvineOrder(alpha)`
- **Description**: Root-first C-vine order from a Gumbel-θ matrix (alpha). Heuristic: pick node with max total |tau| as root, then greedily add the next.
- **Category**: Data Preparation & Utilities
- **Inputs**:
  - `alpha`
- **Outputs**:
  - `order`

#### tsEVstatistics
- **Signature**: `[EVmeta,EVdata,isValid] = tsEVstatistics(pointData, varargin)`
- **Description**: Evangelos Voukouvalas, Michalis Vousdoukas 2015 gevMaxima can be annual or monthly. annual by default
- **Category**: Data Preparation & Utilities
- **Inputs**:
  - `pointData`
  - `varargin`
- **Outputs**:
  - `EVmeta`
  - `EVdata`
  - `isValid`

#### tsEasyParseNamedArgs
- **Signature**: `[argStruct] = tsEasyParseNamedArgs(args, argStruct)`
- **Category**: Data Preparation & Utilities
- **Inputs**:
  - `args`
  - `argStruct`
- **Outputs**:
  - `argStruct`

#### tsEmpirical
- **Signature**: `[ C ] = tsEmpirical( U )`
- **Description**: tsEmpirical Empirical copula [ C ] = tsEmpirical( U ) retuns a variable C containing the empirical copula corresponding with the U variable including the uniform variates
- **Category**: Data Preparation & Utilities
- **Inputs**:
  - `U`
- **Outputs**:
  - `C`

#### tsEnsemble
- **Signature**: `[nonStationaryEvaParamsEns, stationaryTransformDataEns]  = tsEnsemble( nonStatEvaParamsArray, stationaryTransformDataArray )`
- **Description**: calls tsEnsembleEvaParams and tsEnsembleStatTransfData. ASSUMING THAT ALL THE SERIES HAVE THE SAME LENGTH AND UNIFORM TIME STAMPS
- **Category**: Data Preparation & Utilities
- **Inputs**:
  - `nonStatEvaParamsArray`
  - `stationaryTransformDataArray`
- **Outputs**:
  - `nonStationaryEvaParamsEns`
  - `stationaryTransformDataEns`

#### tsEnsembleEvaParams
- **Signature**: `nonStatEvaParamsEnsemble = tsEnsembleEvaParams( nonStatEvaParamsArray )`
- **Description**: From a cell array of non nonStatEvaParams computes the average nonStatEvaParamsArray is a cell array of nonStatEvaParams, the object type returned by tsEvaNonStationary ASSUMING THAT ALL THE SERIES HAVE THE SAME LENGTH AND UNIFORM TIME STAMPS
- **Category**: Data Preparation & Utilities
- **Inputs**:
  - `nonStatEvaParamsArray`
- **Outputs**:
  - `nonStatEvaParamsEnsemble`

#### tsEnsembleStatTransfData
- **Signature**: `[ stationaryTransformDataEnsemble ] = tsEnsembleStatTransfData( stationaryTransformDataArray )`
- **Description**: From a cell array of non stationaryTransformData computes the average stationaryTransformDataArray is a cell array of stationaryTransformData, the transformation object type returned by tsEvaNonStationary ASSUMING THAT ALL THE SERIES HAVE THE SAME LENGTH AND UNIFORM TIME STAMPS
- **Category**: Data Preparation & Utilities
- **Inputs**:
  - `stationaryTransformDataArray`
- **Outputs**:
  - `stationaryTransformDataEnsemble`

#### tsEstimateConfidenceIntervalOfRL
- **Signature**: `[ lowCI, highCI ] = tsEstimateConfidenceIntervalOfRL( rl, stdErr, p )`
- **Description**: % tsEstimateConfidenceIntervalOfRL: estimates the confidence interval of the return level (rl), assuming that the uncertainty distribution around rl is a lognormal with standard deviation stdErr %
- **Category**: Data Preparation & Utilities
- **Inputs**:
  - `rl`
  - `stdErr`
  - `p`
- **Outputs**:
  - `lowCI`
  - `highCI`

#### tsEvaFillSeries
- **Signature**: `[ filledTimeStamps, filledSeries, dt ] = tsEvaFillSeries( timeStamps, series )`
- **Description**: ensuring monotonic time vector and no dubplicates. In case of dubplicate time stamps the highest value is considered
- **Category**: Data Preparation & Utilities
- **Inputs**:
  - `timeStamps`
  - `series`
- **Outputs**:
  - `filledTimeStamps`
  - `filledSeries`
  - `dt`

#### tsEvaGetReturnPeriodOfLevelGEV
- **Signature**: `[retPer, exceedProb] = tsEvaGetReturnPeriodOfLevelGEV( epsilon, sigma, mu, retLev)`
- **Description**: GEV
- **Category**: Data Preparation & Utilities
- **Inputs**:
  - `epsilon`
  - `sigma`
  - `mu`
  - `retLev`
- **Outputs**:
  - `retPer`
  - `exceedProb`

#### tsEvaGetReturnPeriodOfLevelGPD
- **Signature**: `[retPer, exceedProb] = tsEvaGetReturnPeriodOfLevelGPD( epsilon, sigma, pPeak, retLev)`
- **Category**: Data Preparation & Utilities
- **Inputs**:
  - `epsilon`
  - `sigma`
  - `pPeak`
  - `retLev`
- **Outputs**:
  - `retPer`
  - `exceedProb`

#### tsEvaGetTimeStep
- **Signature**: `dt = tsEvaGetTimeStep( times )`
- **Category**: Data Preparation & Utilities
- **Inputs**:
  - `times`
- **Outputs**:
  - `dt`

#### tsEvaSampleData
- **Signature**: `pointData = tsEvaSampleData(ms, varargin)`
- **Description**: args.pcts = [50 70 85:2:99 99.1:0.1:99.9 99.91:0.01:99.99];
- **Category**: Data Preparation & Utilities
- **Inputs**:
  - `ms`
  - `varargin`
- **Outputs**:
  - `pointData`

#### tsGet
- **Signature**: `rho = tsGet()`
- **Category**: Data Preparation & Utilities
- **Outputs**:
  - `rho`

#### tsGetNumberPerYear
- **Signature**: `nperYear=tsGetNumberPerYear(ms,locs)`
- **Description**: function nperYear=tsGetNumberPerYear(ms,locs) Gives number of events per year from a time series ms (sd and values vertical vectors) and sporadic indices of events (locs) Michalis Vousdoukas 2015
- **Category**: Data Preparation & Utilities
- **Inputs**:
  - `ms`
  - `locs`
- **Outputs**:
  - `nperYear`

#### tsGetPOT
- **Signature**: `[POTdata]=tsGetPOT(ms,pcts,desiredEventsPerYear, varargin)`
- **Description**: function POTdata=tsGetPOT(ms,pcts,desiredEventsPerYear) Gets POT using an automatic threshold such that the mean number of events per year is equal to desiredEventsPerYear
- **Category**: Data Preparation & Utilities
- **Inputs**:
  - `parameter`: ms  data in two columns sdate and values. IT IS ASSUMED
  - `parameter`: pcts    vestor of percentiles tested
  - `parameter`: desiredEventsPerYear    mean number of events per year
  - `parameter`: Michalis Vousdoukas, Evangelos Voukouvalas, Lorenzo Mentaschi 2015
- **Outputs**:
  - `POTdata`

#### tsGetReturnPeriodOfLevel
- **Signature**: `[ myRetPeriod, myRetPeriodCISup, myRetPeriodCIInf ] = tsGetReturnPeriodOfLevel( retPeriod, retLevel, retLevError, myLevel, varargin )`
- **Description**: given a list of retPeriod and corresponding retLevel with error, estimates the return period for a level myLevel, and the related confidence interval.
- **Category**: Data Preparation & Utilities
- **Inputs**:
  - `retPeriod`
  - `retLevel`
  - `retLevError`
  - `myLevel`
  - `varargin`
- **Outputs**:
  - `myRetPeriod`
  - `myRetPeriodCISup`
  - `myRetPeriodCIInf`

#### tsInterp1Extrap
- **Signature**: `nvals = tsInterp1Extrap(X,V,Xq,logExtrap)`
- **Description**: function nvals=interp1Extrap(Tr,values,RP,logExtrap) This function interpolates a series and extrapolates linearly or exponeniantly if necessary. Can also handle matrices (only for dependent variable values) Will also sort the Tr series if not monotonically increasing Tr-     dependent variable
- **Category**: Data Preparation & Utilities
- **Inputs**:
  - `X`
  - `V`
  - `Xq`
  - `logExtrap`
- **Outputs**:
  - `nvals`

#### tsLinearExtrapolation
- **Signature**: `[nx,ny,addedy]=tsLinearExtrapolation(x,y,extrax,npoints)`
- **Description**: Linearly extrapolates time series x,y to points extrax, according to the slope of the first/last number of points of the series (npoints) Michalis Vousdoukas 2016
- **Category**: Data Preparation & Utilities
- **Inputs**:
  - `x`
  - `y`
  - `extrax`
  - `npoints`
- **Outputs**:
  - `nx`
  - `ny`
  - `addedy`

#### tsLoglogExtrapolation
- **Signature**: `[nx,ny,addedy]=tsLoglogExtrapolation(x,y,extrax,npoints)`
- **Description**: Extrapolates exponentiallytime series x,y to points extrax, according to the slope of the first/last number of points of the series (npoints) Michalis Vousdoukas 2016
- **Category**: Data Preparation & Utilities
- **Inputs**:
  - `x`
  - `y`
  - `extrax`
  - `npoints`
- **Outputs**:
  - `nx`
  - `ny`
  - `addedy`

#### tsPseudoObservations
- **Signature**: `[ U ] = tsPseudoObservations( X )`
- **Category**: Data Preparation & Utilities
- **Inputs**:
  - `X`
- **Outputs**:
  - `U`

#### tsRankmax
- **Signature**: `[ R ] = tsRankmax( X )`
- **Description**: tsRankmax Returns vector of one-based ranks for each element
- **Category**: Data Preparation & Utilities
- **Inputs**:
  - `X`
- **Outputs**:
  - `R`

#### tsRemoveConstantSubseries
- **Signature**: `cleaned_series = tsRemoveConstantSubseries( srs, stackedValuesCount )`
- **Category**: Data Preparation & Utilities
- **Inputs**:
  - `srs`
  - `stackedValuesCount`
- **Outputs**:
  - `cleaned_series`

#### tsRoundSDate
- **Signature**: `[sdround,dvecm,sdunique,dvunique]=tsRoundSDate(sd,sd_precision)`
- **Description**: function [sdround,dvec,sdunique,dvunique]=tsRoundSDate(sd,sd_precision) ROunds sd values according to desired precision, years, months, days, etc : 1-4 Michalis Vousdoukas 2015
- **Category**: Data Preparation & Utilities
- **Inputs**:
  - `sd`
  - `sd_precision`
- **Outputs**:
  - `sdround`
  - `dvecm`
  - `sdunique`
  - `dvunique`

#### tsSameValuesSegmentation
- **Signature**: `[inds,rinds]=tsSameValuesSegmentation(iii,varargin)`
- **Description**: function [inds,rinds]=tsSameValuesSegmentation(iii,val) separates segments of same value val
- **Category**: Data Preparation & Utilities
- **Inputs**:
  - `iii`: the series
  - `parameter`: optional variable-the value
- **Outputs**:
  - `inds`: cell with the values of each continuous set
  - `rinds`: cell with indexes of each continuous set
  - `output`: Michalis Vousdoukas 2009

#### tsTimeSeriesToPointData
- **Signature**: `pointData = tsTimeSeriesToPointData( ms, potThreshold, potThresholdError )`
- **Description**: tsTimeSeriesToPointData: given a ms produces a structure pointData like the one produced by tsEvaSampleData
- **Category**: Data Preparation & Utilities
- **Inputs**:
  - `ms`
  - `potThreshold`
  - `potThresholdError`
- **Outputs**:
  - `pointData`

#### tsYear
- **Signature**: `year = tsYear(timeStamp)`
- **Category**: Data Preparation & Utilities
- **Inputs**:
  - `timeStamp`
- **Outputs**:
  - `year`

