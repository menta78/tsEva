# tsEva #

This MATLAB toolbox implements the Transformed-Stationary (TS) methodology for non-stationary Extreme Value Analysis (EVA), originally described in Mentaschi et al. (2016) for the univariate case and extended to the multivariate case in Bahmanpour et al. (2026).

## Univariate Analysis (tsEVA 1.0)

The univariate approach consists in (i) transforming a non-stationary time series into a stationary one to which stationary EVA theory can be applied; and (ii) reverse-transforming the result into a non-stationary extreme value distribution. This covers non-stationary GEV and GPD estimation, return level computation, and seasonal extreme analysis.

## Multivariate Analysis (tsEVA 2.0)

tsEVA 2.0 extends the framework to multivariate compound events through time-varying copulas (Bahmanpour et al., 2026). The key steps are:

1. **Marginal transformation:** each univariate marginal is transformed to stationarity using the TS approach of Mentaschi et al. (2016).
2. **Dependence modeling:** the stationary margins are coupled via time-dependent copulas, capturing evolving joint dependence structures under non-stationarity.
3. **Montecarlo simulation:** allows Montecarlo-based analysis of hazard, quantificaion of the goodness of fit, estimation of multivariate return levels.
4. **Multivariate Gumbel C-vine:** a flexible Archimedean vine copula structure is available for modeling higher-dimensional dependence, replacing Gaussian-only assumptions with a realistic, tail-sensitive framework.


This toolbox is free of external dependencies and contains calls for both statistical analysis and graphical rendering of results.


## How to start ##

Subdirectory `test` contains sample scripts that illustrate the features of the toolbox.

### Univariate examples

* `examplesMonovariate/exampleGenerateSeriesEVAGraphs.m`: analyzes a time series of residual water levels modeled at the Hebrides islands on a 30-year time horizon. The minimum time window for the detection of non-stationarity is 6 years. The script estimates non-stationary GEV and GPD, and plots 2D and 3D distributions and return levels. The analysis is then repeated including seasonal extremes.

* `examplesMonovariate/exampleGenerateSeriesEVAGraphs_ciPercentile.m`: same analysis as above, but estimates time-varying confidence intervals using a moving 98.5 percentile instead of a moving standard deviation. **Using the moving percentile is useful because it is more sensitive to changes in the extremes.** The drawback is that the confidence interval tends to be broader.

* `examplesMonovariate/exampleSPISeries.m`: analyzes a time series of projected SPI (Standardized Precipitation Index) with peaks at least 5 months apart. The concept of annual maxima is meaningless here, so only the GPD approach is applied. The minimum time window for non-stationarity detection is 50 years (total horizon: 130 years).

### Multivariate examples
The multivariate examples allow to reconstruct the case studies analyzed in Bahmanpour et al. 2026, plus a couple of other examples. In particular:

* `examplesCopula\caseStudy01.m`: reproduces figure 2 of Bahmanpour et al. 2026. A POT Gumbel bivariate analysis of river discharge and significant wave height.

* `examplesCopula\caseStudy02.m`: reproduces figure 3 of Bahmanpour et al. 2026. A POT Gumbel trivariate analysis of extreme significant wave height in three locations.

* `examplesCopula\caseStudy03.m`: reproduces figure 4 of Bahmanpour et al. 2026. A GEV Gumbel bivariate analysis of extreme temperature and SPEI.

* the other three scripts in examplesCopula exemplify how the different functions of tsEva can be used separately to build customized a more customized analysis.

## References ##

* [Mentaschi, L., Vousdoukas, M., Voukouvalas, E., Sartini, L., Feyen, L., Besio, G., and Alfieri, L.: The transformed-stationary approach: a generic and simplified methodology for non-stationary extreme value analysis, Hydrol. Earth Syst. Sci., 20, 3527–3547, doi:10.5194/hess-20-3527-2016, 2016](http://www.hydrol-earth-syst-sci.net/20/3527/2016/)

* [Bahmanpour, M. H., Tilloy, A., Vousdoukas, M., Federico, I., Coppini, G., Feyen, L., and Mentaschi, L.: Transformed-stationary EVA 2.0: a generalized framework for non-stationary multivariate extremes analysis, Hydrol. Earth Syst. Sci., 30, 2301–2314, https://doi.org/10.5194/hess-30-2301-2026, 2026]([http://www.hydrol-earth-syst-sci.net/20/3527/2016/](https://hess.copernicus.org/articles/30/2301/2026/))


## AI Assistant ##

An agentic assistant (Tessa M) is available to help users navigate the toolbox and implement their analysis. 
Tessa M is available as a [custom GPT](https://chatgpt.com/g/g-69a72fb764048191ae68e2201f7741a3-tessa-m) and as a github customized copilot.

Tessa M can be woken up in VS Code github copilot extension after opening this repository. In the github chat use this prompt:
"Hello Tessa-M, please wake up and read carefully your prompts and instructions in tsEva/.github/copilot/agents/tessa-M.json"

The best way to use Tessa M on Vs Code is by using the MATLAB extension for Vs Code provided by Mathworks, and letting Tessa M act in agent mode.

Tessa M assists with running analyses, interpreting outputs, and troubleshooting common errors in tsEVA 2.0.
