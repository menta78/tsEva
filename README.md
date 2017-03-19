# tsEva #

This MATLAB toolbox contains an implementation of the Transformed-Stationary (TS) methodology for non-stationary EVA as described in Mentaschi et al. (2016). In synthesis this approach consists in (i) transforming a non-stationary time series into a stationary one to which the stationary EVA theory can be applied; and (ii) reverse-transforming the result into a non-stationary extreme value distribution.

This toolbox is free of external dependencies, and contains calls for both the statistical analysis and the graphical rendering of the results.

A version working on MATLAB R2014b is available at the branch 0.1_R2014b of this project, and can be downloaded [here](https://github.com/menta78/tsEva/archive/0.1_R2014b.zip).


# How to start #

Subdirectory "test" contains sample scripts that illustrate the features of the toolbox. The main sample scripts are:

* test/sampleGenerateSeriesEVAGraphs.m: this script analyzes a time series of residual water levels modeled at the Hebrides islands on a time horizon of 30 years. 
The minimum time window for the detection of non-stationarity is 6 years (i.e. the time window below which the statistics of the series are considered stationary). 
The script estimates the non stationary GEV and GPD, and plots the results. 
The following plots are produced: the 2d plot of the non-stationary GEV and GPD, the 3d plot of the GEV, the return levels for both. 
The whole analysis is then repeated including the estimation of the seasonal extremes.

* test/sampleGenerateSeriesEVAGraphs_ciPercentile.m: this script analyzes the same time series as the previous script, but estimates the time-varying confidence interval using a moving 98.5 percentile, instead of a moving standard deviation.
**Using the moving percentile instead of the moving standard deviation is useful because the moving percentile is more sensitive to changes of the extremes**. 
The drawback is that the confidence interval tends to be broader.

* test/sampleSPISeries.m: this script analyzes a time series of projected SPI (Standardized Precipitation Index).
In this series the peaks are distant at least 5 months one from the others. 
Therefore the concept of annual maxima is meaningless, and the sole GPD approach is applied. 
The minimum time window for the detection of non-stationarity is 50 years (the whole time horizon of the time series is 130 years).


## References ##
[Mentaschi, L., Vousdoukas, M., Voukouvalas, E., Sartini, L., Feyen, L., Besio, G., and Alfieri, L.: The transformed-stationary approach: a generic and simplified methodology for non-stationary extreme value analysis, Hydrol. Earth Syst. Sci., 20,3527-3547, doi:10.5194/hess-20-3527-2016, 2016](http://www.hydrol-earth-syst-sci.net/20/3527/2016/)
