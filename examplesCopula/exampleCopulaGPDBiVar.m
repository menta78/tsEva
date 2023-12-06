% an example where dependence structure of two different variables (storm surges at two different
% locations (Marshall-North and Marshall-South stations) are tested againt
% each other using different copula functions
% if timeWindow is selected larger than duration of time series (i.e., ~ 30
% years), a stationary (or time-invariant) copula will be used (i.e, one
% copula only). If the timeWindow is selected to be less than duration of
% time series (i.e., 10 years), a time-varying copula will be plotted
% Copula could be of type Gaussian, t, Frank, Gumbel and Clayton

% M. H. Bahmanpour, 2023

close all;
clc
clearvars;
addpath('../');

load('testCopulaCompoundGPDbivariate')  %data is SS (storm surge) from Marshall Islan in two different locations

thresholdPercentiles=[99,99]; %threshold levels in each series for sampling of data
ciPercentile = 99;  %to be used in tsEvaNonStationary
potPercentiles=[99]; %o be used in tsEvaNonStationary; better be set to only one value but a range is also possible

timeWindow = 365*15; %if time window (in days) is equal or larger than
% length of time series, copula will be of stationary type; a very small
% time window may result in error (e.g., less than 3 years)

minPeakDistanceInDaysforjointpeaks=[3,3]; %minimum distance between monovariate peaks of each series
maxDistanceMultivariatePeaksInDays=[3]; %can either take one value or has to have a format and size matching size(nchoosek([1:numvar],2),1) where numvar is 2 in bivariate case, 3 in trivariate case, and so on
minPeakDistanceInDaysfornonstationary=[3,3]; %used for sampling peaks inside tsevavnonstationary

[copulaAnalysis] = tsCopulaCompoundGPD(timeAndSeries1(:,1), ...
    [timeAndSeries1(:,2),timeAndSeries2(:,2)], ...
    thresholdPercentiles, ...
    minPeakDistanceInDaysforjointpeaks, ...
    maxDistanceMultivariatePeaksInDays, ...
    'copulaFamily','Gaussian',...
    'transfType','trendlinear','timeWindow',timeWindow,...
    'ciPercentile',ciPercentile,'potPercentiles',potPercentiles,'minPeakDistanceInDays',minPeakDistanceInDaysfornonstationary);

nResample=10000; %used for Monte-Carlo simulation

[resampleLevel, resampleProb] = tsCopulaCompoundGPDMontecarlo(copulaAnalysis,...
    nResample,'timeIndex',1);

if iscell(copulaAnalysis.jointExtremeMonovariateProb)
    figHnd = tsCopulaTimeVaryingPlot(resampleLevel, copulaAnalysis, 'xlbl', 'Marshall-north', 'ylbl', 'Marshall-south','numberofverticalpanels',3,'numberofhorizontalpanels',3);
    tsCopulaPlotJointReturnPeriod(copulaAnalysis,'marginalDistributions','gp','plotType','AND')% plot type: AND, OR

else
    peakIndices=copulaAnalysis.jointExtremeIndices;
    yMax=copulaAnalysis.jointExtremes; %non-stationary peaks
    tMax=copulaAnalysis.jointExtremeTimeStamps;
    thresholdC=copulaAnalysis.thresholdPotNS; %non-stationary threshold parameters
    figure

    plot(datetime(datevec(timeAndSeries1(:,1))),timeAndSeries1(:,2))
    hold on
    plot(datetime(datevec(timeAndSeries2(:,1))),timeAndSeries2(:,2))

    plot(datetime(datevec(timeAndSeries1(:,1))),thresholdC(:,1))
    plot(datetime(datevec(timeAndSeries2(:,1))),thresholdC(:,2))
    plot(datetime(datevec(tMax(:,1))),yMax(:,1),'.r')
    plot(datetime(datevec(tMax(:,2))),yMax(:,2),'.k')

    legend('Series1','Series2','Threshold-1','Threshold-2','Peaks-1','Peaks-2')
    title('Non-stationary series')

    figHnd = tsCopulaTimeVaryingPlot(resampleLevel, copulaAnalysis, 'xlbl', 'Marshall-north', 'ylbl', 'Marshall-south','numberofverticalpanels',1,'numberofhorizontalpanels',1);

    tsCopulaPlotJointReturnPeriod(copulaAnalysis,'marginalDistributions','gp','plotType','AND')% plot type: AND, OR
end