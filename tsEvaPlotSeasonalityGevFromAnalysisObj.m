function phandles = tsEvaPlotSeasonalityGevFromAnalysisObj( extremesRange, referenceYear, nonStationaryEvaParams, stationaryTransformData, varargin )
    timeStamps = stationaryTransformData.timeStamps;
    nonStatSrs = stationaryTransformData.nonStatSeries;
    epsilon = nonStationaryEvaParams(1).parameters.epsilon;
    sigma = nonStationaryEvaParams(1).parameters.sigma;
    mu = nonStationaryEvaParams(1).parameters.mu;
    monthlyMaxIndexes = nonStationaryEvaParams(1).objs.monthlyMaxIndexes;
    series = nonStatSrs;
    trend = stationaryTransformData.trendSeriesNonSeasonal;
    stddev = stationaryTransformData.stdDevSeriesNonSeasonal;
    phandles = tsEvaPlotSeasonalityGev(extremesRange, referenceYear, timeStamps,...
        epsilon, sigma, mu, monthlyMaxIndexes, series, trend, stddev, varargin);
end

