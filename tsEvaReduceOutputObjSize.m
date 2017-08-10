function [ redNonStatEvaParams, redStatTransData ] = tsEvaReduceOutputObjSize( nonStationaryEvaParams, stationaryTransformData, newTimeStamps )
  origTimeStamps = stationaryTransformData.timeStamps;
  tsIndxs = knnsearch(origTimeStamps, newTimeStamps);

  redStatTransData = stationaryTransformData;
  %redStatTransData.stationarySeries = redStatTransData.stationarySeries(tsIndxs);
  %redStatTransData.nonStatSeries = redStatTransData.nonStatSeries(tsIndxs);
  redStatTransData.statsTimeStamps = redStatTransData.timeStamps(tsIndxs);
  redStatTransData.trendSeries = redStatTransData.trendSeries(tsIndxs);
  redStatTransData.trendSeriesNonSeasonal = redStatTransData.trendSeriesNonSeasonal(tsIndxs);
  redStatTransData.stdDevSeries = redStatTransData.stdDevSeries(tsIndxs);
  redStatTransData.stdDevSeriesNonSeasonal = redStatTransData.stdDevSeriesNonSeasonal(tsIndxs);
  redStatTransData.statSer3Mom = redStatTransData.statSer3Mom(tsIndxs);
  redStatTransData.statSer4Mom = redStatTransData.statSer4Mom(tsIndxs);
  
  redNonStatEvaParams = nonStationaryEvaParams;
  if ~isempty(redNonStatEvaParams(1).parameters)
    redNonStatEvaParams(1).parameters.sigma = redNonStatEvaParams(1).parameters.sigma(tsIndxs);
    redNonStatEvaParams(1).parameters.mu = redNonStatEvaParams(1).parameters.mu(tsIndxs);
    redNonStatEvaParams(1).paramErr.sigmaErrFit = redNonStatEvaParams(1).paramErr.sigmaErrFit(tsIndxs);
    redNonStatEvaParams(1).paramErr.sigmaErrTransf = redNonStatEvaParams(1).paramErr.sigmaErrTransf(tsIndxs);
    redNonStatEvaParams(1).paramErr.sigmaErr = redNonStatEvaParams(1).paramErr.sigmaErr(tsIndxs);
    redNonStatEvaParams(1).paramErr.muErrFit = redNonStatEvaParams(1).paramErr.muErrFit(tsIndxs);
    redNonStatEvaParams(1).paramErr.muErrTransf = redNonStatEvaParams(1).paramErr.muErrTransf(tsIndxs);
    redNonStatEvaParams(1).paramErr.muErr = redNonStatEvaParams(1).paramErr.muErr(tsIndxs);
    redNonStatEvaParams(1).objs = [];
  end
  if ~isempty(redNonStatEvaParams(2).parameters)
    redNonStatEvaParams(2).parameters.sigma = redNonStatEvaParams(2).parameters.sigma(tsIndxs);
    redNonStatEvaParams(2).parameters.threshold = redNonStatEvaParams(2).parameters.threshold(tsIndxs);
    redNonStatEvaParams(2).paramErr.sigmaErrFit = redNonStatEvaParams(2).paramErr.sigmaErrFit(tsIndxs);
    redNonStatEvaParams(2).paramErr.sigmaErrTransf = redNonStatEvaParams(2).paramErr.sigmaErrTransf(tsIndxs);
    redNonStatEvaParams(2).paramErr.sigmaErr = redNonStatEvaParams(2).paramErr.sigmaErr(tsIndxs);
    redNonStatEvaParams(2).paramErr.thresholdErrFit = redNonStatEvaParams(2).paramErr.thresholdErrFit;
    redNonStatEvaParams(2).paramErr.thresholdErrTransf = redNonStatEvaParams(2).paramErr.thresholdErrTransf(tsIndxs);
    redNonStatEvaParams(2).paramErr.thresholdErr = redNonStatEvaParams(2).paramErr.thresholdErr(tsIndxs);
    redNonStatEvaParams(2).objs = [];
  end
  
end

