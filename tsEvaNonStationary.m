function [nonStationaryEvaParams, stationaryTransformData, isValid] = tsEvaNonStationary( timeAndSeries, timeWindow, varargin )
% tsEvaNonStationary:
% performs the TS EVA analysis as described by Mentaschi et al 2016.
%
% input parameters:
%   timeAndSeries: array with shape nx2, with the time stamps in the first
%            column and the values in the second.
%   timeWindow: time window for the transformation expressed in days.
%
%   (some) label parameters:
%         transfType: can assume values 
%                    1) 'trend': long term variability. The trend is computed 
%                            with a running mean, the ci with the running standard deviation.
%                    2) 'seasonal': long term + seasonal variability. The trend is computed 
%                            with a running mean, the ci with the running
%                            standard deviation.
%                    3) 'trendCIPercentile': long term variability. The trend is computed 
%                            with a running mean, the ci with the running xx percentile.
%                            Using this option the argument ciPercentile is
%                            mandatory.
% 
%
%% sample calls
% nonStatEvaParams = tsEvaNonStationary(ms, timeWindow, 'potPercentiles',[95], 'minPeakDistanceInDays', 3)
%     samples POT data using a fixed 95 percentile threshold, with peaks at
%     a minimum distance of 3 days, looking for a threshold so that we have an average
%     of 5 events every year.
% nonStatEvaParams = tsEvaNonStationary(ms, timeWindow, 'minPeakDistanceInDays', 3, 'desiredeventsperyear', 6)
%     samples POT data looking for a threshold so that we have an average
%     of 6 events every year.
% nonStatEvaParams = tsEvaNonStationary(ms, timeWindow, 'minPeakDistanceInDays', 3, 'trasftype', 'trendCIPercentile', 'ciPercentile', 99)
%     for the transformation uses instead of the moving standard deviation,
%     the moving 99th percentile.
%% %%%%%%%%%%%%%

args.transfType = 'trend';
args.minPeakDistanceInDays = -1;
args.ciPercentile = NaN;
args.potEventsPerYear = 5;
args.evdType = {'GEV', 'GPD'};
args.gevType = 'GEV'; % can be 'GEV' or 'Gumbel'
args = tsEasyParseNamedArgs(varargin, args);
minPeakDistanceInDays = args.minPeakDistanceInDays;
ciPercentile = args.ciPercentile;
transfType = args.transfType;
evdType = args.evdType;
gevType = args.gevType;

if ~( strcmpi(transfType, 'trend') || strcmpi(transfType, 'seasonal') || strcmpi(transfType, 'trendCIPercentile') || strcmpi(transfType, 'seasonalCIPercentile') )
  error('nonStationaryEvaJRCApproach: transfType can be in (trend, seasonal, trendCIPercentile)');
end
if minPeakDistanceInDays == -1
  error('label parameter ''minPeakDistanceInDays'' must be set')
end
   
nonStationaryEvaParams = [];
stationaryTransformData = [];

timeStamps = timeAndSeries(:, 1);
series = timeAndSeries(:, 2);

if strcmpi(transfType, 'trend')
  disp('evalueting long term variations of extremes');
  trasfData = tsEvaTransformSeriesToStationaryTrendOnly( timeStamps, series, timeWindow, varargin{:} );
  gevMaxima = 'annual';
  potEventsPerYear = 5;
elseif strcmpi(transfType, 'seasonal')
  disp('evalueting long term an seasonal variations of extremes');
  trasfData = tsEvaTransformSeriesToStationaryMultiplicativeSeasonality( timeStamps, series, timeWindow, varargin{:}  );
  gevMaxima = 'monthly';
  potEventsPerYear = 12;
elseif strcmpi(transfType, 'trendCIPercentile') 
  if isnan(ciPercentile)
    error('For trendCIPercentile transformation the label parameter ''cipercentile'' is mandatory');
  end
  disp(['evalueting long term variations of extremes using the ' num2str(ciPercentile) 'th percentile']);
  trasfData = tsEvaTransformSeriesToStationaryTrendOnly_ciPercentile( timeStamps, series, timeWindow, ciPercentile, varargin{:} );
  gevMaxima = 'annual';
  potEventsPerYear = 5;
elseif strcmpi(transfType, 'seasonalCIPercentile') 
  if isnan(ciPercentile)
    error('For seasonalCIPercentile transformation the label parameter ''cipercentile'' is mandatory');
  end
  disp(['evalueting long term variations of extremes using the ' num2str(ciPercentile) 'th percentile']);
  trasfData = tsEvaTransformSeriesToStatSeasonal_ciPercentile( timeStamps, series, timeWindow, ciPercentile, varargin{:} );
  gevMaxima = 'monthly';
  potEventsPerYear = 12;
end
if args.potEventsPerYear ~= -1
  potEventsPerYear = args.potEventsPerYear;
end
ms = cat(2, trasfData.timeStamps, trasfData.stationarySeries);
%dt = trasfData.timeStamps(2) - trasfData.timeStamps(1);
dt = tsEvaGetTimeStep(trasfData.timeStamps);
minPeakDistance = minPeakDistanceInDays/dt;

%% estimating the non stationary EVA parameters
fprintf('\n');
disp('Executing stationary eva')
pointData = tsEvaSampleData(ms, 'meanEventsPerYear', potEventsPerYear, varargin{:});
evaAlphaCI = .68; % in a gaussian approximation alphaCI~68% corresponds to 1 sigma confidence
[~, eva, isValid] = tsEVstatistics(pointData, 'alphaci', evaAlphaCI, ...
       'gevmaxima', gevMaxima, 'gevType', gevType, 'evdType', evdType);
if ~isValid
  return;
end
eva(2).thresholdError = pointData.POT.thresholdError;
fprintf('\n');

% !!! Assuming a Gaussian approximation to compute the standard errors for
% the GEV parameters
if ~isempty(eva(1).parameters)
  epsilonGevX = eva(1).parameters(1);
  errEpsilonX = epsilonGevX - eva(1).paramCIs(1, 1);
  sigmaGevX = eva(1).parameters(2);
  errSigmaGevX = sigmaGevX - eva(1).paramCIs(1, 2);
  muGevX = eva(1).parameters(1, 3);
  errMuGevX = muGevX - eva(1).paramCIs(1, 3);

  fprintf('\n');
  disp('Transforming to non stationary eva ...')
  epsilonGevNS = epsilonGevX;
  errEpsilonGevNS = errEpsilonX;
  sigmaGevNS = trasfData.stdDevSeries*sigmaGevX;
  % propagating the errors on stdDevSeries and sigmaGevX to sigmaGevNs.
  % err(sigmaNs) = sqrt{ [sigmaX*err(stdDev)]^2 + [stdDev*err(sigmaX)]^2 }
  % the error on sigmaGevNs is time dependant.
  errSigmaGevFit = trasfData.stdDevSeries.*errSigmaGevX;
  errSigmaGevTransf = sigmaGevX*trasfData.stdDevError;
  errSigmaGevNS = (  errSigmaGevTransf.^2   +  errSigmaGevFit.^2  ).^.5;
  muGevNS = trasfData.stdDevSeries*muGevX + trasfData.trendSeries;
  % propagating the errors on stdDevSeries, trendSeries and sigmaGevX to muGevNS.
  % err(muNs) = sqrt{ [muX*err(stdDev)]^2 + [stdDev*err(muX)]^2 + err(trend)^2 }
  % the error on muGevNS is time dependant.
  errMuGevFit = trasfData.stdDevSeries.*errMuGevX;
  errMuGevTransf = (  (muGevX*trasfData.stdDevError).^2 + trasfData.trendError^2  ).^.5;
  errMuGevNS = (  errMuGevTransf.^2   +  errMuGevFit.^2  ).^.5;
  gevParams.epsilon = epsilonGevNS;
  gevParams.sigma = sigmaGevNS;
  gevParams.mu = muGevNS;
  if strcmpi(gevMaxima, 'annual')
    gevParams.timeDelta = 365.25;
    gevParams.timeDeltaYears = 1;
  elseif strcmpi(gevMaxima, 'monthly')
    gevParams.timeDelta = 365.25/12.;
    gevParams.timeDeltaYears = 1/12.;
  end
  
  gevParamStdErr.epsilonErr = errEpsilonGevNS;
  
  gevParamStdErr.sigmaErrFit = errSigmaGevFit;
  gevParamStdErr.sigmaErrTransf = errSigmaGevTransf;
  gevParamStdErr.sigmaErr = errSigmaGevNS;
  
  gevParamStdErr.muErrFit = errMuGevFit;
  gevParamStdErr.muErrTransf = errMuGevTransf;
  gevParamStdErr.muErr = errMuGevNS;
  
  gevObj.method = eva(1).method;
  gevObj.parameters = gevParams;
  gevObj.paramErr = gevParamStdErr;
  gevObj.stationaryParams = eva(1);
  gevObj.objs.monthlyMaxIndexes = pointData.monthlyMaxIndexes;
else
  gevObj.method = eva(1).method;
  gevObj.parameters = [];
  gevObj.paramErr = [];
  gevObj.stationaryParams = [];
  gevObj.objs.monthlyMaxIndexes = [];
end

%% estimating the non stationary GPD parameters
% !!! Assuming a Gaussian approximation to compute the standard errors for
% the GPD parameters
if ~isempty(eva(2).parameters)
  epsilonPotX = eva(2).parameters(2);
  errEpsilonPotX = epsilonPotX - eva(2).paramCIs(1, 2);
  sigmaPotX = eva(2).parameters(1);
  errSigmaPotX = sigmaPotX - eva(2).paramCIs(1, 1);
  thresholdPotX = eva(2).parameters(3);
  errThresholdPotX = eva(2).thresholdError;
  nPotPeaks = eva(2).parameters(5);
  percentilePotX = eva(2).parameters(6);
  % 72 is the minumum interval in time steps used by procedure
  % tsGetPOTAndRlargest, when it calls findpeaks.
  dtPeaks = minPeakDistance/2;
  dtPotX = (timeStamps(end) - timeStamps(1))/length(series)*dtPeaks;

  epsilonPotNS = epsilonPotX;
  errEpsilonPotNS = errEpsilonPotX;
  sigmaPotNS = sigmaPotX*trasfData.stdDevSeries;
  % propagating the errors on stdDevSeries and sigmaPotX to sigmaPotNs.
  % err(sigmaNs) = sqrt{ [sigmaX*err(stdDev)]^2 + [stdDev*err(sigmaX)]^2 }
  % the error on sigmaGevNs is time dependant.
  errSigmaPotFit = trasfData.stdDevSeries.*errSigmaPotX;
  errSigmaPotTransf = sigmaPotX*trasfData.stdDevError;
  errSigmaPotNS = (  errSigmaPotTransf.^2   +  errSigmaPotFit.^2  ).^.5;
  thresholdPotNS = thresholdPotX*trasfData.stdDevSeries + trasfData.trendSeries;
  % propagating the errors on stdDevSeries and trendSeries to thresholdPotNs.
  % err(thresholdPotNs) = sqrt{ [thresholdPotX*err(stdDev)]^2 + err(trend)^2 }
  % the error on thresholdPotNs is constant.
  thresholdErrFit = 0;
  thresholdErrTransf = ( (trasfData.stdDevSeries*errThresholdPotX).^2 + (thresholdPotX*trasfData.stdDevError).^2  +  trasfData.trendError^2  ).^.5;
  thresholdErr = thresholdErrTransf;
  potParams.epsilon = epsilonPotNS;
  potParams.sigma = sigmaPotNS;
  potParams.threshold = thresholdPotNS;
  potParams.percentile = percentilePotX;
  potParams.timeDelta = dtPotX;
  potParams.timeDeltaYears = dtPotX/365.2425;
  potParams.timeHorizonStart = min(trasfData.timeStamps);
  potParams.timeHorizonEnd = max(trasfData.timeStamps);
  potParams.nPeaks = nPotPeaks;

  potParamStdErr.epsilonErr = errEpsilonPotNS;

  potParamStdErr.sigmaErrFit = errSigmaPotFit;
  potParamStdErr.sigmaErrTransf = errSigmaPotTransf;
  potParamStdErr.sigmaErr = errSigmaPotNS;

  potParamStdErr.thresholdErrFit = thresholdErrFit;
  potParamStdErr.thresholdErrTransf = thresholdErrTransf;
  potParamStdErr.thresholdErr = thresholdErr;

  potObj.method = eva(2).method;
  potObj.parameters = potParams;
  potObj.paramErr = potParamStdErr;
  potObj.stationaryParams = eva(2);
  potObj.objs = [];
else
  potObj.method = eva(2).method;
  potObj.parameters = [];
  potObj.paramErr = [];
  potObj.stationaryParams = [];
  potObj.objs = [];
end

%% setting output objects
clear nonStationaryEvaParams;
nonStationaryEvaParams(1) = gevObj;
nonStationaryEvaParams(2) = potObj;

stationaryTransformData = trasfData;

end

