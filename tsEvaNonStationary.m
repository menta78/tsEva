function [nonStationaryEvaParams, stationaryTransformData] = tsEvaNonStationary( timeAndSeries, timeWindow, varargin )
%% sample calls
% nonStatEvaParams = nonStationaryEvaJRCApproach(ms, timeWindow, 'pcts',[95], 'minPeakDistance', 72)
%     samples POT data using a fixed 95 percentile threshold, with peaks at
%     a minimum distance of 
% nonStatEvaParams = nonStationaryEvaJRCApproach(ms, timeWindow, 'desiredeventsperyear', 6)
%     samples POT data looking for a threshold so that we have an average
%     of 6 events every year.
% nonStatEvaParams = nonStationaryEvaJRCApproach(ms, timeWindow)
%     samples POT data looking for a threshold so that we have an average
%     of 5 events every year.
%% %%%%%%%%%%%%%

args.transfType = 'trend';
args.minPeakDistanceInDays = -1;
args = tsEasyParseNamedArgs(varargin, args);
minPeakDistanceInDays = args.minPeakDistanceInDays;
transfType = args.transfType;
if ~( strcmpi(transfType, 'trend') || strcmpi(transfType, 'seasonal') || strcmpi(transfType, 'seasonalAdditive') )
    error('nonStationaryEvaJRCApproach: transfType can be in (trend, seasonal, seasonalAdditive)');
end
if strcmpi(transfType, 'seasonalAdditive')
    error('nonStationaryEvaJRCApproach: transfType==seasonalAdditive is not yet supported')
end
if minPeakDistanceInDays == -1
    error('label parameter ''minPeakDistanceInDays'' must be set')
end
    

timeStamps = timeAndSeries(:, 1);
series = timeAndSeries(:, 2);

if strcmpi(transfType, 'trend')
    trasfData = tsEvaTransformSeriesToStationaryTrendOnly( timeStamps, series, timeWindow );
    gevMaxima = 'annual';
    potEventsPerYear = 5;
elseif strcmpi(transfType, 'seasonal')
    trasfData = tsEvaTransformSeriesToStationaryMultiplicativeSeasonality( timeStamps, series, timeWindow );
    gevMaxima = 'monthly';
    potEventsPerYear = 12;
end
ms = cat(2, trasfData.timeStamps, trasfData.stationarySeries);
%dt = trasfData.timeStamps(2) - trasfData.timeStamps(1);
dt = tsEvaGetTimeStep(trasfData.timeStamps);
minPeakDistance = minPeakDistanceInDays/dt;

%% estimating the non stationary GEV parameters
fprintf('\n');
disp('Executing the stationary eva')
pointData = tsEvaSampleData(ms, 'meanEventsPerYear', potEventsPerYear, varargin{:});
evaAlphaCI = .68; % in a gaussian approximation alphaCI~68% corresponds to 1 sigma confidence
[~, eva] = tsEVstatistics(pointData, 'alphaci', evaAlphaCI, 'gevmaxima', gevMaxima);
eva(2).thresholdError = pointData.POT.thresholdError;
fprintf('\n');

% !!! Assuming a Gaussian approximation to compute the standard errors for
% the GEV parameters
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
errSigmaGevNS = (  (sigmaGevX*trasfData.stdDevError).^2   +  (trasfData.stdDevSeries.*errSigmaGevX).^2  ).^.5;
muGevNS = trasfData.stdDevSeries*muGevX + trasfData.trendSeries;
% propagating the errors on stdDevSeries, trendSeries and sigmaGevX to muGevNS.
% err(muNs) = sqrt{ [muX*err(stdDev)]^2 + [stdDev*err(muX)]^2 + err(trend)^2 }
% the error on muGevNS is time dependant.
errMuGevNS = (  (muGevX*trasfData.stdDevError).^2   +  (trasfData.stdDevSeries.*errMuGevX).^2   +  trasfData.trendError^2  ).^.5;
gevParams.epsilon = epsilonGevNS;
gevParams.sigma = sigmaGevNS;
gevParams.mu = muGevNS;
gevParamStdErr.epsilonErr = errEpsilonGevNS;
gevParamStdErr.sigmaErr = errSigmaGevNS;
gevParamStdErr.muErr = errMuGevNS;
gevObj.method = eva(1).method;
gevObj.parameters = gevParams;
gevObj.paramErr = gevParamStdErr;
gevObj.stationaryParams = eva(1);
gevObj.objs.monthlyMaxIndexes = pointData.monthlyMaxIndexes;

%% estimating the non stationary GPD parameters
% !!! Assuming a Gaussian approximation to compute the standard errors for
% the GEV parameters
epsilonPotX = eva(2).parameters(2);
errEpsilonPotX = epsilonPotX - eva(2).paramCIs(1, 2);
sigmaPotX = eva(2).parameters(1);
errSigmaPotX = sigmaPotX - eva(2).paramCIs(1, 1);
thresholdPotX = eva(2).parameters(3);
errThresholdPotX = eva(2).thresholdError;
nnSrs = trasfData.stationarySeries(~isnan(trasfData.stationarySeries));
thresholdExceedProbability = sum(nnSrs > thresholdPotX)/length(nnSrs); % should be 1 - percentile
percentilePotX = eva(2).parameters(6);
% 72 is the minumum interval in time steps used by procedure
% tsGetPOTAndRlargest, when it calls findpeaks.
% If it finds 2 peaks at a distance of 72 hours, it means that there is
% something in the middle
dtPeaks = minPeakDistance/2;
dtPotX = (timeStamps(end) - timeStamps(1))/length(series)*dtPeaks;

epsilonPotNS = epsilonPotX;
errEpsilonPotNS = errEpsilonPotX;
sigmaPotNS = sigmaPotX*trasfData.stdDevSeries;
% propagating the errors on stdDevSeries and sigmaPotX to sigmaPotNs.
% err(sigmaNs) = sqrt{ [sigmaX*err(stdDev)]^2 + [stdDev*err(sigmaX)]^2 }
% the error on sigmaGevNs is time dependant.
errSigmaPotNS = (  (sigmaPotX*trasfData.stdDevError).^2   +  (trasfData.stdDevSeries.*errSigmaPotX).^2  ).^.5;
thresholdPotNS = thresholdPotX*trasfData.stdDevSeries + trasfData.trendSeries;
% propagating the errors on stdDevSeries and trendSeries to thresholdPotNs.
% err(thresholdPotNs) = sqrt{ [thresholdPotX*err(stdDev)]^2 + err(trend)^2 }
% the error on thresholdPotNs is constant.
thresholdErr = ( (trasfData.stdDevSeries*errThresholdPotX).^2 + (thresholdPotX*trasfData.stdDevError).^2  +  trasfData.trendError^2  ).^.5;
potParams.epsilon = epsilonPotNS;
potParams.sigma = sigmaPotNS;
potParams.threshold = thresholdPotNS;
potParams.percentile = percentilePotX;
potParams.timeDelta = dtPotX;
potParams.timeDeltaYears = dtPotX/365.2425;
potParamStdErr.epsilonErr = errEpsilonPotNS;
potParamStdErr.sigmaErr = errSigmaPotNS;
potParamStdErr.thresholdErr = thresholdErr;
potObj.method = eva(2).method;
potObj.parameters = potParams;
potObj.paramErr = potParamStdErr;
potObj.stationaryParams = eva(2);
potObj.objs = [];

%% setting output objects
nonStationaryEvaParams(1) = gevObj;
nonStationaryEvaParams(2) = potObj;

stationaryTransformData = trasfData;

end

