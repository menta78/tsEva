function [stationaryEvaParams, isValid] = tsEvaStationary( timeAndSeries, varargin )
% tsEvaStationary:
% executes a regular stationary EVA on timeAndSeries.
% stationaryEvaParams includes the parameters estimated for GEV and GPD
%% %%%%%%%%%%%%%

args.minPeakDistanceInDays = -1;
args.potEventsPerYear = 5;
args.gevMaxima = 'annual';
args.gevType = 'GEV'; % can be 'GEV' or 'Gumbel'
args.doSampleData = true;
args.potThreshold = nan;
args.evdType = {'GEV', 'GPD'};
args = tsEasyParseNamedArgs(varargin, args);
minPeakDistanceInDays = args.minPeakDistanceInDays;
if minPeakDistanceInDays == -1
    error('label parameter ''minPeakDistanceInDays'' must be set');
end
potEventsPerYear = args.potEventsPerYear;
gevMaxima = args.gevMaxima;
gevType = args.gevType;
doSampleData = args.doSampleData;
evdType = args.evdType;
potThreshold = args.potThreshold;
   
stationaryEvaParams = [];

fprintf('\n');
disp('Executing stationary eva');
if doSampleData
  % removing nan
  cnd = ~isnan(timeAndSeries(:,2));
  timeAndSeries = timeAndSeries(cnd, :);
  pointData = tsEvaSampleData(timeAndSeries, 'meanEventsPerYear', potEventsPerYear, varargin{:});
else
  if isnan(potThreshold)
    error('if doSampleData==false, than the POT is not performed, and you need to provide a value for the threshold.');
  end
  pointData = tsTimeSeriesToPointData(timeAndSeries, potThreshold, 0);
end
evaAlphaCI = .68; % in a gaussian approximation alphaCI~68% corresponds to 1 sigma confidence
[~, eva, isValid] = tsEVstatistics(pointData, 'alphaci', evaAlphaCI, ...
       'gevmaxima', gevMaxima, 'gevType', gevType, 'evdType', evdType, varargin{:});
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

  gevParams.epsilon = epsilonGevX;
  gevParams.sigma = sigmaGevX;
  gevParams.mu = muGevX;
  if strcmpi(gevMaxima, 'annual')
    gevParams.timeDelta = 365.25;
    gevParams.timeDeltaYears = 1;
  elseif strcmpi(gevMaxima, 'monthly')
    gevParams.timeDelta = 365.25/12.;
    gevParams.timeDeltaYears = 1/12.;
  end
  
  gevParamStdErr.epsilonErr = errEpsilonX;
  gevParamStdErr.sigmaErr = errSigmaGevX;
  gevParamStdErr.muErr = errMuGevX;
  
  gevObj.method = eva(1).method;
  gevObj.parameters = gevParams;
  gevObj.paramErr = gevParamStdErr;
  gevObj.objs.monthlyMaxIndexes = pointData.monthlyMaxIndexes;
else
  gevObj.method = eva(1).method;
  gevObj.parameters = [];
  gevObj.paramErr = [];
  gevObj.objs.monthlyMaxIndexes = [];
end

%% estimating the non stationary GPD parameters
% !!! Assuming a Gaussian approximation to compute the standard errors for
% the GEV parameters
epsilonPotX = eva(2).parameters(2);
errEpsilonPotX = epsilonPotX - eva(2).paramCIs(1, 2);
sigmaPotX = eva(2).parameters(1);
errSigmaPotX = sigmaPotX - eva(2).paramCIs(1, 1);
thresholdPotX = eva(2).parameters(3);
errThresholdPotX = eva(2).thresholdError;
percentilePotX = eva(2).parameters(6);
nPotPeaks = eva(2).parameters(5);

timeStamps = timeAndSeries(:,1);
dt = tsEvaGetTimeStep(timeStamps);
minPeakDistance = minPeakDistanceInDays/dt;
% 72 is the minumum interval in time steps used by procedure
% tsGetPOTAndRlargest, when it calls findpeaks.
% If it finds 2 peaks at a distance of 72 hours, it means that dtPeaks = 72
% hours
dtPeaks = minPeakDistance;
dtSample =  (timeStamps(end) - timeStamps(1))/length(timeStamps);
dtPotX = max(dtSample, dtSample*dtPeaks);

potParams.epsilon = epsilonPotX;
potParams.sigma = sigmaPotX;
potParams.threshold = thresholdPotX;
potParams.percentile = percentilePotX;
potParams.timeDelta = dtPotX;
potParams.timeDeltaYears = dtPotX/365.2425;
potParams.timeHorizonStart = min(timeStamps);
potParams.timeHorizonEnd = max(timeStamps);
potParams.nPeaks = nPotPeaks;

potParamStdErr.epsilonErr = errEpsilonPotX;
potParamStdErr.sigmaErr = errSigmaPotX;
potParamStdErr.thresholdErr = errThresholdPotX;

potObj.method = eva(2).method;
potObj.parameters = potParams;
potObj.paramErr = potParamStdErr;
potObj.objs = pointData.POT.ipeaks;

%% setting output objects
clear stationaryEvaParams;
stationaryEvaParams(1) = gevObj;
stationaryEvaParams(2) = potObj;

end

