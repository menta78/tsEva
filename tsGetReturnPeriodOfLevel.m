function [ myRetPeriod, myRetPeriodError, myRetPeriodErrorSup, myRetPeriodErrorInf ] = tsGetReturnPeriodOfLevel( retPeriod, retLevel, retLevError, myLevel, varargin )
% given a list of retPeriod and corresponding retLevel with error,
% estimates the return period for a level myLevel, and the related error.

args.lowerBoundTo0 = true;
%args.interpMethod = 'linear';
args.logExtrap = true;
args = tsEasyParseNamedArgs(varargin, args);
lowerBoundTo0 = args.lowerBoundTo0;
%interpMethod = args.interpMethod;
logExtrap = args.logExtrap;

bndTo0 = lowerBoundTo0 && retPeriod(1);
if bndTo0 ~= 0
  retPeriod_ = [0; sort(retPeriod)];
  retLevel_ = [0; sort(retLevel)];
  retLevError_ = [0; sort(retLevError)];
  rlLow = retLevel_ - retLevError_;
  rlHigh = retLevel_ + retLevError_;
else
  retPeriod_ = retPeriod;
  retLevel_ = retLevel;
  rlLow = retLevel_ - retLevError;
  rlHigh = retLevel_ + retLevError;
end
ll = length(rlLow);
rlLowDiff = [1; diff(rlLow)];
rlLowLastIncIndx = sum(rlLowDiff > 0);
if rlLowLastIncIndx <= ll
  rlLow(rlLowLastIncIndx:ll) = retLevel_(rlLowLastIncIndx:ll) - retLevError_(rlLowLastIncIndx);
end

% myRetPeriod = interp1(retLevel_, retPeriod_, myLevel, interpMethod, 'extrap');
% myRetPeriodSup = interp1(rlLow, retPeriod_, myLevel, interpMethod, 'extrap');
% myRetPeriodInf = interp1(rlHigh, retPeriod_, myLevel, interpMethod, 'extrap');
% myRetPeriodError = (myRetPeriodSup - myRetPeriodInf)/2.;

myRetPeriod = tsInterp1Extrap(retLevel_, retPeriod_, myLevel, logExtrap);
myRetPeriodSup = tsInterp1Extrap(rlLow, retPeriod_, myLevel, logExtrap);
myRetPeriodInf = tsInterp1Extrap(rlHigh, retPeriod_, myLevel, logExtrap);
myRetPeriodError = (myRetPeriodSup - myRetPeriodInf)/2.;

myRetPeriodErrorSup = myRetPeriodSup - myRetPeriod;
myRetPeriodErrorInf = myRetPeriod - myRetPeriodInf;

end

