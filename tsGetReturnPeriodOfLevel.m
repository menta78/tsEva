function [ myRetPeriod, myRetPeriodError ] = tsGetReturnPeriodOfLevel( retPeriod, retLevel, retLevError, myLevel, varargin )
% given a list of retPeriod and corresponding retLevel with error,
% estimates the return period for a level myLevel, and the related error.

args.lowerBoundTo0 = true;
args.interpMethod = 'spline';
args = tsEasyParseNamedArgs(varargin, args);
lowerBoundTo0 = args.lowerBoundTo0;
interpMethod = args.interpMethod;

bndTo0 = lowerBoundTo0 && retPeriod(1);
if bndTo0 ~= 0
  retPeriod_ = [0; sort(retPeriod)];
  retLevel_ = [0; sort(retLevel)];
  retLevError_ = [0; sort(retLevError)];
  rlLow = retLevel_ - retLevError_;
  rlLow(rlLow < 0) = 0;
  rlHigh = retLevel_ + retLevError_;
else
  retPeriod_ = retPeriod;
  retLevel_ = retLevel;
  rlLow = retLevel_ - retLevError;
  rlHigh = retLevel_ + retLevError;
end

myRetPeriod = interp1(retLevel_, retPeriod_, myLevel, interpMethod, 'extrap');
myRetPeriodSup = interp1(rlLow, retPeriod_, myLevel, interpMethod, 'extrap');
myRetPeriodInf = interp1(rlHigh, retPeriod_, myLevel, interpMethod, 'extrap');
myRetPeriodError = (myRetPeriodSup - myRetPeriodInf)/2.;

end

