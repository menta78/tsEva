function [ myRetPeriod, myRetPeriodCISup, myRetPeriodCIInf ] = tsGetReturnPeriodOfLevel( retPeriod, retLevel, retLevError, myLevel, varargin )
  % given a list of retPeriod and corresponding retLevel with error,
  % estimates the return period for a level myLevel, and the related confidence interval.

  args.lowerBoundTo0 = true;
  args.logExtrap = true;
  args.cpP = .68; % by default estimating 1-sigma confidence interval
  args = tsEasyParseNamedArgs(varargin, args);
  logExtrap = args.logExtrap;
  cpP = args.cpP;
  lowerBoundTo0 = args.lowerBoundTo0;

  [rlLow, rlHigh] = tsEstimateConfidenceIntervalOfRL(retLevel, retLevError, cpP);
  
  if lowerBoundTo0
    retPeriod = [0; retPeriod];
    retLevel = [0; retLevel];
    rlLow = [0; rlLow];
    rlHigh = [0; rlHigh];
  end
  
  myRetPeriod = tsInterp1Extrap(retLevel, retPeriod, myLevel, logExtrap)';
  myRetPeriodCISup = tsInterp1Extrap(rlLow, retPeriod, myLevel, logExtrap)';
  myRetPeriodCIInf = tsInterp1Extrap(rlHigh, retPeriod, myLevel, logExtrap)';

end

