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
  
  myLevNonNan = ~isnan(myLevel);
  myRetPeriod_ = tsInterp1Extrap(retLevel, retPeriod, myLevel(myLevNonNan), logExtrap)';
  myRetPeriodCISup_ = tsInterp1Extrap(rlLow, retPeriod, myLevel(myLevNonNan), logExtrap)';
  myRetPeriodCIInf_ = tsInterp1Extrap(rlHigh, retPeriod, myLevel(myLevNonNan), logExtrap)';
  
  myRetPeriod = ones(size(myLevel))*nan;
  myRetPeriod(myLevNonNan) = myRetPeriod_;
  
  myRetPeriodCISup = ones(size(myLevel))*nan;
  myRetPeriodCISup(myLevNonNan) = myRetPeriodCISup_;
  
  myRetPeriodCIInf = ones(size(myLevel))*nan;
  myRetPeriodCIInf(myLevNonNan) = myRetPeriodCIInf_;

end

