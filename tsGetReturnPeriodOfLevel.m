function [ myRetPeriod, myRetPeriodCISup, myRetPeriodCIInf ] = tsGetReturnPeriodOfLevel( retPeriod, retLevel, retLevError, myLevel, varargin )
  % given a list of retPeriod and corresponding retLevel with error,
  % estimates the return period for a level myLevel, and the related confidence interval.

  args.logExtrap = true;
  args.cpP = .68; % by default estimating 1-sigma confidence interval
  args = tsEasyParseNamedArgs(varargin, args);
  logExtrap = args.logExtrap;
  cpP = args.cpP;

  [rlLow, rlHigh] = tsEstimateConfidenceIntervalOfRL(retLevel, retLevError, cpP);
  
  myLevNonNan = ~isnan(myLevel) & ~isinf(myLevel);

  if logExtrap
    myRetPeriod_ = tsInterp1Extrap(retLevel, retPeriod, myLevel(myLevNonNan), logExtrap)';
    try
      myRetPeriodCISup_ = tsInterp1Extrap(rlLow, retPeriod, myLevel(myLevNonNan), logExtrap)';
    catch
      myRetPeriodCISup_ = interp1(rlLow, retPeriod, myLevel(myLevNonNan));
      myRetPeriodCISup_(isnan(myRetPeriodCISup_)) = inf;
    end
    try
      myRetPeriodCIInf_ = tsInterp1Extrap(rlHigh, retPeriod, myLevel(myLevNonNan), logExtrap)';
    catch
      myRetPeriodCIInf_ = interp1(rlHigh, retPeriod, myLevel(myLevNonNan));
    end
  else
    myRetPeriod_ = interp1(retLevel, retPeriod, myLevel(myLevNonNan));
    myRetPeriodCISup_ = interp1(rlLow, retPeriod, myLevel(myLevNonNan));
    myRetPeriodCIInf_ = interp1(rlHigh, retPeriod, myLevel(myLevNonNan));
  end    
  
  myRetPeriod = ones(size(myLevel))*nan;
  myRetPeriod(myLevNonNan) = myRetPeriod_;
  
  myRetPeriodCISup = ones(size(myLevel))*nan;
  myRetPeriodCISup(myLevNonNan) = myRetPeriodCISup_;
  
  myRetPeriodCIInf = ones(size(myLevel))*nan;
  myRetPeriodCIInf(myLevNonNan) = myRetPeriodCIInf_;

end

