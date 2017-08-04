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

  try
    myRetPeriod_ = tsInterp1Extrap(retLevel, retPeriod, myLevel(myLevNonNan), logExtrap)';
  catch
    disp('    ... log extrapolation of return failed for some levels. Setting it to infinite.');
    myRetPeriod_ = interp1(retLevel, retPeriod, myLevel(myLevNonNan));
    myRetPeriod_(isnan(myRetPeriod_)) = inf;
  end
  
  try
    myRetPeriodCISup_ = tsInterp1Extrap(rlLow, retPeriod, myLevel(myLevNonNan), logExtrap)';
  catch
    disp('    ... log extrapolation of return period superior bar failed for some levels. Setting there an infinite superior ci.');
    myRetPeriodCISup_ = interp1(rlLow, retPeriod, myLevel(myLevNonNan));
    myRetPeriodCISup_(isnan(myRetPeriodCISup_)) = inf;
  end
  
  try
    myRetPeriodCIInf_ = tsInterp1Extrap(rlHigh, retPeriod, myLevel(myLevNonNan), logExtrap)';
  catch
    disp('    ... log extrapolation of return period superior bar failed for some levels. Setting there an infinite superior ci.');
    myRetPeriodCIInf_ = interp1(rlHigh, retPeriod, myLevel(myLevNonNan));
    myRetPeriodCIInf_(isnan(myRetPeriodCIInf_)) = inf;
  end
  
  myRetPeriod = ones(size(myLevel))*nan;
  myRetPeriod(myLevNonNan) = myRetPeriod_;
  
  myRetPeriodCISup = ones(size(myLevel))*nan;
  myRetPeriodCISup(myLevNonNan) = myRetPeriodCISup_;
  
  myRetPeriodCIInf = ones(size(myLevel))*nan;
  myRetPeriodCIInf(myLevNonNan) = myRetPeriodCIInf_;

end

