function pointData = tsEvaSampleData(ms, varargin)
  pctsDesired=[90 95 99 99.9];

  args.meanEventsPerYear = 5;
  %args.pcts = [50 70 85:2:99 99.1:0.1:99.9 99.91:0.01:99.99];
  args.potPercentiles = [50 70 85:2:97];
  args = tsEasyParseNamedArgs(varargin, args);
  meanEventsPerYear = args.meanEventsPerYear;
  potPercentiles = args.potPercentiles;

  [POTData,rLargestData]=tsGetPOTAndRlargest(ms, potPercentiles, meanEventsPerYear, varargin{:});

  vals=prctile(ms(:,2), pctsDesired);
  percentiles.precentiles=pctsDesired;
  percentiles.values=vals;
  
  pointData.completeSeries = ms;
  pointData.POT = POTData;
  pointData.Rlargest = rLargestData;
  if ~isempty(rLargestData)
    pointData.annualMax = rLargestData.peaks(:,1)';
    pointData.annualMaxTimeStamp = rLargestData.sdpeaks(:,1)';
  else
    pointData.annualMax = [];
    pointData.annualMaxTimeStamp = [];
  end
  [pointData.monthlyMax, pointData.monthlyMaxTimeStamp, pointData.monthlyMaxIndexes] = tsEvaComputeMonthlyMaxima(ms);
  
  yrs = unique(tsYear(ms(:,1)));
  yrs = yrs - min(yrs);
  pointData.years = (nanmin(yrs):1:nanmax(yrs))';
  
  pointData.Percentiles = percentiles;
end

