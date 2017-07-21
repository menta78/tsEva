function pointData = tsTimeSeriesToPointData( ms, potThreshold, potThresholdError )

% tsTimeSeriesToPointData: given a ms produces a structure pointData
% like the one produced by tsEvaSampleData

  pointData.completeSeries = ms;
  
  POTdata.threshold = potThreshold;
  POTdata.thresholdError = 0;
  POTdata.percentile = 0;
  POTdata.peaks = ms(:,2);
  POTdata.ipeaks=1:size(ms, 1);
  POTdata.sdpeaks=ms(:,1);  
  pointData.POT = POTdata;
  
  [pointData.annualMax, pointData.annualMaxTimeStamp, pointData.annualMaxIndexes] = tsEvaComputeAnnualMaxima(ms);
  [pointData.monthlyMax, pointData.monthlyMaxTimeStamp, pointData.monthlyMaxIndexes] = tsEvaComputeMonthlyMaxima(ms);
  
  yrs = unique(tsYear(ms(:,1)));
  yrs = yrs - min(yrs);
  pointData.years = (nanmin(yrs):1:nanmax(yrs))';
  
  pointData.Percentiles = [0];



end

