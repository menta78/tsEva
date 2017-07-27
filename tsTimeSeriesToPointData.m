function pointData = tsTimeSeriesToPointData( ms, potThreshold, potThresholdError )

% tsTimeSeriesToPointData: given a ms produces a structure pointData
% like the one produced by tsEvaSampleData

  ms1 = ms(ms(:,2) > potThreshold, :);
  percentile = (1 - size(ms1, 1)/size(ms, 1))*100;

  pointData.completeSeries = ms1;
  
  POTdata.threshold = potThreshold;
  POTdata.thresholdError = potThresholdError;
  POTdata.percentile = percentile;
  POTdata.peaks = ms1(:,2);
  POTdata.ipeaks=1:size(ms1, 1);
  POTdata.sdpeaks=ms1(:,1);  
  pointData.POT = POTdata;
  
  [pointData.annualMax, pointData.annualMaxTimeStamp, pointData.annualMaxIndexes] = tsEvaComputeAnnualMaxima(ms1);
  [pointData.monthlyMax, pointData.monthlyMaxTimeStamp, pointData.monthlyMaxIndexes] = tsEvaComputeMonthlyMaxima(ms1);
  
  yrs = unique(tsYear(ms1(:,1)));
  yrs = yrs - min(yrs);
  pointData.years = (nanmin(yrs):1:nanmax(yrs))';
  
  pointData.Percentiles = [0];



end

