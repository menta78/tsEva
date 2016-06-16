function [ filledTimeStamps, filledSeries, dt ] = tsEvaFillSeries( timeStamps, series )
  indxs = ~isnan(series);
  timeStamps = timeStamps(indxs);
  series = series(indxs);

  mint = min(timeStamps);
  maxt = max(timeStamps);
  dt = min(diff(timeStamps));
  filledTimeStamps = (mint:dt:maxt)';
  filledSeries = interp1(timeStamps, series, filledTimeStamps, 'nearest');
  filledSeries = tsRemoveConstantSubseries(filledSeries, 4);
end

