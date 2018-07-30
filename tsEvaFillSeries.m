function [ filledTimeStamps, filledSeries, dt ] = tsEvaFillSeries( timeStamps, series )
  indxs = ~isnan(series);
  timeStamps = timeStamps(indxs);
  series = series(indxs);
  
  %   ensuring monotonic time vector and no dubplicates. In case of dubplicate time stamps the highest value is considered
  [newTs,~,idx] = unique(sort(timeStamps));
  newSeries = accumarray(idx,series,[],@max);
  
  mint = min(newTs);
  maxt = max(newTs);
  dt = min(diff(newTs));
  if dt >= 350 && dt <= 370
    % this is an annual series
    mindtVec = datevec(mint);
    mindtY = mindtVec(1);
    maxdtVec = datevec(maxt);
    maxdtY = maxdtVec(1);
    years = (mindtY:maxdtY)';
    dtvec = [years, ones(size(years)), ones(size(years))];
    filledTimeStamps = datenum(dtvec);
  elseif dt >= 28 && dt <= 31
    % this is a monthly series
    mindtVec = datevec(mint);
    mindtY = mindtVec(1);
    maxdtVec = datevec(maxt);
    maxdtY = maxdtVec(1);
    years = (mindtY:maxdtY);
    months = 1:12;
    [ymtx, mmtx] = meshgrid(years, months);
    ys = ymtx(:);
    ms = mmtx(:);
    dtvec = [ys, ms, ones(size(ys))];
    filledTimeStamps = datenum(dtvec);
  else  
    filledTimeStamps = (mint:dt:maxt)';
  end
  

  
  filledSeries = interp1(newTs, newSeries, filledTimeStamps, 'nearest');
  filledSeries = tsRemoveConstantSubseries(filledSeries, 4);
end

