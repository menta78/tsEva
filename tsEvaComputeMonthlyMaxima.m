function [monthlyMax, monthlyMaxDate, monthlyMaxIndx] = tsEvaComputeMonthlyMaxima(timeAndSeries)

  timeStamps = timeAndSeries(:,1);
  srs = timeAndSeries(:,2);

  tmvec = datevec(timeStamps);
  yrs = tmvec(:,1);
  mnts = tmvec(:,2);
  mnttmvec = [yrs mnts];
  valsIndxs = 1:numel(srs);

  function maxIndx = findMax(indxs)
    [~, maxSubIndx] = max(srs(indxs));
    maxIndx = indxs(maxSubIndx);
  end

  monthlyMaxIndx = accumarray(mnttmvec, valsIndxs, [], @findMax);
  monthlyMaxIndx = sort(monthlyMaxIndx(monthlyMaxIndx ~= 0));
  monthlyMax = srs(monthlyMaxIndx);
  monthlyMaxDate = timeStamps(monthlyMaxIndx);

end

