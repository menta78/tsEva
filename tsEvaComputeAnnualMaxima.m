function [annualMax, annualMaxDate, annualMaxIndx] = tsEvaComputeAnnualMaxima(timeAndSeries)

timeStamps = timeAndSeries(:,1);
srs = timeAndSeries(:,2);

tmvec = datevec(timeStamps);
years = tmvec(:,1);
srsIndices = 1:numel(srs);

function maxIndx = findMax(subIndxs)
  [~, subIndxMaxIndx] = max(srs(subIndxs));
  maxIndx = subIndxs(subIndxMaxIndx);
end

annualMaxIndx = accumarray(years, srsIndices, [], @findMax);
annualMaxIndx = annualMaxIndx(annualMaxIndx ~= 0);
annualMax = srs(annualMaxIndx);
annualMaxDate = timeStamps(annualMaxIndx);

end