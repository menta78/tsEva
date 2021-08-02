function [annualMax, annualMaxDate, annualMaxIndx] = tsEvaComputeAnnualMaximaMtx(timeStamps, srs)

if size(timeStamps, 1) ~= size(srs, 1)
  error('tsEvaComputeAnnualMtxMaxima: the 1st size of srs should equal to the size of timeStamps');
end

tmvec = datevec(double(timeStamps));
years = tmvec(:,1);

yu = unique(years);
nyrs = length(yu);

srsSize = size(srs);
amxSize = [nyrs, srsSize(2:end)];
annualMax = ones(amxSize)*nan;
annualMaxIndx = ones(amxSize)*nan;
lastMxYrIndx = 0;
for iy = 1:nyrs
  yr = yu(iy);
  indx = years == yr;
  [ymx, indxmx] = max(srs(indx, :, :), [], 1);
  annualMax(iy, :, :) = squeeze(ymx);
  annualMaxIndx(iy, :, :) = indxmx + lastMxYrIndx;
  lastMxYrIndx = lastMxYrIndx + sum(indx);
end

if nargout > 1
  tmstmpmtx = timeStamps(:, ones(amxSize(2), 1), ones(amxSize(3), 1));
  annualMaxDate = tmstmpmtx(annualMaxIndx);
end

end