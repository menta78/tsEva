function [monthlyMax, monthlyMaxDate, monthlyMaxIndx] = tsEvaComputeMonthlyMaxima(timeAndSeries)

timeStamps = timeAndSeries(:,1);
srs = timeAndSeries(:,2);
lsrs = length(srs);
indexes = (1:1:lsrs)';

dtvec = datevec(timeStamps);
dtvecMonth = dtvec(:,1:2);
dtvecMonth = cat(2, dtvecMonth, ones(length(dtvec), 1));
tsMonth = datenum(dtvecMonth);
tsMonthUnique = unique(tsMonth);

dataByMonth = cellfun(@(tsmonth) srs(tsMonth == tsmonth), num2cell(tsMonthUnique), 'uniformoutput', false);
indxByMonth = cellfun(@(tsmonth) indexes(tsMonth == tsmonth), num2cell(tsMonthUnique), 'uniformoutput', false);
monthlyMax = cellfun(@(mnthData) max(mnthData), dataByMonth);
monthlyMaxIndx0 = cellfun(@(imonth) find(dataByMonth{imonth} == monthlyMax(imonth)), num2cell(1:1:length(monthlyMax)), 'uniformoutput', false);
monthlyMaxIndx0(cellfun(@(c) isempty(c), monthlyMaxIndx0)) = {1};
monthlyMaxIndx0 = [monthlyMaxIndx0{:}];
monthlyMaxIndx = cellfun(@(imonth) indxByMonth{imonth}(monthlyMaxIndx0(imonth)), num2cell(1:1:length(monthlyMax)));

monthlyMaxDate = timeStamps(monthlyMaxIndx);

monthlyMax = monthlyMax';
monthlyMaxDate = monthlyMaxDate';
monthlyMaxIndx = monthlyMaxIndx';

