function [ rnmn ] = tsEvaNanRunningVariance( series, windowSize )
%!!! series must be 0 averaged!!

minNThreshold = 1;

rnmn = zeros([length(series), 1])*nan;
dx = ceil(windowSize/2);
l = length(series);
smsq = 0;
n = 0;
for ii = 1:l
    minindx = max(ii - dx, 1);
    maxindx = min(ii + dx, l);
    if (ii == 1)
        subSqSrs = series(minindx:maxindx).^2;
        smsq = nansum(subSqSrs);
        n = sum(~isnan(subSqSrs));
    else
        if minindx > 1
            sprev = series(minindx - 1);
            if ~isnan(sprev)
                smsq = max(0, smsq - sprev^2);
                n = n - 1;
            end
        end
        if maxindx < l
            snext = series(maxindx + 1);
            if ~isnan(snext)
                smsq = smsq + snext^2;
                n = n + 1;
            end
        end
    end
      
    if n > minNThreshold
        rnmn(ii) = smsq/n;
    else
        rnmn(ii) = nan;
    end
end

end

