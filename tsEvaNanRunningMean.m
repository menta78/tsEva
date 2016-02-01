function [ rnmn ] = tsEvaNanRunningMean( series, windowSize )

minNThreshold = 1;

rnmn = zeros([length(series), 1])*nan;
dx = ceil(windowSize/2);
l = length(series);
sm = 0;
n = 0;
for ii = 1:l
    minindx = max(ii - dx, 1);
    maxindx = min(ii + dx, l);
    if (ii == 1)
        subsrs = series(minindx:maxindx);
        sm = nansum(subsrs);
        n = sum(~isnan(subsrs));
    else
        if minindx > 1
            sprev = series(minindx - 1);
            if ~isnan(sprev)
                sm = sm - sprev;
                n = n - 1;
            end
        end
        if maxindx < l
            snext = series(maxindx + 1);
            if ~isnan(snext)
                sm = sm + snext;
                n = n + 1;
            end
        end
    end
      
    if n > minNThreshold
        rnmn(ii) = sm/n;
    else
        rnmn(ii) = nan;
    end
end

end

