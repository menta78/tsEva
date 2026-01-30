function [ rnmn, rnvar, rn3mom, rn4mom ] = tsEvaNanRunningStatistics( series, windowSize )
% tsEvaNanRunningStatistics: returns the the moving statistical momentums
% to the forth.
% rnmn: running mean
% rnvar: running variance
% rn3mom: running third statistical momentum
% rn4mom: running fourth statistical momentum

% at both extremeties of the series, half windowSize is used which
% gradually increases to reach windowSize; once windowSize is reached,
% windowSize is rolled throghout the series

minNThreshold = 1;

rnmn = tsEvaNanRunningMean(series, windowSize);
rnvar = zeros([length(series), 1])*nan;
rn3mom = zeros([length(series), 1])*nan;
rn4mom = zeros([length(series), 1])*nan;

dx = ceil(windowSize/2);
l = length(series);
sm = 0;
smsq = 0;
sm3pw = 0;
sm4pw = 0;
n = 0;
for ii = 1:l
    minindx = max(ii - dx, 1);
    maxindx = min(ii + dx, l);
    if (ii == 1)
        subsrs = series(minindx:maxindx);
        subsrsMean = rnmn(1);
        subSqSrs = (subsrs - subsrsMean).^2;
        sub3pwSrs = (subsrs - subsrsMean).^3;
        sub4pwSrs = (subsrs - subsrsMean).^4;
        smsq = nansum(subSqSrs);
        sm3pw = nansum(sub3pwSrs);
        sm4pw = nansum(sub4pwSrs);
        n = sum(~isnan(subSqSrs));
    else
        if minindx > 1
            sprev = series(minindx - 1) - rnmn(minindx - 1);
            if ~isnan(sprev)
                smsq = max(0, smsq - sprev^2);
                sm3pw = sm3pw - sprev^3;
                sm4pw = max(0, sm4pw - sprev^4);
                n = n - 1;
            end
        end
        if maxindx < l
            snext = series(maxindx + 1) - rnmn(minindx + 1);
            if ~isnan(snext)
                sm = sm + snext;
                smsq = smsq + snext^2;
                sm3pw = sm3pw + snext^3;
                sm4pw = sm4pw + snext^4;
                n = n + 1;
            end
        end
    end
      
    if n > minNThreshold
        rnvar(ii) = smsq/n;
        rn3mom(ii) = sm3pw/n;
        rn4mom(ii) = sm4pw/n;
    else
        rnvar(ii) = nan;
        rn3mom(ii) = nan;
        rn4mom(ii) = nan;
    end
end

end

