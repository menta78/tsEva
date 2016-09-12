function [ rnprcnt, stdError ] = tsEvaNanRunningPercentile( series, windowSize, percent, varargin )
% tsEvaNanRunningPercentile:
% computes a runnig percentile for a given series,
% using a window with a size given by windowSize.
%
% Input parameters:
% series: input series.
% windowSize: size of the window for the running percentile. Cannot be < 1000
% percent: percent level to which the percentile is compute.
% label parameters:
%    percentDelta: delta for the computation of a percentile interval
%                  around the requested percentage. If for example
%                  percent==90 and percentDelta==1, then the 89th, 90th and
%                  91st percentiles are computed. Default value: 1 if
%                  windowSize > 2000, 2 if 2000 > windowsize > 1000. 
%    nLowLimit: minimum number of non nan elements for a window for
%               percentile computation.
%
% Output parameters:
% rnprcnt: approximated running percentile.
%
% How it works:
% let's suppose that percent == 90.
% For the first window we compute the right percentile using matlab
% function prctile, for percentages 89, 90, 91.
% Then for each step, we update these percentages on the basis
% of the quitting values and incoming values,
% and interpolate an approximated percentile for the requested percentage.

  if windowSize > 2000
    args.percentDelta = 1.;
  elseif windowSize > 1000
    args.percentDelta = 2.;
  elseif windowSize > 100
    args.percentDelta = 5.;
  else
    error('window size cannot be less than 100');
  end
  args.nLowLimit = 100;
  args = tsEasyParseNamedArgs(varargin, args);
  percentDelta = args.percentDelta;
  nLowLimit = args.nLowLimit;

  percentP = percent + percentDelta;
  if percentP > 100
    error(['max percent: ' num2str(100 - percentDelta)]);
  end

  percentM = percent - percentDelta;
  if percentM < 0
    error(['min percent: ' num2str(percentDelta)]);
  end

  rnprcnt0 = zeros([length(series), 1])*nan;
  dx = ceil(windowSize/2);
  l = length(series);
  
  %% initializing probObj
  minindx = 1;
  maxindx = min(1 + dx, l);
  subsrs = series(minindx:maxindx);
  probObj = initPercentiles(subsrs, percentM, percent, percentP);
  rnprcnt0(1) = probObj.t;

  for ii = 2:l
    minindx = max(ii - dx, 1);
    maxindx = min(ii + dx, l);
    if minindx > 1
      sprev = series(minindx - 1);

      %% removing element and reviewing probability
      sprevNotNan = ~isnan(sprev);
      if sprevNotNan
        Nold = probObj.N;
        Nnew = probObj.N - 1;

        nle = sprev < probObj.tM;
        probObj.probM = (probObj.probM * Nold - nle) / Nnew;
        probObj.percentM = probObj.probM * 100;

        nle = sprev < probObj.t;
        probObj.prob = (probObj.prob * Nold - nle) / Nnew;
        probObj.percent = probObj.prob * 100;

        nle = sprev < probObj.tP;
        probObj.probP = (probObj.probP * Nold - nle) / Nnew;
        probObj.percentP = probObj.probP * 100;

        probObj.N = Nnew;
      end

    end

    if maxindx < l
      snext = series(maxindx + 1);

      %% adding element and reviewing probability    
      snextNotNan = ~isnan(snext);
      if snextNotNan
        Nold = probObj.N;
        Nnew = probObj.N + 1;

        nle = snext < probObj.tM;
        probObj.probM = (probObj.probM * Nold + nle) / Nnew;
        probObj.percentM = probObj.probM * 100;

        nle = snext < probObj.t;
        probObj.prob = (probObj.prob * Nold + nle) / Nnew;
        probObj.percent = probObj.prob * 100;

        nle = snext < probObj.tP;
        probObj.probP = (probObj.probP * Nold + nle) / Nnew;
        probObj.percentP = probObj.probP * 100;

        probObj.N = Nnew;
      end

    end

    cout1 = probObj.percentM > percent;
    cout2 = probObj.percentP < percent;
    outOfInterval = cout1 || cout2;
    if outOfInterval
      %disp('reinit');
      subsrs = series(minindx:maxindx);
      probObj = initPercentiles(subsrs, percentM, percent, percentP);
    end

    if probObj.N > nLowLimit
      %% interpolating percentile
      if percent == probObj.percentM
        prcntii = probObj.tM;
      elseif (probObj.percentM < percent) && (percent < probObj.percent)
        %prcntii = intp(probObj.percentM, probObj.percent, probObj.tM, probObj.t, percent);
        h1 = probObj.percent - percent;
        h2 = percent - probObj.percentM;
        prcntii = (h1*probObj.tM + h2*probObj.t)/(h1 + h2);
      elseif percent == probObj.percent
        prcntii = probObj.t;
      elseif (probObj.percent < percent) && (percent < probObj.percentP)
        %prcntii = intp(probObj.percent, probObj.percentP, probObj.t, probObj.tP, percent);
        h1 = probObj.percentP - percent;
        h2 = percent - probObj.percent;
        prcntii = (h1*probObj.t + h2*probObj.tP)/(h1 + h2);
      elseif percent == probObj.percentP
        prcntii = probObj.tP;
      end
      %%
  
      rnprcnt0(ii) = prcntii;
    else
      probObj.isNull = true;
      %rnprcnt(ii) = nan; rnprcnt is already initialized to nan
    end
  end

  % smoothing output
  rnprcnt = tsEvaNanRunningMean(rnprcnt0, windowSize);
  stdError = nanstd(rnprcnt0 - rnprcnt);
end


function probObj = initPercentiles(subsrs, percentM, percent, percentP)
%  coder.inline('always');
  probObj.percentM = percentM;
  probObj.percent = percent;
  probObj.percentP = percentP;
  
  probObj.probM = percentM/100.;
  probObj.prob = percent/100.;
  probObj.probP = percentP/100.;
  
  probObj.tM = prctile(subsrs, percentM);
  probObj.t = prctile(subsrs, percent);
  probObj.tP = prctile(subsrs, percentP);
  
  probObj.N = sum(~isnan(subsrs));
end