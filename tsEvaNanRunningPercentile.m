function [ rnprcnt ] = tsEvaNanRunningPercentile( series, windowSize, percent, varargin )

  if windowSize > 2000
    args.percentDelta = 1.;
  elseif windowSize > 1000
    args.percentDelta = 2.;
  else
    error('window size cannot be less than 1000');
  end
  args.nLowLimit = 1000;
  args.smoothoutput = true;
  args = tsEasyParseNamedArgs(varargin, args);
  percentDelta = args.percentDelta;
  nLowLimit = args.nLowLimit;
  smoothoutput = args.smoothoutput;

  percentP = percent + percentDelta;
  if percentP > 100
    error(['max percent: ' num2str(100 - percentDelta)]);
  end

  percentM = percent - percentDelta;
  if percentM < 0
    error(['min percent: ' num2str(percentDelta)]);
  end

  rnprcnt = zeros([length(series), 1])*nan;
  dx = ceil(windowSize/2);
  l = length(series);

  probObj.isNull = true;
  for ii = 1:l
    minindx = max(ii - dx, 1);
    maxindx = min(ii + dx, l);
    if probObj.isNull
      subsrs = series(minindx:maxindx);
      probObj = initPercentiles(subsrs, percentM, percent, percentP);
    else
      if minindx > 1
        sprev = series(minindx - 1);
        
        %% removing element and reviewing probability
        if ~isnan(sprev)
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
        if ~isnan(snext)
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
  
      rnprcnt(ii) = prcntii;
    else
      probObj.isNull = true;
      %rnprcnt(ii) = nan; rnprcnt is already initialized to nan
    end
  end

  if smoothoutput
    % smoothing
    rnprcnt = tsEvaNanRunningMean(rnprcnt, windowSize);
  end
end


function probObj = initPercentiles(subsrs, percentM, percent, percentP)
%  coder.inline('always');
  probObj.isNull = false;

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