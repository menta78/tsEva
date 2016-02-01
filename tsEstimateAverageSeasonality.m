function averageSeasonalitySeries = tsEstimateAverageSeasonality( timeStamps, seasonalitySeries )
  avgYearLength = 365.2425;
  nMonthInYear = 12;
  avgMonthLength = avgYearLength/nMonthInYear;

  firstTmStmp = timeStamps(1);
  lastTmStmp = timeStamps(end);
  monthTmStampStart = firstTmStmp:avgMonthLength:lastTmStmp;
  monthTmStampEnd = monthTmStampStart + avgMonthLength;
  
  grpdSsn_ = cellfun(@(i) nanmean(seasonalitySeries((monthTmStampStart(i) <= timeStamps) & (timeStamps < monthTmStampEnd(i)))), num2cell(1:length(monthTmStampStart)), 'uniformoutput', false);
  nYears = ceil(length(grpdSsn_)/nMonthInYear);
  grpdSsn = nan*ones(nYears*nMonthInYear, 1);
  grpdSsn(1:length(grpdSsn_)) = [grpdSsn_{:}];
  
  grpdSsnMtx = reshape(grpdSsn, 12, []);
  mnSsn_ = nanmean(grpdSsnMtx, 2);
  
  %estimating the first 2 fourier components
  it = (1:nMonthInYear)';
  x = it/6.*pi;
  dx = pi/6.;
  a0 = mean(mnSsn_);
  a1 = 1/pi * sum(cos(x).*mnSsn_)*dx;
  b1 = 1/pi * sum(sin(x).*mnSsn_)*dx;
  a2 = 1/pi * sum(cos(2*x).*mnSsn_)*dx;
  b2 = 1/pi * sum(sin(2*x).*mnSsn_)*dx;
  a3 = 1/pi * sum(cos(3*x).*mnSsn_)*dx;
  b3 = 1/pi * sum(sin(3*x).*mnSsn_)*dx;
  
  mnSsn = a0  +  ( a1*cos(x) + b1*sin(x) )  +  ( a2*cos(2*x) + b2*sin(2*x) )  +  ( a3*cos(3*x) + b3*sin(3*x) );
  monthAvgMtx = mnSsn(:,ones(nYears,1));
  monthAvgVec = reshape(monthAvgMtx, nMonthInYear*nYears, []);

  imnth = (0:length(monthAvgVec)-1)';
  avgTmStamp = firstTmStmp + avgMonthLength/2. + imnth*avgMonthLength;
  
  %adding first and last times
  monthAvgVec = cat(1, monthAvgVec(1), monthAvgVec, monthAvgVec(end));
  avgTmStamp = cat(1, firstTmStmp, avgTmStamp, max(monthTmStampEnd));
  
  averageSeasonalitySeries = interp1(avgTmStamp, monthAvgVec, timeStamps);
end

