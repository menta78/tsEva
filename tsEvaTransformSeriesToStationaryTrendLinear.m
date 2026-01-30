function [trasfData] = tsEvaTransformSeriesToStationaryTrendLinear( timeStamps, series, timeWindow, percentile, varargin )
% this function first calculates linear trend of the series and also return
% linear trend of percetile series (at any particular percentile level set
% by the user) through a fitting algorithm

args.extremeLowThreshold = -Inf;
args = tsEasyParseNamedArgs(varargin, args);
extremeLowThreshold = args.extremeLowThreshold;

disp('computing the trend series using linear regression of yearly-averaged values...');
disp(['computing the percentile series using linear regression of yearly ' num2str(percentile) 'th percentile values...']);

%remove nan values, add nan where constant values are found, omit
%duplicates and use a uniformly increasing time stamp
[filledTimeStamps, filledSeries, dt ] = tsEvaFillSeries( timeStamps, series );

% nRunMn is used for calculating of statistical moments (useful in
% assessing stationarity)
nRunMn = ceil(timeWindow/dt); 

% Selection of yearly segments and pre-assigning
timeVec=datevec(filledTimeStamps);
[~,iAA,~]=unique(timeVec(:,1),'stable'); 
yearlyAveragedSeries=zeros(length(iAA),1);
yearlyMaxSeries=zeros(length(iAA),1);
percentileSeries=zeros(length(iAA),1);
filledTimeStampsYearly=zeros(length(iAA),1);
iAA(end+1)=length(filledTimeStamps)+1;

% loop through all yearly segments and calculate mean and percentile series
for ij=1:length(iAA)-1

    filledTimeStampsYearly(ij)=mean(filledTimeStamps(iAA(ij):iAA(ij+1)-1));

    filledSeriesYear=filledSeries(iAA(ij):iAA(ij+1)-1);

    yearlyAveragedSeries(ij)=mean(filledSeriesYear,'omitnan');
     yearlyMaxSeries(ij)=max(filledSeriesYear,[],'omitnan');
    percentileSeries(ij)=prctile(filledSeriesYear, percentile);
    
end

% perform linear regression 
idGood=find(~isnan(yearlyAveragedSeries));
[p,S] = polyfit(filledTimeStampsYearly(idGood),yearlyAveragedSeries(idGood),1);
[p1,S1] = polyfit(filledTimeStampsYearly(idGood),percentileSeries(idGood),1);

% assess Mann_Kendall analysis on percentile series
[~,p_value]=tsMann_Kendall(percentileSeries(idGood),0.05);
[~,pValueChangeAnnual]=tsMann_Kendall(yearlyMaxSeries(idGood),0.05);
% expand linear regression model to the entire length of series; also assess
% error in regression 
trendSeries=p(1).*filledTimeStamps+p(2);
[~,trendErrorSeries] = polyval(p,filledTimeStamps,S);
trendError=mean(trendErrorSeries);
% perfrom detrending
filledSeries(filledSeries < extremeLowThreshold) = nan;
detrendSeries = filledSeries - trendSeries;

%expand linear regression model to the entire length of series 
PercentileSeriesTotal=p1(1).*filledTimeStamps+p1(2);
% calculate std dev series 
stdDevSeries=(PercentileSeriesTotal-trendSeries);

% calculae stationary series and its p-value
statSeries=detrendSeries./(stdDevSeries);
for ij=1:length(iAA)-1

    statSeriesYear=statSeries(iAA(ij):iAA(ij+1)-1);
   
    percentileSeriesStat(ij)=prctile(statSeriesYear, percentile);
    
end
[~,p_valueStat]=tsMann_Kendall(percentileSeriesStat,0.05);


% assess tendencies
percentChangeTrend=((trendSeries(end)-trendSeries(1))/abs(trendSeries(1)))*100;
percentChangePercentile=((PercentileSeriesTotal(end)-PercentileSeriesTotal(1))/abs(PercentileSeriesTotal(1)))*100;

%assess error in regression model used for percentiles
[~,ErrorPercentile] = polyval(p1,filledTimeStamps,S1);

%assess error in std dev series (combination of error in assessing trend
%and error in assessing percentile series)
stdErr=sqrt((ErrorPercentile).^2+(trendErrorSeries).^2);

%assess thrid and fourth moment statistics of stationary series (also pass
%it through a running-mean algorithm)
[~, ~, statSer3Mom, statSer4Mom] = tsEvaNanRunningStatistics(statSeries, nRunMn);
statSer3Mom = tsEvaNanRunningMean(statSer3Mom, ceil(nRunMn));
statSer4Mom = tsEvaNanRunningMean(statSer4Mom, ceil(nRunMn));
%preparing the output
trasfData.runningStatsMulteplicity = 0;
trasfData.stationarySeries = statSeries;   
trasfData.trendSeries = trendSeries;  
trasfData.trendSeriesNonSeasonal = trendSeries;
trasfData.trendError = trendError; 
trasfData.stdDevSeries = stdDevSeries;
trasfData.stdDevSeriesNonSeasonal = stdDevSeries;
trasfData.stdDevError = stdErr;
trasfData.timeStamps = filledTimeStamps;  
trasfData.nonStatSeries = filledSeries;
trasfData.statSer3Mom = statSer3Mom;
trasfData.statSer4Mom = statSer4Mom;
trasfData.pValueChange=p_value;
trasfData.pValueChangeStat=p_valueStat;
trasfData.pValueChangeAnnual=pValueChangeAnnual;
trasfData.percentChangeTrend=percentChangeTrend;
trasfData.percentChangePercentile=percentChangePercentile;

end
