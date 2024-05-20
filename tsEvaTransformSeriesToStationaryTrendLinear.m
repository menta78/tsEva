function [trasfData] = tsEvaTransformSeriesToStationaryTrendLinear( timeStamps, series, timeWindow, percentile, varargin )
% this function first calculates linear trend of the series and also return
% linear trend of percetile series (at any particular percentile level set
% by the user) through a fitting algorithm

args.extremeLowThreshold = -Inf;
args = tsEasyParseNamedArgs(varargin, args);
extremeLowThreshold = args.extremeLowThreshold;

disp('computing the linear trend ...');
disp(['computing the slowly varying ' num2str(percentile) 'th percentile ...']);

[filledTimeStamps, filledSeries, dt ] = tsEvaFillSeries( timeStamps, series );
nRunMn = ceil(timeWindow/dt); %this line is kept just for sake of code consistency

% calculate linear trend
timeVec=datevec(filledTimeStamps);
[~,ia,~]=unique(timeVec(:,1),'stable'); %this ensures selection of yearly data
yearlyAveragedSeries=zeros(length(ia),1);
percentileSeries=zeros(length(ia),1);
ia(end+1)=length(filledTimeStamps)+1;
for ij=1:length(ia)-1

    filledSeriesYear=filledSeries(ia(ij):ia(ij+1)-1);
    averagedValYear = mean(filledSeriesYear,'omitnan');
    yearlyAveragedSeries(ij)=averagedValYear;

    percentileValYear = prctile(filledSeriesYear, percentile);
    percentileSeries(ij)=percentileValYear;
    
end
ia=ia(1:end-1);
tvec=datevec(filledTimeStamps(ia));
tvec=[tvec(:,1),7*ones(size(tvec,1),1),ones(size(tvec,1),1),zeros(size(tvec,1),3)];
timeStampsYearly=datenum(tvec);
[p,S] = polyfit(timeStampsYearly,yearlyAveragedSeries,1);
[p1,S1] = polyfit(timeStampsYearly,percentileSeries,1);
trendSeries=p(1).*filledTimeStamps+p(2);
[~,trendErrorSeries] = polyval(p,filledTimeStamps,S);
trendError=mean(trendErrorSeries);
filledSeries(filledSeries < extremeLowThreshold) = nan;
detrendSeries = filledSeries - trendSeries;

PercentileSeriesTotal=p1(1).*filledTimeStamps+p1(2);
cy=(PercentileSeriesTotal-trendSeries);

statSeries=detrendSeries./(cy);
[~,p_value]=tsMann_Kendall(percentileSeries,0.05);
stdDevSeries = cy;
[~,ErrorPercentile] = polyval(p1,filledTimeStamps,S1);
stdErr=sqrt((ErrorPercentile).^2+(trendError).^2);


[~, ~, statSer3Mom, statSer4Mom] = tsEvaNanRunningStatistics(statSeries, nRunMn);
statSer3Mom = tsEvaNanRunningMean(statSer3Mom, ceil(nRunMn));
statSer4Mom = tsEvaNanRunningMean(statSer4Mom, ceil(nRunMn));

trasfData.runningStatsMulteplicity = 0;
trasfData.stationarySeries = statSeries;   
trasfData.trendSeries = trendSeries;  
trasfData.trendSeriesNonSeasonal = trendSeries;
trasfData.trendError = trendError; 
trasfData.stdDevSeries = cy;
trasfData.stdDevSeriesNonSeasonal = stdDevSeries;
trasfData.stdDevError = stdErr;
trasfData.timeStamps = filledTimeStamps;  
trasfData.nonStatSeries = filledSeries;
trasfData.statSer3Mom = statSer3Mom;
trasfData.statSer4Mom = statSer4Mom;
trasfData.pValueChange=p_value;
end
