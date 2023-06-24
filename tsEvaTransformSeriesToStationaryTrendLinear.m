function [trasfData] = tsEvaTransformSeriesToStationaryTrendLinear( timeStamps, series, timeWindow, percentile, varargin )
% this function first calculates linear trend of the series and also return
% linear trend of percetile series (at any particular percentile level set
% by the user) through a fitting algorithm
disp('computing the trend ...');
[statSeries, trendSeries, filledTimeStamps, filledSeries, nRunMn,trendError] = tsEvaDetrendTimeSeries(timeStamps, series, timeWindow, varargin{:});

disp(['computing the slowly varying ' num2str(percentile) 'th percentile ...']);
TimeVec=datevec(filledTimeStamps);
[~,ia,~]=unique(TimeVec(:,1),'stable'); %this ensures selection of yearly data

PercentileSeries=zeros(length(ia),1);
ia(end+1)=length(filledTimeStamps)+1;
for ij=1:length(ia)-1
    filledSeriesYear=filledSeries(ia(ij):ia(ij+1)-1);
    PercentileValYear = prctile(filledSeriesYear, percentile);
    PercentileSeries(ij)=PercentileValYear;
end
ia=ia(1:end-1);
tvec=datevec(filledTimeStamps(ia));
tvec=[tvec(:,1),7*ones(length(tvec),1),1*ones(length(tvec),1),zeros(length(tvec),3)];
TimeStampsciPercentile=datenum(tvec);
[p,S] = polyfit(TimeStampsciPercentile,PercentileSeries,1);
PercentileSeriesTotal=p(1).*filledTimeStamps+p(2);
cy=(PercentileSeriesTotal-trendSeries); 
statSeries=statSeries./(cy);
[~,p_value]=tsMann_Kendall(PercentileSeries,0.05);
meanPerc = nanmean(PercentileSeriesTotal);
%normalizing to standard deviation (just to be able to make acceptable graphs with the scripts of this library)
stdDevSeries = cy;%percentileSeries/meanPerc.*cy;
[~,ErrorPercentile] = polyval(p,filledTimeStamps,S);
stdErr=sqrt((ErrorPercentile).^2+(trendError).^2);

% statSeries = statSeries./stdDevSeries;
[~, ~, statSer3Mom, statSer4Mom] = tsEvaNanRunningStatistics(statSeries, nRunMn);
statSer3Mom = tsEvaNanRunningMean(statSer3Mom, ceil(nRunMn));
statSer4Mom = tsEvaNanRunningMean(statSer4Mom, ceil(nRunMn));

% N is the size of each sample used to compute the average
N = nRunMn;
% the error on the trend is computed as the error on the average:
%  error(average) = stdDev/sqrt(N)
% trendError = nanmean(stdDevSeries)/N^.5;

trasfData.runningStatsMulteplicity = 0;%nRunMn;
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
