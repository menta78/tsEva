function [ detrendSeries, trendSeries, filledTimeStamps, filledSeries, nRunMn, TrendError ] = tsEvaDetrendTimeSeries( timeStamps, series, timeWindow, varargin )
% function to fill series and calculate detrendSeries. If transfType is set
% to trendlinear, linear trend of series is calculated using a simple
% fitting function. Else, trend of series is calculated using a running-mean formula

args.extremeLowThreshold = -Inf;
args.transfType = 'trend';
args = tsEasyParseNamedArgs(varargin, args);
extremeLowThreshold = args.extremeLowThreshold;
transfType=args.transfType;
[filledTimeStamps, filledSeries, dt ] = tsEvaFillSeries( timeStamps, series );

if strcmp(transfType,'trendlinear')
    nRunMn = ceil(timeWindow/dt); %this line is kept just for sake of code consistency
    TimeVec=datevec(filledTimeStamps);
    [~,ia,~]=unique(TimeVec(:,1),'stable'); %this ensures selection of yearly data
    YearlyAveragedSeries=zeros(length(ia),1);
    ia(end+1)=length(filledTimeStamps)+1;
    for ij=1:length(ia)-1
        filledSeriesYear=filledSeries(ia(ij):ia(ij+1)-1);
        AveragedValYear = nanmean(filledSeriesYear);
        YearlyAveragedSeries(ij)=AveragedValYear;
    end
    ia=ia(1:end-1);
    tvec=datevec(filledTimeStamps(ia));
    tvec=[tvec(:,1),7*ones(length(tvec),1),ones(length(tvec),1),zeros(length(tvec),3)];
    TimeStampsYearly=datenum(tvec);
    [p,S] = polyfit(TimeStampsYearly,YearlyAveragedSeries,1);
    trendSeries=p(1).*filledTimeStamps+p(2);
    [~,TrendErrorSeries] = polyval(p,filledTimeStamps,S);
    TrendError=mean(TrendErrorSeries);
else
    [trendSeries, filledTimeStamps, filledSeries, nRunMn] = tsEvaRunningMeanTrend(filledTimeStamps, filledSeries, timeWindow);
    
end
statSeries = filledSeries;
statSeries(statSeries < extremeLowThreshold) = nan;
detrendSeries = statSeries - trendSeries;
end

