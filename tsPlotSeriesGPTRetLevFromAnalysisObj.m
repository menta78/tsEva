%%This function plots original series plus time-varying return levels as
%%well as peaks

function phandles = tsPlotSeriesGPTRetLevFromAnalysisObj( nonStationaryEvaParams, stationaryTransformData, varargin )
args.plotPercentile = -1;
args = tsEasyParseNamedArgs(varargin, args);
plotPercentile = args.plotPercentile;
  
timestamps = stationaryTransformData.timeStamps;
series = stationaryTransformData.nonStatSeries;
trend = stationaryTransformData.trendSeries;
stdDev = stationaryTransformData.stdDevSeries;
if isfield(stationaryTransformData, 'statsTimeStamps')
  statsTimeStamps = stationaryTransformData.statsTimeStamps;
else
  statsTimeStamps = timestamps;
end

args.legendLocation = 'northwest';
args.rlevel=1;
args.pks=1;

args = tsEasyParseNamedArgs(varargin, args);
rlevel=args.rlevel

% phandles = tsEvaPlotSeriesTrendStdDev(timestamps, series, trend,stdDev, 'statsTimeStamps', statsTimeStamps, varargin{:});

% 
% if plotPercentile ~= -1
%   prcntile = tsEvaNanRunningPercentile(series, stationaryTransformData.runningStatsMulteplicity, plotPercentile);
%   figure(phandles{1});
%   hold on;
%   hndl = plot(timestamps, prcntile);
%   phandles{length(phandles) + 1} = hndl;
% end

hold off;

figure
set(gcf,'position',[87.4000  342.0000  960.6000  420.0000])
plot(timestamps,series)
hold on
% line([timestamps(1) timestamps(end)],rlevel(1)*ones(1,2),'Color','r')
plot(timestamps,rlevel(:,1),'Color','r')
plot(timestamps,rlevel(:,2),'Color','g')
plot(timestamps,rlevel(:,3),'Color','b')
plot(timestamps,rlevel(:,4),'Color','k')
% line([timestamps(1) timestamps(end)],rlevel(2)*ones(1,2),'Color','g')
% line([timestamps(1) timestamps(end)],rlevel(3)*ones(1,2),'Color','b')
% line([timestamps(1) timestamps(end)],rlevel(4)*ones(1,2),'Color','k')
pks=args.pks;
phandles=plot(timestamps(pks),series(pks),'*')
%  xlim([timestamps(1)-6*30 timestamps(end)+48*30])
%  ylim([min(series)-0.5 max(series)+1.2])
 datetick('x',12)
  legend('Series','5','10','30','100','peaks','location','bestoutside')
xlabel('Date (time)')
ylabel('TWL (m)')
end

