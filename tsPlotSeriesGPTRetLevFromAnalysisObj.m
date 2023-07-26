%%This function plots original series plus time-varying return levels as
%%well as peaks

function phandles = tsPlotSeriesGPTRetLevFromAnalysisObj( nonStationaryEvaParams, stationaryTransformData, varargin )

timestamps = stationaryTransformData.timeStamps;
series = stationaryTransformData.nonStatSeries;

args.legendLocation = 'northwest';
args.rlevel=1;
args.pks=1;

% computing the return levels
epsilon = nonStationaryEvaParams(2).parameters.epsilon;
sigma = nonStationaryEvaParams(2).parameters.sigma;
threshold = nonStationaryEvaParams(2).parameters.threshold;
thStart =  nonStationaryEvaParams(2).parameters.timeHorizonStart;
thEnd = nonStationaryEvaParams(2).parameters.timeHorizonEnd;
timeHorizonInYears = (thEnd - thStart)/365.2425;
nPeaks = nonStationaryEvaParams(2).parameters.nPeaks;
epsilonStdErr = nonStationaryEvaParams(2).paramErr.epsilonErr;
sigmaStdErr = nonStationaryEvaParams(2).paramErr.sigmaErr;
thresholdStdErr = nonStationaryEvaParams(2).paramErr.thresholdErr;
returnPeriods = [5, 10, 30, 100];
[rlevel, retrunLevelsErrs] = tsEvaComputeReturnLevelsGPD(epsilon, sigma, threshold, epsilonStdErr, sigmaStdErr, thresholdStdErr, nPeaks, timeHorizonInYears, returnPeriods);

phandles{1} = figure;
set(gcf,'position',[87.4000  342.0000  960.6000  420.0000])
plot(timestamps,series);
hold on;
% line([timestamps(1) timestamps(end)],rlevel(1)*ones(1,2),'Color','r')
phandles{2} = plot(timestamps,rlevel(:,1),'Color','r');
phandles{3} = plot(timestamps,rlevel(:,2),'Color','g');
phandles{4} = plot(timestamps,rlevel(:,3),'Color','b');
phandles{5} = plot(timestamps,rlevel(:,4),'Color','k');
pks = nonStationaryEvaParams(2).objs.peakIndexes;
phandles{6} = plot(timestamps(pks),series(pks),'*');
datetick('x',12);
legend('Series','5','10','30','100','peaks','location','bestoutside');
xlabel('Date (time)');
ylabel('TWL (m)');
end

