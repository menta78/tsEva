%%This function plots original series plus time-varying GEV return levels as
%%well as annual maxima (block maxima)

function phandles = tsPlotSeriesYearMaxGEVRetLevFromAnalysisObj( nonStationaryEvaParams, stationaryTransformData, varargin )

timestamps = stationaryTransformData.timeStamps;
series = stationaryTransformData.nonStatSeries;

args.legendLocation = 'northwest';
args.rlevel = 1;
args.annualMax = 1;
args = tsEasyParseNamedArgs(varargin, args);

% computing the return levels
epsilon = nonStationaryEvaParams(1).parameters.epsilon;
sigma = nonStationaryEvaParams(1).parameters.sigma;
mu = nonStationaryEvaParams(1).parameters.mu;
epsilonStdErr = nonStationaryEvaParams(1).paramErr.epsilonErr;
sigmaStdErr = nonStationaryEvaParams(1).paramErr.sigmaErr;
muStdErr = nonStationaryEvaParams(1).paramErr.muErr;
returnPeriods = [5, 10, 30, 100];
[rlevel, returnLevelsErrs] = tsEvaComputeReturnLevelsGEV(epsilon, sigma, mu, epsilonStdErr, sigmaStdErr, muStdErr, returnPeriods);

phandles{1} = figure;
set(gcf,'position',[87.4000  342.0000  960.6000  420.0000])
plot(timestamps, series);
hold on;

% Plot return levels
phandles{2} = plot(timestamps, rlevel(:,1), 'Color', 'r');
phandles{3} = plot(timestamps, rlevel(:,2), 'Color', 'g');
phandles{4} = plot(timestamps, rlevel(:,3), 'Color', 'b');
phandles{5} = plot(timestamps, rlevel(:,4), 'Color', 'k');

% Plot annual maxima (block maxima)
annualMaxIdx = nonStationaryEvaParams(1).objs.annualMaxIndexes;
phandles{6} = plot(timestamps(annualMaxIdx), series(annualMaxIdx), '*', 'MarkerSize', 8);

datetick('x', 12);
legend('Series', '5-yr', '10-yr', '30-yr', '100-yr', 'Annual Maxima', 'Location', args.legendLocation);
xlabel('Date (time)');
ylabel('Value');

end
