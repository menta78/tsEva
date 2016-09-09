function phandles = tsEvaPlotSeriesTrendStdDev( timeStamps, series, trend, stdDev, varargin )

args.confidenceAreaColor = [189 252	201]/255*.8; % light green
args.confidenceBarColor = [34 139 34]/255; % dark green
args.seriesColor = [1, .5, .5];
args.trendColor = 'k';
args.xlabel = '';
args.ylabel = 'level (m)';
args.minYear = -7000;
args.maxYear = 11000;
args.title = '';
args.axisFontSize = 22;
args.labelFontSize = 28;
args.titleFontSize = 30;
args.legendLocation = 'northwest';
args.dateformat = 'yyyy';
args.figPosition = [0, 0, 1300, 700] + 10;
args.verticalRange = [];
args.statsTimeStamps = timeStamps;
args.xtick = [];

args = tsEasyParseNamedArgs(varargin, args);

statsTimeStamps = args.statsTimeStamps;
minTS = datenum([args.minYear, 1, 1]);
maxTS = datenum([args.maxYear, 1, 1]);
series = series( (timeStamps >= minTS) & (timeStamps <= maxTS) );
trend = trend( (statsTimeStamps >= minTS) & (statsTimeStamps <= maxTS) );
stdDev = stdDev( (statsTimeStamps >= minTS) & (statsTimeStamps <= maxTS) );
timeStamps = timeStamps( (timeStamps >= minTS) & (timeStamps <= maxTS) );

upCI = trend + stdDev;
downCI = trend - stdDev;

f = figure;
f.Position = args.figPosition;
phandles{1} = f;

phandles{2} = plot(timeStamps, series, 'color', args.seriesColor, 'linewidth', .5);
datetick('x', args.dateformat);
if ~isempty(args.xtick)
    set(gca, 'xtick', args.xtick);
    set(gca, 'xticklabel', datestr(args.xtick, args.dateformat));
end
xlim([min(timeStamps) max(timeStamps)]);
hold on;

xcibar = cat(1, statsTimeStamps, flip(statsTimeStamps));
ycibar = cat(1, upCI, flip(downCI));
fl = fill(xcibar, ycibar, args.confidenceAreaColor);
fl.FaceAlpha = .6;
fl.EdgeColor = 'none';
phandles{3} = fl;

%h = area(timeStamps, cat(2, downCI, upCI - downCI), 'linestyle', 'none');
%h(1).FaceColor = 'none';
%h(2).FaceColor = args.confidenceAreaColor;
%alpha(.5);

phandles{4} = plot(statsTimeStamps, trend, 'color', args.trendColor, 'linewidth', 3);
phandles{5} = plot(statsTimeStamps, upCI, 'color', args.confidenceBarColor, 'linewidth', 2);
phandles{6} = plot(statsTimeStamps, downCI, 'color', args.confidenceBarColor, 'linewidth', 2);
grid on;

if ~isempty(args.verticalRange)
    ylim(args.verticalRange);
end

legend([phandles{2} phandles{4} phandles{6}], {'Series' 'Trend' 'Std dev'}, 'fontsize', args.labelFontSize, 'location', args.legendLocation);

set(gca, 'fontsize', args.axisFontSize);

xlabel(args.xlabel, 'fontsize', args.labelFontSize);
ylabel(args.ylabel, 'fontsize', args.labelFontSize);

if ~isempty(args.title)
    title(args.title, 'fontsize', args.titleFontSize);
end
set(f, 'paperpositionmode', 'auto');

hold off;

end

