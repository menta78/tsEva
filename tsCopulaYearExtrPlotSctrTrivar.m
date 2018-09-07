function handles = tsCopulaYearExtrPlotSctrTrivar(resampleLevel, yMaxLevel, varargin)

args.xlbl = 'X';
args.ylbl = 'Y';
args.zlbl = 'Z';
args.xaxiscolor = [.7, .4, 0.05];
args.yaxiscolor = [0.3, 0.4, 0.1];
args.zaxiscolor = [0, 0.5, 0.5];
args.figPosition = [50, 50, 850, 850];
args.fontSize = 12;
args.azimuth = -74;
args.elevation = 35;
args.title = '';
args.titleFontSize = 15;
args.markerTranspAlpha = .3;
args = tsEasyParseNamedArgs(varargin, args);
xlbl = args.xlbl;
ylbl = args.ylbl;
zlbl = args.zlbl;
figPosition = args.figPosition;
fontSize = args.fontSize;
azimuth = args.azimuth;
elevation = args.elevation;
xaxiscolor = args.xaxiscolor;
yaxiscolor = args.yaxiscolor;
zaxiscolor = args.zaxiscolor;
title = args.title;
titleFontSize = args.titleFontSize;
markerTranspAlpha = args.markerTranspAlpha;

if (size(resampleLevel, 2) ~= 3) || (size(yMaxLevel, 2) ~= 3)
  error(['tsCopulaYearExtrPlotSctrTrivar: resampleLevel and yMaxLevel must be Nx3 arrays']);
end

fig = figure('position', figPosition, 'color', 'w');
rsmplSctr = scatter3(resampleLevel(:, 1), resampleLevel(:, 2), resampleLevel(:, 3), 'MarkerFaceAlpha', markerTranspAlpha, 'MarkerEdgeAlpha', markerTranspAlpha);
hold on;
ymaxSctr = scatter3(yMaxLevel(:,1), yMaxLevel(:,2), yMaxLevel(:, 3), 'markerfacecolor', 'r');
view(azimuth, elevation);

grid on;

axsc3d = gca;
axsc3d.FontSize = fontSize;
axsc3d.Position = [.35, .4, .5, .6];
axsc3d.XAxis.Color = xaxiscolor;
axsc3d.YAxis.Color = yaxiscolor;
axsc3d.ZAxis.Color = zaxiscolor;

lgnd = legend([ymaxSctr, rsmplSctr], {'Yearly Maxima', ['Joint Distribution' newline 'Montecarlo']}, 'fontsize', fontSize, 'box', 'off');
lgnd.Position = [.08, .85, .25, .12];

axHist1 = axes('position', [.55, .11, .4, .25]);
h1 = histogram(resampleLevel(:, 1));
axHist1.XAxis.Color = xaxiscolor;
axHist1.YAxis.Visible = 'off';
axHist1.FontSize = fontSize;
xlm = xlim;
ylm = ylim;
xxlbl = (xlm(1) + xlm(2))/2;
yylbl = ylm(1) + .9*(ylm(2) - ylm(1));
text(xxlbl, yylbl, xlbl, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline', 'fontsize', fontSize, 'color', xaxiscolor);
box off;

axHist2 = axes('position', [.1, .15, .4, .25]);
h2 = histogram(resampleLevel(:, 2));
axHist2.XAxis.Color = yaxiscolor;
axHist2.YAxis.Visible = 'off';
axHist2.FontSize = fontSize;
xlm = xlim;
ylm = ylim;
xxlbl = (xlm(1) + xlm(2))/2;
yylbl = ylm(1) + .9*(ylm(2) - ylm(1));
text(xxlbl, yylbl, ylbl, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline', 'fontsize', fontSize, 'color', yaxiscolor);
box off;


axHist3 = axes('position', [.04, .45, .25, .35]);
h3 = histogram(resampleLevel(:, 3));
axHist3.XAxis.Color = yaxiscolor;
axHist3.YAxis.Visible = 'off';
axHist3.XAxis.Visible = 'off';
axHist3.FontSize = fontSize;
h3.Orientation = 'horizontal';
axHist3.XDir = 'reverse';
axHist3.XTick = [];
xlm = xlim;
ylm = ylim;
xxlbl = (xlm(1) + xlm(2))/2;
yylbl = ylm(1) + .9*(ylm(2) - ylm(1));
text(xxlbl, yylbl, zlbl, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline', 'fontsize', fontSize, 'color', zaxiscolor);
box off;

if ~strcmpi(title, '')
  axTitle = axes('position', [0, 0, 1, .07]);
  axTitle.XAxis.Visible = 'off';
  axTitle.YAxis.Visible = 'off';
  box off;
  text(.5, 1, title, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'fontsize', titleFontSize);
else
  axTitle = [];
end

set(fig, 'paperpositionmode', 'auto');

handles = [fig, ymaxSctr, rsmplSctr, h1, h2, h3, axTitle];
