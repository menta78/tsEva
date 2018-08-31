function [handles] = tsCopulaYearExtrScatter2D(resampleLevel, yMaxLevel, varargin)

args.xlbl = 'X';
args.ylbl = 'Y';
args.figPosition = [100, 100, 800, 800];
args.fontSize = 15;
args = tsEasyParseNamedArgs(varargin, args);
xlbl = args.xlbl;
ylbl = args.ylbl;
figPosition = args.figPosition;
fontSize = args.fontSize;

if (size(resampleLevel, 2) ~= 2) || (size(yMaxLevel, 2) ~= 2)
  error(['tsCopulaYearExtrScatter2D: resampleLevel and yMaxLevel must be Nx2 arrays']);
end

fig = figure('position', figPosition);
rsmplSctr = scatterhist(resampleLevel(:,1), resampleLevel(:,2), 'direction', 'out');
rsmplSctrObj = findall(rsmplSctr(1), 'type', 'line');
rsmplDist1Ax = rsmplSctr(2);
rsmplDist2Ax = rsmplSctr(3);
hold on;
ymaxSctr = scatter(yMaxLevel(:,1), yMaxLevel(:,2), 'markerfacecolor', 'r');
xlabel(xlbl);
ylabel(ylbl);
grid on;
ax = gca;
ax.FontSize = fontSize;
lgnd = legend([ymaxSctr, rsmplSctrObj], {'Yearly Maxima', ['Joint Distribution' newline 'Montecarlo']}, 'fontsize', fontSize);
lgndPos = lgnd.Position;
newHigh = lgndPos(4)*2;
lgnd.Position = [rsmplDist2Ax.Position(1), rsmplDist1Ax.Position(2) + rsmplDist1Ax.Position(4) - newHigh, lgndPos(3), newHigh];
set(fig, 'paperpositionmode', 'auto');

handles = [fig, rsmplSctrObj, ymaxSctr, lgnd, rsmplSctr(1), rsmplSctr(2), rsmplSctr(3)];
