function [axxArray] = tsCopulaPlotTrivariate(copulaAnalysis, monteCarloAnalysis, varargin)

% tsCopulaPlotTrivariate  plotting joint peaks and Monte-Carlo resampled
%values

% [axxArray] = tsCopulaPlotTrivariate(copulaAnalysis,gofStatistics,varargin)
% used for trivariate plotting of joint extremes and time-variation of the
% coupling parameter



% input:
%  copulaAnalysis          - a variable of type structure provided by
%                              calling tsCopulaCompoundGPD and
%                              tsCopulaCompoundGPDMontecarlo functions
%                              first
%  gofStatistics           - a variable of type structure provided as the
%                             output of tsCopulaGOFNonStat function


% output:
%  axxArray:                - a set of handles concerning the generated plots;
%
%
%
%
%

% M.H.Bahmanpour 2025

%REFERENCES

% [1] Bahmanpour, M.H., Mentaschi, L., Tilloy, A., Vousdoukas, M.,
%     Federico, I., Coppini, G., and Feyen, L., 2025,
%     Transformed-Stationary EVA 2.0: A Generalized Framework for
%     Non-stationary Joint Extreme Analysis (submitted to Hydrology and
%     Earth System Sciences; Feb 2025)

% setting the default parameters

args.xlbl = 'Date (time)';
args.fontSize = 12;
args.varLabels = ["Var1","Var2","Var3"];
args.gofStatistics = [];
args = tsEasyParseNamedArgs(varargin, args);

xlbl = args.xlbl;
fontSize = args.fontSize;
varLabels = args.varLabels;
gofStatistics = args.gofStatistics;

labelMark = (["(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)"]);
methodology = copulaAnalysis.methodology;

%% Initialize subplot manager and figure, set figure and axes properties;
%keep axes handles as an array
rt = 1;
b0 = 27;
l0 = 27;
spMan = tsLcSubplotManager(b0, l0, 'CellXSize', round(27*rt), 'CellYSize', round(27*rt), 'gap', [0 0]);
spMan.initFigure;

%% Define axes layout parameters
h = 5;
b = repmat(7.5,1,11); b(1)=11.25;b(5)=11.25;
h0 = [5,12,19,26,5,12,19,26,12,19,26];
b0 = [2,2,2,2,16,11.5,11.5,11.5,21,21,21];
% initialize empty array of axes
axxArray = gobjects(1,0);

%%
% Retrieve original series and timestamps
marginalAnalysis = copulaAnalysis.marginalAnalysis;
nonStatSeriesCell = cellfun(@(x) x{2}.nonStatSeries, marginalAnalysis, 'UniformOutput', false);
timeStampsCell = cellfun(@(x) x{2}.timeStamps, marginalAnalysis, 'UniformOutput', false);
nonStatSeries = cell2mat(nonStatSeriesCell);
timeStamps = cell2mat(timeStampsCell);
pval = cellfun(@(x) x{2}.pValueChange, copulaAnalysis.marginalAnalysis);
pvalStat = cellfun(@(x) x{2}.pValueChangeStat, copulaAnalysis.marginalAnalysis);
mc = monteCarloAnalysis.monteCarloRsmpl;
timeStampsByTimeWindow = copulaAnalysis.copulaParam.timeStampsByTimeWindow;
couplingParam = copulaAnalysis.copulaParam.rho;
couplingParamRaw = copulaAnalysis.copulaParam.rhoRaw; % used to compute the Mann-Kendall test
ttRho = copulaAnalysis.copulaParam.rhoTimeStamps;
t1xStrt = datestr(timeStampsByTimeWindow{1}(1),'yyyy');
t2xStrt = datestr(timeStampsByTimeWindow{1}(end),'yyyy');
t1xEnd=datestr(timeStampsByTimeWindow{end}(1),'yyyy');
t2xEnd=datestr(timeStampsByTimeWindow{end}(end),'yyyy');

couplingParamMat = (cell2mat(cellfun(@(x) x(find(tril(x,-1))), couplingParam, 'UniformOutput', 0)))';
couplingParamMat = arrayfun(@(col) couplingParamMat(:, col), 1:size(couplingParamMat, 2), 'UniformOutput', false);
couplingParamMat = cell2mat(couplingParamMat);  % Convert cell array back to matrix

jointExtremes = copulaAnalysis.jointExtremes;
couplingParamMatRaw = (cell2mat(cellfun(@(x) x(find(tril(x,-1))), couplingParamRaw, 'UniformOutput', 0)))';
couplingParamMatRaw = arrayfun(@(col) couplingParamMatRaw(:, col), 1:size(couplingParamMatRaw, 2), 'UniformOutput', false);
couplingParamMatRaw = cell2mat(couplingParamMatRaw);  % Convert cell array back to matrix
family = copulaAnalysis.copulaParam.family;

if strcmpi(family, 'gaussian')
    cplSymbol='\rho';
else
    cplSymbol='\theta';
end

if strcmp(methodology,'gpd')
    thresholds=copulaAnalysis.thresholdPotNS;
else
    thresholds=[];
end
% Extract peaks and timestamps
if copulaAnalysis.timeVaryingCopula
    yMax = copulaAnalysis.yMax;
    tMax = copulaAnalysis.tMax;
else
    yMax = copulaAnalysis.jointExtremes;
    tMax = copulaAnalysis.jointExtremeTimeStamps;
end

[~, iyMax] = sort(mean(yMax,2), 'ascend');
yMax = yMax(iyMax,:);
tMax = tMax(iyMax,:);

%%
% Plot 3D scatter
axx = spMan.createAxes('sc3', h0(1), b0(1), h(1), b(1));
axes(axx); axxArray(end+1) = axx;
scatter3Plot(yMax,varLabels,labelMark(1),gofStatistics)

%%

% Plot time series
axx = spMan.createAxes('ts1', h0(2), b0(2), h, b(2));
axes(axx); axxArray(end+1) = axx;
plotTimeSeries(timeStamps(:,1),nonStatSeries(:,1),thresholds(:,1),tMax(:,1),yMax(:,1),mean(yMax,2),varLabels(1),...
    fontSize,labelMark(2),...
    xlbl,pval(1),pvalStat(1))

%%

axx = spMan.createAxes('ts2', h0(6), b0(6), h, b(6));
axes(axx); axxArray(end+1) = axx;
plotTimeSeries(timeStamps(:,2),nonStatSeries(:,2),thresholds(:,2),tMax(:,2),yMax(:,2),mean(yMax,2),varLabels(2),...
    fontSize,labelMark(6),...
    xlbl,pval(2),pvalStat(2))
%%

axx = spMan.createAxes('ts3', h0(9), b0(9), h, b(9));
axes(axx); axxArray(end+1) = axx;
plotTimeSeries(timeStamps(:,3),nonStatSeries(:,3),thresholds(:,3),tMax(:,3),yMax(:,3),mean(yMax,2),varLabels(3),...
    fontSize,labelMark(9),...
    xlbl,pval(3),pvalStat(3))
%%
%plot gof series
axx = spMan.createAxes('tsgf', h0(5), b0(5), h, b(5));
axes(axx); axxArray(end+1) = axx;
plotCouplingSeries(ttRho,couplingParamMat,couplingParamMatRaw,fontSize,xlbl,cplSymbol,family,labelMark(5))

%%
% plot begining and edning samples overplotted with MC
axx = spMan.createAxes('sc12', h0(3), b0(3), h, b(3));
axes(axx); axxArray(end+1) = axx;
pairs = nchoosek(1:3, 2);  % generates [1 2; 1 3; 2 3]

scatter2DPlot(mc{1},pairs(1,:),jointExtremes{1},varLabels,...
    labelMark(3),couplingParam{1},t1xStrt,t2xStrt,family,cplSymbol,fontSize)

axx = spMan.createAxes('sc13', h0(7), b0(7), h, b(7));
axes(axx); axxArray(end+1) = axx;
scatter2DPlot(mc{1},pairs(2,:),jointExtremes{1},varLabels,...
    labelMark(7),couplingParam{1},t1xStrt,t2xStrt,family,cplSymbol,fontSize)

axx = spMan.createAxes('sc23', h0(10), b0(10), h, b(10));
axes(axx); axxArray(end+1) = axx;
scatter2DPlot(mc{1},pairs(3,:),jointExtremes{1},varLabels,...
    labelMark(10),couplingParam{1},t1xStrt,t2xStrt,family,cplSymbol,fontSize)


axx = spMan.createAxes('sc12e', h0(4), b0(4), h, b(4));
axes(axx); axxArray(end+1) = axx;
scatter2DPlot(mc{end},pairs(1,:),jointExtremes{end},varLabels,...
    labelMark(4),couplingParam{end},t1xEnd,t2xEnd,family,cplSymbol,fontSize)

axx = spMan.createAxes('sc13e', h0(8), b0(8), h, b(8));
axes(axx); axxArray(end+1) = axx;
scatter2DPlot(mc{end},pairs(2,:),jointExtremes{end},varLabels,...
    labelMark(8),couplingParam{end},t1xEnd,t2xEnd,family,cplSymbol,fontSize)

axx = spMan.createAxes('sc23e', h0(11), b0(11), h, b(11));
axes(axx); axxArray(end+1) = axx;
scatter2DPlot(mc{end},pairs(3,:),jointExtremes{end},varLabels,...
    labelMark(11),couplingParam{end},t1xEnd,t2xEnd,family,cplSymbol,fontSize)

end

function [] = scatter3Plot(X,label,panelLabel,gofStatistics)
scatter3(X(:,1),X(:,2),X(:,3), [], X(:,1), 'filled');
view(57.5, 30);
set(gca, 'XDir', 'reverse')
hold('on');
for iVrbl = 1:3
    axisLabel = label{iVrbl};
    switch iVrbl
        case 1
            xlabel(axisLabel);
        case 2
            ylabel(axisLabel);
        case 3
            zlabel(axisLabel);
    end
end

set(gca,'FontSize', 12);
text(0.05, 0.9, panelLabel, 'Units','normalized');

% Add GOF statistics
stats = {
    sprintf('$\\overline{\\Delta\\rho}_{\\mathrm{S}} = %.2g$', gofStatistics.corrSpearmanSampleDelta),
    sprintf('$\\overline{\\Delta\\tau}_{\\mathrm{K}} = %.2g$', gofStatistics.corrKendallSampleDelta),
    sprintf('$\\overline{S_n} = %.1g$', gofStatistics.snSample)
    };
txth=[];
for iVrbl = 1:length(stats)
    th=text(0.7, 0.46 - 0.1*(iVrbl-1), stats{iVrbl}, 'units','normalized', 'Interpreter','latex', 'FontSize', 14);
    txth=[txth,th];
end
texts = txth;
drawnow;

% combine text extents (still in axes-normalized units)
extents = cell2mat(get(texts,'Extent'));
xMin = min(extents(:,1));
yMin = min(extents(:,2));
xMax = max(extents(:,1)+extents(:,3));
yMax = max(extents(:,2)+extents(:,4));
pad = 0.01;
posAxNorm = [xMin-pad, yMin-pad, (xMax-xMin)+2*pad, (yMax-yMin)+2*pad];

% convert axes-normalized -> figure-normalized coordinates 
ax = ancestor(th,'axes');
fig = ancestor(ax,'figure');

% get pixel positions of axes and figure
axPix  = getpixelposition(ax, true);   % axes position in pixels (relative to figure)
figPix = getpixelposition(fig);        % figure position in pixels

% convert axes-normalized box to pixel coordinates in figure
rectPix = [ axPix(1) + posAxNorm(1)*axPix(3), ...
    axPix(2) + posAxNorm(2)*axPix(4), ...
    posAxNorm(3)*axPix(3), ...
    posAxNorm(4)*axPix(4) ];

% convert pixels -> figure-normalized (0–1)
posFigNorm = [rectPix(1)/figPix(3), rectPix(2)/figPix(4), ...
    rectPix(3)/figPix(3), rectPix(4)/figPix(4)];

% draw annotation rectangle in the right place
annotation('rectangle', posFigNorm, ...
    'EdgeColor','k','LineWidth',1.5,'FaceColor','none');

end

function [] = plotTimeSeries (tt,Series,thresholdPotNS,tMax,yMax,CC,varLabels,fontSize,labelMark,xlbl,...
    pval,pvalStat)

plot(datetime(datevec(tt)), Series); hold on;
if ~isempty(thresholdPotNS)
    plot(datetime(datevec(tt)), thresholdPotNS,'LineWidth',2)
end
scatter(datetime(datevec(tMax)), yMax, [], CC, 'filled');
ylabel(varLabels);
set(gca, 'FontSize', fontSize);
text(0.05, 0.9, labelMark, 'Units','normalized');
xlabel(xlbl);
grid on;
text(0.1675, 0.8356, ['{\it p-value}_{nonStat}= ', sprintf('%0.3g', pval)],'units','normalized')
text(0.1675, 0.9356, ['{\it p-value}_{Stat}= ', sprintf('%0.3g', pvalStat)],'units','normalized')
end

function []=plotCouplingSeries(ttRho,couplingParamMat,couplingParamMatRaw,fontSize,xlbl,cplSymbol,family,labelMark)

plot(datetime(datevec(ttRho)), couplingParamMat,'LineWidth',1);
grid on;

pairs = nchoosek(1:3, 2);  % generates [1 2; 1 3; 2 3]

for i=1:3

    pairLabel = sprintf('%d-%d', pairs(i,1), pairs(i,2));

    [~, p_value] = tsMann_Kendall(couplingParamMatRaw(:, i), 0.05);

    text(0.2682, 0.4755+0.1*(i-1), ...
        ['$p$-value$_{', cplSymbol, ',', pairLabel, '}= ', sprintf('%0.3g', p_value), '$'], ...
        'Units', 'normalized', 'HorizontalAlignment', 'left', ...
        'Interpreter', 'latex');

end
set(gca, 'FontSize', fontSize);
xlabel(xlbl);
legendLabels = arrayfun(@(i) sprintf('$%s_{%d-%d}$', ...
    cplSymbol, pairs(i,1), pairs(i,2)), 1:size(pairs,1), 'UniformOutput', false);

% Add legend with LaTeX interpreter
legend(legendLabels, 'Interpreter', 'latex', 'Location', 'best');

% Extract family name from cell
familyname = family{1};

ylbStr = "$" + cplSymbol + "_{" + familyname + "}$";

% Use in ylabel
ylabel(ylbStr, 'Interpreter','latex', 'Rotation',0, ...
       'HorizontalAlignment','right', 'FontSize', fontSize);
text(0.05, 0.9, labelMark, 'Units','normalized');

end

function [] = scatter2DPlot(mc,pairs,yMax,varLabels,...
    labelMark,couplingParam,t1xStrt,t2xStrt,family,cplSymbol,fontSize)

sc_ = scatter(mc(:,pairs(1)), mc(:,pairs(2)), 10, 'filled', 'MarkerFaceAlpha', 0.3);
set(sc_, 'LineWidth', 1, 'Marker', 'o',...
    'MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor',[0.65,0.65,0.65],...
    'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.6);

hold on;

scatter(yMax(:,pairs(1)),yMax(:,pairs(2)), [],...
    yMax(:,1), 'filled');

xlabel(varLabels{pairs(1)});
ylabel(varLabels{pairs(2)});

set(gca, 'FontSize', fontSize);

text(0.05, 0.9, labelMark, 'Units','normalized');
grid on;

pr = couplingParam(pairs(1), pairs(2));

% Escape cplSymbol and convert to string

% Convert other parts to string scalars as well
line1 = string(t1xStrt) + " - " + string(t2xStrt);
line2 = string(family) + " ($" + cplSymbol + " = " + num2str(pr,'%.2f') + "$)";

% Set title with LaTeX interpreter
ttl = title([line1; line2], 'FontSize', fontSize, 'Interpreter', 'latex');
ttl.VerticalAlignment = 'top';

axis('tight')
set(gca,'fontsize',fontSize)
end