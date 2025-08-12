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

if isempty(gofStatistics)
    [gofStatistics] = tsCopulaGOFNonStat(copulaAnalysis, monteCarloAnalysis);
end

labelMark = (['_','(b)','(c)','(d)','(a)','(f)','(g)','(h)','(e)','(i)','(j)','(k)']);
methodology = copulaAnalysis.methodology;

axxArray = [];
rt = 1; hf = 27; wf = 27;
spMan = tsLcSubplotManager(hf, wf, 'CellXSize', round(27*rt), 'CellYSize', round(27*rt), 'gap', [0 0]);
spMan.initFigure;

h = repmat(5,1,12);
b = repmat(7.5,1,12); b(9)=15;
h0 = [5,12,19,26,5,12,19,26,5,12,19,26];
b0 = [2,2,2,2,11.5,11.5,11.5,11.5,21,21,21,21];
b0(5)=2; b0(9)=11.5;
hNone = gobjects(1,1);

for ij=1:length(h0)
    axx = hNone;
    if ij ~= 1
        axx = spMan.createAxes(num2str(ij), h0(ij), b0(ij), h(ij), b(ij));
    end
    axxArray = [axxArray, axx];
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

% Plot 3D scatter
scatter3(axxArray(5), yMax(:,1), yMax(:,2), yMax(:,3), [], yMax(:,1), 'filled');
view(57.5, 30);
hold(axxArray(5), 'on');
for iVrbl = 1:3
    axisLabel = varLabels{iVrbl};
    switch iVrbl
        case 1
            xlabel(axxArray(5), axisLabel);
        case 2
            ylabel(axxArray(5), axisLabel);
        case 3
            zlabel(axxArray(5), axisLabel);
    end
end

set(axxArray(5), 'FontSize', fontSize);
text(axxArray(5), 0.05, 0.9, labelMark(5), 'Units','normalized');

% Add GOF statistics
stats = {
    sprintf('$\\overline{\\Delta\\rho}_{\\mathrm{Spearman}} = %.2g$', gofStatistics.corrSpearmanSampleDelta),
    sprintf('$\\overline{\\Delta\\tau}_{\\mathrm{Kendall}} = %.2g$', gofStatistics.corrKendallSampleDelta),
    sprintf('$\\overline{S_n} = %.1g$', gofStatistics.snSample)
};
for iVrbl = 1:length(stats)
    text(axxArray(5), 0.65, 0.4 - 0.1*(iVrbl-1), stats{iVrbl}, 'units','normalized', 'Interpreter','latex', 'FontSize', 14);
end

% Retrieve original series and timestamps
marginalAnalysis = copulaAnalysis.marginalAnalysis;
nonStatSeriesCell = cellfun(@(x) x{2}.nonStatSeries, marginalAnalysis, 'UniformOutput', false);
timeStampsCell = cellfun(@(x) x{2}.timeStamps, marginalAnalysis, 'UniformOutput', false);
nonStatSeries = cell2mat(nonStatSeriesCell);
timeStamps = cell2mat(timeStampsCell);
pval = cellfun(@(x) x{2}.pValueChange, copulaAnalysis.marginalAnalysis);

iAx=[2,6,10];
for iVrbl = 1:3
    axes(axxArray(iAx(iVrbl)));
    plot(datetime(datevec(timeStamps(:,iVrbl))), nonStatSeries(:,iVrbl)); hold on;
    if strcmpi(methodology,'gpd')
        plot(datetime(datevec(timeStamps(:,iVrbl))), copulaAnalysis.thresholdPotNS(:,iVrbl),'LineWidth',2)
    end
    scatter(datetime(datevec(tMax(:,iVrbl))), yMax(:,iVrbl), [], mean(yMax,2), 'filled');
    ylabel(varLabels{iVrbl});
    set(gca, 'FontSize', fontSize);
    text(0.05, 0.9, labelMark(iAx(iVrbl)), 'Units','normalized');
    if iVrbl == 3
        xlabel(xlbl);
    end
    grid on;
    text(0.5, 0.1, ['{\it p-value}_{MK}= ', sprintf('%0.3g', pval(iVrbl))],'units','normalized')
end

yMaxLevel=copulaAnalysis.jointExtremes;
mc = monteCarloAnalysis.monteCarloRsmpl;
pairs = [1 2; 1 3; 2 3];

timeStampsByTimeWindow = copulaAnalysis.copulaParam.timeStampsByTimeWindow;
couplingParam = copulaAnalysis.copulaParam.rho;
couplingParamRaw = copulaAnalysis.copulaParam.rhoRaw; % used to compute the Mann-Kendall test
t1xStrt = datestr(timeStampsByTimeWindow{1}(1),'yyyy');
t2xStrt = datestr(timeStampsByTimeWindow{1}(end),'yyyy');
t1xEnd=datestr(timeStampsByTimeWindow{end}(1),'yyyy');
t2xEnd=datestr(timeStampsByTimeWindow{end}(end),'yyyy');

family = copulaAnalysis.copulaParam.family;
if iscell(family), family = family{1}; end
if strcmpi(family, "gaussian")
    cplSymbol='\rho'; 
else 
    cplSymbol='\alpha'; 
end;

iAx = [3, 7, 11];
for iCpl = 1:3
    axes(axxArray(iAx(iCpl)));
    sc_ = scatter(mc{1}(:,pairs(iCpl,1)), mc{1}(:,pairs(iCpl,2)), 10, 'filled', 'MarkerFaceAlpha', 0.3);
    set(sc_, 'LineWidth', 1, 'Marker', 'o',...
        'MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor',[0.65,0.65,0.65],...
        'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.6);
    hold on;
    scatter(yMaxLevel{1}(:,pairs(iCpl,1)), ...
        yMaxLevel{1}(:,pairs(iCpl,2)), [],...
        yMaxLevel{1}(:,1), 'filled');
    xlabel(varLabels{pairs(iCpl,1)});
    ylabel(varLabels{pairs(iCpl,2)});
    set(gca, 'FontSize', fontSize);
    text(0.05, 0.9, labelMark(5 + iCpl), 'Units','normalized');
    grid on;
    pr = couplingParam{1}(pairs(iCpl, 1), pairs(iCpl, 2));
    ttl = title( {[t1xStrt,' - ',t2xStrt]; [family ' (' cplSymbol ' = ',num2str(pr, '%.2f'),')']}, ...
        'fontsize', fontSize);
    ttl.VerticalAlignment = 'top';
end

iAx = [4, 8, 12];
for iCpl = 1:3
    axes(axxArray(iAx(iCpl)));
    sc_ = scatter(mc{end}(:,pairs(iCpl,1)), mc{end}(:,pairs(iCpl,2)), 10, 'filled', 'MarkerFaceAlpha', 0.3);
    set(sc_, 'LineWidth', 1, 'Marker', 'o',...
        'MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor',[0.65,0.65,0.65],...
        'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.6);
    hold on;
    scatter(yMaxLevel{end}(:,pairs(iCpl,1)), ...
        yMaxLevel{end}(:,pairs(iCpl,2)), [],...
        yMaxLevel{end}(:,1), 'filled');
    xlabel(varLabels{pairs(iCpl,1)});
    ylabel(varLabels{pairs(iCpl,2)});
    set(gca, 'FontSize', fontSize);
    text(0.05, 0.9, labelMark(5 + iCpl), 'Units','normalized');
    grid on;
    pr = couplingParam{end}(pairs(iCpl, 1), pairs(iCpl, 2));
    ttl = title( {[t1xEnd,' - ',t2xEnd]; [family ' (' cplSymbol ' = ',num2str(pr, '%.2f'),')']}, ...
        'fontsize', fontSize);
    ttl.VerticalAlignment = 'top';
end

ttRho = copulaAnalysis.copulaParam.rhoTimeStamps;

axes(axxArray(9));
couplingParamMat = (cell2mat(cellfun(@(x) x(find(tril(x,-1))), couplingParam, 'UniformOutput', 0)))';
couplingParamMat = arrayfun(@(col) couplingParamMat(:, col), 1:size(couplingParamMat, 2), 'UniformOutput', false);
couplingParamMat = cell2mat(couplingParamMat);  % Convert cell array back to matrix
plot(datetime(datevec(ttRho)), couplingParamMat);
grid on;

couplingParamMat = (cell2mat(cellfun(@(x) x(find(tril(x,-1))), couplingParamRaw, 'UniformOutput', 0)))';
couplingParamMat = arrayfun(@(col) couplingParamMat(:, col), 1:size(couplingParamMat, 2), 'UniformOutput', false);
couplingParamMat = cell2mat(couplingParamMat);  % Convert cell array back to matrix
for i = 1:3
    [~, p_value] = tsMann_Kendall(couplingParamMat(:, i), 0.05);
    text(0.35, 0.1*i, ...
        ['{\it p-value}_{MK,1-2}= ',sprintf('%0.3g',p_value)], ...
        'units', 'normalized', 'HorizontalAlignment', 'left');
end

for iAx = 2:length(h0)
    set(axxArray(iAx),'FontSize',fontSize);
    axis(axxArray(iAx),'tight');
end

end