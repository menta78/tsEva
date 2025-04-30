% Simplified main function structure for tsCopulaPlotBivariate
function [axxArray] = tsCopulaPlotBivariate(copulaAnalysis, gofStatistics, varargin)

    %% Parse optional arguments
    args.xlbl = 'Date (time)';
    args.ylbl = {'X1', "X2"};
    args.fontSize = 14;
    args.rpPlot=[];
    
    args = tsEasyParseNamedArgs(varargin, args);
    
    xlbl = args.xlbl;
    ylbl = args.ylbl;
    fontSize = args.fontSize;
    rpPlot = args.rpPlot;

    %% Unpack key variables

    % Extract peaks and timestamps
    if copulaAnalysis.timeVaryingCopula
        yMax = copulaAnalysis.yMax;
        tMax = copulaAnalysis.tMax;
    else
        yMax = copulaAnalysis.jointExtremes;
        tMax = copulaAnalysis.jointExtremeTimeStamps;
    end
    scatterColor = mean(yMax, 2);
    
    [~, iyMax] = sort(scatterColor, 'ascend');
    yMax = yMax(iyMax,:);
    tMax = tMax(iyMax,:);
    scatterColor = scatterColor(iyMax);

    % Retrieve original series and timestamps
    marginalAnalysis = copulaAnalysis.marginalAnalysis;
    nonStatSeriesCell = cellfun(@(x) x{2}.nonStatSeries, marginalAnalysis, 'UniformOutput', false);
    timeStampsCell = cellfun(@(x) x{2}.timeStamps, marginalAnalysis, 'UniformOutput', false);
    nonStatSeries = cell2mat(nonStatSeriesCell);
    timeStamps = cell2mat(timeStampsCell);
    pval = cellfun(@(x) x{2}.pValueChange, copulaAnalysis.marginalAnalysis);
    thresholdPotNS = copulaAnalysis.thresholdPotNS;

    timeStampsByTimeWindow = copulaAnalysis.copulaParam.timeStampsByTimeWindow;
    ttRho = copulaAnalysis.copulaParam.rhoTimeStamps;
    couplingParam = copulaAnalysis.copulaParam.rho;

    %% Initialize subplot manager and figure
    rt = 1;
    b0 = 27;
    l0 = 27;
    spMan = tsLcSubplotManager(b0, l0, 'CellXSize', round(27*rt), 'CellYSize', round(27*rt), 'gap', [0 0]);
    spMan.initFigure;

    %% Define axes layout parameters
    h = [7,7,7,11,11];
    b = [13,13,13,11,11];
    h0 = [7,16,25,11.5,24.5];
    b0 = [2,2,2,17.5,17.5];

    axxArray = gobjects(1,0); % initialize empty array of axes

    %% Plot 1: time series 1 (panel a)
    axx = spMan.createAxes('ts1', h0(1), b0(1), h(1), b(1));
    axes(axx); axxArray(end+1) = axx;
    plotTimeSeries(timeStamps(:,1), nonStatSeries(:,1), thresholdPotNS(:,1), ...
        tMax(:,1), yMax(:,1), scatterColor, ...
        pval(1), xlbl, ylbl{1}, fontSize);

    %% Plot 2: time series 2 (panel b)
    axx = spMan.createAxes('ts2', h0(2), b0(2), h(2), b(2));
    axes(axx); axxArray(end+1) = axx;
    plotTimeSeries(timeStamps(:,2), nonStatSeries(:,2), thresholdPotNS(:,2), ...
        tMax(:,2), yMax(:,2), scatterColor, ...
        pval(2), xlbl, ylbl{2}, fontSize);

    %% Plot the time varying (or not varying) coupling parameters (panel c)
    axx = spMan.createAxes('copulaParam', h0(3), b0(3), h(3), b(3));
    axes(axx); axxArray(end+1) = axx;
    plotCopulaAnalysis(copulaAnalysis, gofStatistics, fontSize);

    %% Plot Montecarlo at the beginning of the series (panel d)
    axx = spMan.createAxes('mc1', h0(4), b0(4), h(4), b(4));
    axes(axx); axxArray(end+1) = axx;    
    plotMonteCarlo(copulaAnalysis, 1, fontSize);

    %% Plot Montecarlo at the end of the series (panel d)
    axx = spMan.createAxes('mc2', h0(4), b0(4), h(4), b(4));
    axes(axx); axxArray(end+1) = axx;
    plotMonteCarlo(copulaAnalysis, length(copulaAnalysis.monteCarloRsmpl), fontSize);

%% Plot Time Series with Thresholds
function plotTimeSeries(timeStamps, nonStatSeries, threshold, ...
        tPeaks, peaks, scatterColor, pval, xlbl, ylbl, fontSize)
    hold on; box on;
    plot(timeStamps, nonStatSeries);
    plot(timeStamps, threshold, '--r');
    xlim( [min(timeStamps), max(timeStamps)] );
    cmap = jet(256);
    minC = min(scatterColor);
    maxC = max(scatterColor);
    colorIdx = round(1 + (scatterColor - minC) / (maxC - minC) * (size(cmap,1)-1));
    colorIdx = max(1, min(colorIdx, size(cmap,1))); % Clamp values to valid index range
    scatter(tPeaks, peaks, [], cmap(colorIdx,:), 'filled');
    datetick('x','yyyy','keeplimits');
    xlabel(xlbl, 'FontSize', fontSize);
    ylabel(ylbl, 'FontSize', fontSize);
    text(0.5,0.1,['{\it p-value}_{MK}= ',sprintf('%0.3g', pval)], 'units', 'normalized');
    grid on;
end

%% plot the time-varying (or fixed in time) copula parameters and GOF
function plotCopulaAnalysis(copulaAnalysis, gofStatistics, fontSize)
    ttRho = copulaAnalysis.copulaParam.rhoTimeStamps;
    couplingParam = copulaAnalysis.copulaParam.rho;
    corrSpearmanSamplex = gofStatistics.corrSpearmanSamplex;
    corrSpearmanMontex = gofStatistics.corrSpearmanMontex;
    x11 = datetime(datevec(ttRho));
    y11 = [cell2mat(couplingParam);cell2mat(couplingParam)];
    x22 = [datetime(datevec(ttRho)),datetime(datevec(ttRho))];
    y22 = [[corrSpearmanSamplex{:};corrSpearmanSamplex{:}],[corrSpearmanMontex{:};corrSpearmanMontex{:}]];
    [hAx,hLine1,hLine2] = plotyy(x11,y11,x22,y22);
    set(hAx,{'ycolor'},{'k';'k'})
    hLine1.LineWidth = 1;
    hLine2(1).LineWidth = 1;
    hLine2(2).LineWidth = 1;
    str = lower(char(copulaFamily));
    idx = regexp([' ' str],'(?<=\s+)\S','start')-1;
    str(idx) = upper(str(idx));
    yll = ylabel(hAx(1),sprintf('\\theta_{%s}', str));

    %positionLabel2=[positionLabel2;yll.Position];
    %yLabel2=[yLabel2,yll];
end

%% Plots the montecarlo for the time timeIndex + the & return period
function plotMonteCarlo(copulaAnalysis, timeIndex, fontSize)
    monteCarloRsmpl = copulaAnalysis.monteCarloRsmpl{timeIndex};

end
