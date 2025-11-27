function [axxArray] = tsCopulaPlotBivariate(copulaAnalysis, monteCarloAnalysis, varargin)
% tsCopulaPlotBivariate  plotting of joint peaks fitted by a copula

% [axxArray] = tsCopulaPlotBivariate(copulaAnalysis,gofStatistics,varargin)
% first two subplots show time series of input series along with joint
% peaks. The third subplot shows gof data. The last two subplots show
% scatter plots using the joint peaks and the fitted copula.



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

% M.H.Bahmanpour, 2025

%REFERENCES

% [1] Bahmanpour, M. H., Tilloy, A., Vousdoukas, M., Federico, I.,
% Coppini, G., Feyen, L., and Mentaschi, L.: Transformed-Stationary EVA 2.0:
% A Generalized Framework for Non-Stationary Joint Extremes Analysis,
% EGUsphere [preprint], https://doi.org/10.5194/egusphere-2025-843, 2025

% setting the default parameters

args.xlbl = 'Date (time)';
args.ylbl = {'Y1','Y2'};
args.fontSize = 14;
args.gofStatistics = [];
args.retPerAnalysis = [];

args = tsEasyParseNamedArgs(varargin, args);

xlbl = args.xlbl;
ylbl = args.ylbl;
fontSize = args.fontSize;
% set some parameters
labelMark=(["(a)","(b)","(c)","(d)","(e)"]);
gofStatistics = args.gofStatistics;
retPerAnalysis = args.retPerAnalysis;


if isempty(retPerAnalysis)
    retPerAnalysis = tsCopulaComputeBivarRP(copulaAnalysis, monteCarloAnalysis);
end


%% Unpack key variables
% Extract peaks and timestamps
yMax=copulaAnalysis.yMax;
tMax=copulaAnalysis.tMax;
nWindow=length(copulaAnalysis.jointExtremes);
methodology=copulaAnalysis.methodology; %***

% since we use non-stationary peaks as a way of color-coding the data,
% we'll have them sorted at this stage;

scatterColor = mean(yMax, 2);

[~, iyMax] = sort(scatterColor, 'descend');
yMax = yMax(iyMax,:);
tMax = tMax(iyMax,:);
scatterColor = scatterColor(iyMax);

% % Retrieve original series and timestamps, thresholds and pvals
marginalAnalysis=copulaAnalysis.marginalAnalysis;
nonStatSeriesCell=cellfun(@(x) x{2}.nonStatSeries,marginalAnalysis,'UniformOutput',0);
timeStampsCell=cellfun(@(x) x{2}.timeStamps,marginalAnalysis,'UniformOutput',0);
nonStatSeries=cell2mat(nonStatSeriesCell);
timeStamps=cell2mat(timeStampsCell);
thresholdPotNS=copulaAnalysis.thresholdPotNS;   %curves representing non-stationary threshold
pval=cellfun(@(x) x{2}.pValueChange,copulaAnalysis.marginalAnalysis);
pvalStat=cellfun(@(x) x{2}.pValueChangeStat,copulaAnalysis.marginalAnalysis);

%% Initialize subplot manager and figure, set figure and axes properties;
%keep axes handles as an array
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
% initialize empty array of axes
axxArray = gobjects(1,0);

%% Plot 1: time series 1 (panel a)
axx = spMan.createAxes('ts1', h0(1), b0(1), h(1), b(1));
axes(axx); axxArray(end+1) = axx;
plotTimeSeries(timeStamps(:,1),nonStatSeries(:,1),methodology,thresholdPotNS(:,1),...
    tMax(:,1),yMax(:,1),ylbl{1},...
    labelMark(1),xlbl,pval(1),pvalStat(1),fontSize,scatterColor);

%% Plot 2: time series 2 (panel b)
axx = spMan.createAxes('ts2', h0(2), b0(2), h(2), b(2));
axes(axx); axxArray(end+1) = axx;
plotTimeSeries(timeStamps(:,2),nonStatSeries(:,2),methodology,thresholdPotNS(:,2),...
    tMax(:,2),yMax(:,2),ylbl{2},...
    labelMark(2),xlbl,pval(2),pvalStat(2),fontSize,scatterColor);

%% Plot 3; goodness-of-fit statistics
axx = spMan.createAxes('gof', h0(3), b0(3), h(3), b(3));
axes(axx); axxArray(end+1) = axx;
gofPlot(copulaAnalysis, gofStatistics)

%% Plot Montecarlo at the beginning of the series (panel d)
axx = spMan.createAxes('mc1', h0(4), b0(4), h(4), b(4));
axes(axx); axxArray(end+1) = axx;
[xll,yll] = plotMonteCarlo(copulaAnalysis, monteCarloAnalysis, ...
    1, scatterColor, yMax, ylbl, '(d)');


%% Plot Montecarlo at the ending of the series (panel e)
axx = spMan.createAxes('mc2', h0(5), b0(5), h(5), b(5));
axes(axx); axxArray(end+1) = axx;
[xll1,yll1] = plotMonteCarlo(copulaAnalysis, monteCarloAnalysis, ...
    nWindow, scatterColor, yMax, ylbl,'(e)');

%% fix limits of panels (d-e)
xlims = ([xll;xll1]);
xlimsnew = [min(xlims(:,1)),max(xlims(:,2))];
set(axxArray(4),'xlim',xlimsnew)
set(axxArray(5),'xlim',xlimsnew)

ylims = ([yll;yll1]);
ylimsnew = [min(ylims(:,1)),max(ylims(:,2))];
set(axxArray(4),'ylim',ylimsnew)
set(axxArray(5),'ylim',ylimsnew)

%% overplot return levels
axes(axxArray(4))
jointRPPlot(retPerAnalysis, 1)
axes(axxArray(5))
jointRPPlot(retPerAnalysis ,nWindow)

end

function plotTimeSeries(timeStamps,nonStatSeries,methodology,thresholdPotNS,tMax,yMax,ylbl,...
    labelMark,xlbl,pval,pvalStat,fontSize,scatterColor)

hold on;box on;
plot(datetime(datevec(timeStamps)),nonStatSeries)
if strcmpi(methodology,'gpd')
    plot(datetime(datevec(timeStamps)),thresholdPotNS,'--r','LineWidth',2)
end
cmap = jet(256);
minC = min(scatterColor);
maxC = max(scatterColor);
colorIdx = round(1 + (scatterColor - minC) / (maxC - minC) * (size(cmap,1)-1));
colorIdx = max(1, min(colorIdx, size(cmap,1))); % Clamp values to valid index range

scatter(datetime(datevec(tMax)),yMax,[], cmap(colorIdx,:),'filled');
ylabel(ylbl,'FontSize',fontSize);

xlabel(xlbl,'FontSize',fontSize)

text(0.05,0.9,labelMark,'Units','normalized')

text(0.15,0.8,['{\it p-value}_{nonStat}= ',sprintf('%0.3g',pval)],'units','normalized')
text(0.15,0.9,['{\it p-value}_{Stat}= ',sprintf('%0.3g',pvalStat)],'units','normalized')

grid on;
end

function [xll,yll] = plotMonteCarlo(copulaAnalysis, monteCarloAnalysis, ...
    jx, scatterColor, yMax, ylbl, labelMark)

copulaFamily=copulaAnalysis.copulaParam.family;
couplingParam=copulaAnalysis.copulaParam.rho;
if isfield(copulaAnalysis.copulaParam,'rhoMean')
    couplingParamMean=copulaAnalysis.copulaParam.rhoMean;
else
    couplingParamMean=nan;
end
if copulaAnalysis.timeVaryingCopula

    yMaxLevel=copulaAnalysis.jointExtremes;
    timeStampsByTimeWindow=copulaAnalysis.copulaParam.timeStampsByTimeWindow;

else
    yMaxLevel={copulaAnalysis.jointExtremes};
    timeStampsByTimeWindow=cellfun(@(x) x{2}.timeStamps,copulaAnalysis.marginalAnalysis,'UniformOutput',0);

end

monteCarloRsmpl = monteCarloAnalysis.monteCarloRsmpl;

cmap = jet(256);
minC = min(scatterColor);
maxC = max(scatterColor);
colorIdx = round(1 + (scatterColor - minC) / (maxC - minC) * (size(cmap,1)-1));
colorIdx = max(1, min(colorIdx, size(cmap,1))); % Clamp values to valid index range

[~,Locb]=ismember(yMaxLevel{jx},yMax,'rows');

scatterMontCarl01=scatter(monteCarloRsmpl{jx}(:,1), monteCarloRsmpl{jx}(:,2));
set(scatterMontCarl01,'LineWidth',1,'Marker','o','MarkerEdgeColor',[0.5,0.5,0.5],...
    'MarkerFaceColor',[0.65,0.65,0.65],...
    'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.6,'HandleVisibility', 'off')
hold on;
scatter(yMaxLevel{jx}(:,1), yMaxLevel{jx}(:,2),[],colorIdx(Locb,:),'filled','HandleVisibility', 'off');

if strcmpi(copulaFamily,'Gaussian')
    par01=round(couplingParam{jx}(2)*100)/100;

elseif strcmpi(copulaFamily,'Frank') || strcmpi(copulaFamily,'Clayton') ||strcmpi(copulaFamily,'Gumbel')
    if isfield(copulaAnalysis.copulaParam,'rhoMean')

        par01=round(couplingParamMean{jx}(1)*100)/100;

    else
        par01=round(couplingParam{jx}(2)*100)/100;

    end
end

t1x=datestr(timeStampsByTimeWindow{jx}(1),'yyyy');
t2x=datestr(timeStampsByTimeWindow{jx}(end),'yyyy');

if strcmpi(copulaFamily,'Gaussian')
    ht=title({[t1x,' - ',t2x];['Gaussian (\rho = ',num2str(par01),')']});

elseif strcmpi(copulaFamily,'Frank') || strcmpi(copulaFamily,'Clayton') ||strcmpi(copulaFamily,'Gumbel')
    ht=title({[t1x,' - ',t2x];[char(copulaFamily),' (\theta = ',num2str(par01),')']});

end

ht.VerticalAlignment='top';


xlabel(ylbl{1});
ylabel(ylbl{2});

% set(axxArray(5),'FontSize',fontSize)
% set(axxArray(4),'FontSize',fontSize)
grid on;
text(0.05,0.9,labelMark,'Units','normalized')
% text(axxArray(jx(2)),0.05,0.9,labelMark(5),'Units','normalized')
xll=xlim;
yll=ylim;

end

function [] = gofPlot(copulaAnalysis,gofStatistics)

copulaFamily=copulaAnalysis.copulaParam.family;
couplingParam=copulaAnalysis.copulaParam.rho;
couplingParamRaw=copulaAnalysis.copulaParam.rhoRaw; % used to compute the Mann-Kendall test
snSample=gofStatistics.snSample;
corrKendallSampleDelta=gofStatistics.corrKendallSampleDelta;
corrSpearmanSampleDelta=gofStatistics.corrSpearmanSampleDelta;

if isfield(copulaAnalysis.copulaParam,'rhoMean')

    couplingParamMean=copulaAnalysis.copulaParam.rhoMean;

end

if copulaAnalysis.timeVaryingCopula

    timeStampsByTimeWindow=copulaAnalysis.copulaParam.timeStampsByTimeWindow;

else

    timeStampsByTimeWindow=cellfun(@(x) x{2}.timeStamps,copulaAnalysis.marginalAnalysis,'UniformOutput',0);

end

corrSpearmanSamplex=gofStatistics.corrSpearmanSamplex;
corrSpearmanMontex=gofStatistics.corrSpearmanMontex;
ttRho=linspace(timeStampsByTimeWindow{1}(1),timeStampsByTimeWindow{end}(end),...
    length(timeStampsByTimeWindow));
ttRhom=(cellfun(@(x) mean(x),timeStampsByTimeWindow));
if strcmpi(copulaFamily,'gaussian')
    if copulaAnalysis.timeVaryingCopula==0

        rho=cellfun(@(x) x(2),couplingParam);
        x11=datetime(datevec(ttRho));
        y11cpar=[rho;rho];
        x22=[datetime(datevec(ttRho)),datetime(datevec(ttRho))];
        y22=[[corrSpearmanSamplex{:};corrSpearmanSamplex{:}],[corrSpearmanMontex{:};corrSpearmanMontex{:}]];
        [hAx,hLine1,hLine2]=plotyy(axxArray(3),x11,y11cpar,x22,y22);
        set(hAx,{'ycolor'},{'k';'k'})
        hLine1.LineWidth=1;
        hLine2(1).LineWidth=1;
        hLine2(2).LineWidth=1;
        str=lower(char(copulaFamily));
        idx=regexp([' ' str],'(?<=\s+)\S','start')-1;
        str(idx)=upper(str(idx));
        yll=ylabel(hAx(1),sprintf('\\rho_{%s}', str));


    elseif copulaAnalysis.timeVaryingCopula==1
        rho=cellfun(@(x) x(2),couplingParam);
        rhoRaw=cellfun(@(x) x(2),couplingParamRaw);
        x11=datetime(datevec(ttRho));
        y11cpar=rho;
        x22=[datetime(datevec(ttRho)),datetime(datevec(ttRho))];
        y22=[cell2mat(corrSpearmanSamplex)', cell2mat(corrSpearmanMontex)'];
        [hAx,hLine1,hLine2]=plotyy(x11,y11cpar,x22,y22);

        if isfield(copulaAnalysis.copulaParam,'rhoMean')
            set(hAx(1),'NextPlot','add')
            plot(hAx(1),datetime(datevec(ttRho)),cell2mat(copulaAnalysis.copulaParam.rhoMean),...
                '--','LineWidth',1)
        end
        set(hAx,{'ycolor'},{'k';'k'})
        hLine1.LineWidth=1;
        hLine2(1).LineWidth=1;
        hLine2(2).LineWidth=1;
        [~,p_value]=tsMann_Kendall(rhoRaw,0.05);

        str=lower(char(copulaFamily));
        idx=regexp([' ' str],'(?<=\s+)\S','start')-1;
        str(idx)=upper(str(idx));
        yll=ylabel(hAx(1),sprintf('\\rho_{%s}', str));
        grid on;

        text(0.65,0.1,['{\it p-value}_{MK}= ',sprintf('%0.3g',p_value)],'units','normalized',...
            'HorizontalAlignment','left')

        xlabel('Date (time)')

        formattedValue = sprintf('%.2g', corrSpearmanSampleDelta); %
        formattedValue2 = sprintf('%.1g', snSample); %
        formattedValue3 = sprintf('%.2g', corrKendallSampleDelta); %
        if copulaAnalysis.timeVaryingCopula
            latexString = sprintf('$\\overline{\\Delta \\rho}_{\\mathrm{Spearman}} = %s$', formattedValue);
            latexString2 = sprintf('$\\overline{S_n} = %s$', formattedValue2);
            latexString3 = sprintf('$\\overline{\\Delta \\tau}_{\\mathrm{Kendall}} = %s$', formattedValue3);
        else
            latexString = sprintf('${\\Delta \\rho}_{\\mathrm{Spearman}} = %s$', formattedValue);
            latexString2 = sprintf('${S_n} = %s$', formattedValue2);
            latexString3 = sprintf('${\\Delta \\tau}_{\\mathrm{Kendall}} = %s$', formattedValue3);
        end
        text(0.65, 0.4, latexString,'units','normalized','HorizontalAlignment','left','Interpreter',...
            'latex');
        text(0.65, 0.3, latexString3,'units','normalized','HorizontalAlignment','left','Interpreter',...
            'latex');
        text(0.65, 0.2, latexString2 ,'units','normalized','HorizontalAlignment','left','Interpreter',...
            'latex');

        % set(axxArray(3),'FontSize',fontSize)

        text(0.05,0.9,'(c)','Units','normalized')

        legend1 = legend('show');
        if isfield(copulaAnalysis.copulaParam,'rhoMean')
            set(legend1,...
                'Position',[0.0647530040053402 0.237293699277636 0.140854472630173 0.0767690253671562],...
                'EdgeColor','none',...
                'Color','none',...
                'AutoUpdate','off','String',{['$',yll.String,'$'],'$\overline{\theta}_{Gumbel}$',...
                ['$','\rho_{Spearman, S}','$'],['$','\rho_{Spearman, MC}','$']},'Interpreter','latex');
        else
            set(legend1,...
                'Position',[0.0647530040053402 0.237293699277636 0.140854472630173 0.0767690253671562],...
                'EdgeColor','none',...
                'Color','none',...
                'AutoUpdate','off','String',{yll.String,"\rho_{Spearman, S}","\rho_{Spearman, MC}"});
        end

    end
elseif strcmpi(copulaFamily,'clayton') || strcmpi(copulaFamily,'gumbel') || strcmpi(copulaFamily,'frank')
    if copulaAnalysis.timeVaryingCopula==0

        x11=datetime(datevec(ttRho));
        y11cpar=[cell2mat(couplingParam);cell2mat(couplingParam)];
        x22=[datetime(datevec(ttRho)),datetime(datevec(ttRho))];
        y22=[[corrSpearmanSamplex{:};corrSpearmanSamplex{:}],[corrSpearmanMontex{:};corrSpearmanMontex{:}]];
        [hAx,hLine1,hLine2]=plotyy(axxArray(3),x11,y11cpar,x22,y22);
        set(hAx,{'ycolor'},{'k';'k'})
        hLine1.LineWidth=1;
        hLine2(1).LineWidth=1;
        hLine2(2).LineWidth=1;
        str=lower(char(copulaFamily));
        idx=regexp([' ' str],'(?<=\s+)\S','start')-1;
        str(idx)=upper(str(idx));
        yll=ylabel(hAx(1),sprintf('\\theta_{%s}', str));

    elseif copulaAnalysis.timeVaryingCopula==1
        x11=datetime(datevec(ttRho));
        y11cpar=cellfun(@(cpar) cpar(1,2), couplingParam);
        x22=[datetime(datevec(ttRho)),datetime(datevec(ttRho))];
        y22=[cell2mat(corrSpearmanSamplex)', cell2mat(corrSpearmanMontex)'];
        [hAx,hLine1,hLine2]=plotyy(x11,y11cpar,x22,y22);

        if isfield(copulaAnalysis.copulaParam,'rhoMean')
            set(hAx(1),'NextPlot','add')
            plot(hAx(1),datetime(datevec(ttRho)),cell2mat(couplingParamMean),'--','LineWidth',1)
        end
        set(hAx,{'ycolor'},{'k';'k'})
        hLine1.LineWidth=1;
        hLine2(1).LineWidth=1;
        hLine2(2).LineWidth=1;
        [~,p_value]=tsMann_Kendall(y11cpar,0.05);
        %ttRhom

        significance_value_tau = 0.05;
        significance_value_ac = 0.05;
        [~, ~, p_value, H] = tsModified_MannKendall_test(ttRho, y11cpar, significance_value_tau, significance_value_ac);

        str=lower(char(copulaFamily));
        idx=regexp([' ' str],'(?<=\s+)\S','start')-1;
        str(idx)=upper(str(idx));
        yll=ylabel(hAx(1),sprintf('\\theta_{%s}', str));
        yll2=ylabel(hAx(2),sprintf('\\rho_{%s}', 'Spearman'));
        yll2P=yll2.Position;
        yll2P(1)=yll2P(1)-0.03*yll2P(1);
        set(yll2,'Position',yll2P)
        grid on;
        if copulaAnalysis.timeVaryingCopula
            text(0.45,0.7,['{\it p-value}_{\theta}= ',sprintf('%0.3g',p_value)],'units','normalized',...
                'HorizontalAlignment','left')
        end
        xlabel('Date (time)')

        formattedValue = sprintf('%.2g', corrSpearmanSampleDelta); %
        formattedValue2 = sprintf('%.1g', snSample); %
        formattedValue3 = sprintf('%.2g', corrKendallSampleDelta); %
        if copulaAnalysis.timeVaryingCopula
            latexString = sprintf('$\\overline{\\Delta \\rho}_{\\mathrm{Spearman}} = %s$', formattedValue);
            latexString2 = sprintf('$\\overline{S_n} = %s$', formattedValue2);
            latexString3 = sprintf('$\\overline{\\Delta \\tau}_{\\mathrm{Kendall}} = %s$', formattedValue3);
        else
            latexString = sprintf('${\\Delta \\rho}_{\\mathrm{Spearman}} = %s$', formattedValue);
            latexString2 = sprintf('${S_n} = %s$', formattedValue2);
            latexString3 = sprintf('${\\Delta \\tau}_{\\mathrm{Kendall}} = %s$', formattedValue3);
        end
        t1=text(0.65, 0.35, latexString,'units','normalized','HorizontalAlignment','left','Interpreter',...
            'latex','FontSize',12);
        t2=text(0.65, 0.25, latexString3,'units','normalized','HorizontalAlignment','left','Interpreter',...
            'latex','FontSize',12);
        t3=text(0.65, 0.15, latexString2 ,'units','normalized','HorizontalAlignment','left','Interpreter',...
            'latex','FontSize',12);
        texts = [t1 t2 t3];
        drawnow;

        % --- combine text extents (still in axes-normalized units) ---
        extents = cell2mat(get(texts,'Extent'));
        xMin = min(extents(:,1));
        yMin = min(extents(:,2));
        xMax = max(extents(:,1)+extents(:,3));
        yMax = max(extents(:,2)+extents(:,4));
        pad = 0.01;
        posAxNorm = [xMin-pad, yMin-pad, (xMax-xMin)+2*pad, (yMax-yMin)+2*pad];

        % --- convert axes-normalized -> figure-normalized coordinates ---
        ax = ancestor(t1,'axes');
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

        % --- draw annotation rectangle in the right place ---
        annotation('rectangle', posFigNorm, ...
            'EdgeColor','k','LineWidth',1.5,'FaceColor','none');

        % set(axxArray(3),'FontSize',fontSize)

        text(0.05,0.9,'(c)','Units','normalized')

        legend1 = legend('show');
        if isfield(copulaAnalysis.copulaParam,'rhoMean')
            set(legend1,...
                'Position',[0.0647530040053402 0.237293699277636 0.140854472630173 0.0767690253671562],...
                'EdgeColor','none',...
                'Color','none',...
                'AutoUpdate','off','String',{['$',yll.String,'$'],'$\overline{\theta}_{Gumbel}$',...
                ['$','\rho_{Spearman, S}','$'],['$','\rho_{Spearman, MC}','$']},'Interpreter','latex');
        else
            set(legend1,...
                'Position',[0.0647530040053402 0.237293699277636 0.140854472630173 0.0767690253671562],...
                'EdgeColor','none',...
                'Color','none',...
                'AutoUpdate','off','String',{yll.String,"\rho_{Spearman, S}","\rho_{Spearman, MC}"});
        end
    end
end
end

function jointRPPlot(rpPlot,jx)
colorChars = 'rgbcmykw';

rpp=[rpPlot.jointAndRP];

cRL=colorChars(1:length(rpp));

for ij=1:length(rpp)
    IUUcell=rpPlot(ij).X;
    IVVCell=rpPlot(ij).Y;
    % construct x, y
    X=cellfun(@(x1) [x1';x1'],IUUcell,'UniformOutput',0);
    Y=cellfun(@(x1) [x1';x1'],IVVCell,'UniformOutput',0);

    [xd,ixd]=unique(X{jx}(:));
    yd=Y{jx}(:);
    yd=yd(ixd);
    [yd,ixd]=unique(yd);
    xd=xd(ixd);

    % Define the new area (xrange and yrange) for extrapolation

    x_range = xlim; % New x-limits for extrapolation
    x_extended=linspace(x_range(1),x_range(end),100);
    y_extended=interp1(xd,yd,x_extended,'linear','extrap');
    plot(x_extended,y_extended,cRL(ij),'LineWidth',2,'DisplayName',[num2str(rpp(ij)),[' - year R.P.']])
    hold on;

end
legend show

end
