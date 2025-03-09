function [axxArray] = tsCopulaPlotBivariate(copulaAnalysis,gofStatistics,varargin)
% tsCopulaPlotBivariate  plotting joint peaks and Monte-Carlo resampled
%values

% [axxArray] = tsCopulaPlotBivariate(copulaAnalysis,gofStatistics,varargin)
% first two subplots show time series of input series along with joint
% peaks. The third subplot shows gof data. The last two subplots show
% scatter plots using the joint peaks and the Monte-Carlo resampled values.



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

% [1] Bahmanpour, M.H., Mentaschi, L., Tilloy, A., Vousdoukas, M.,
%     Federico, I., Coppini, G., and Feyen, L., 2025,
%     Transformed-Stationary EVA 2.0: A Generalized Framework for
%     Non-stationary Joint Extreme Analysis (submitted to Hydrology and
%     Earth System Sciences; Feb 2025)

% setting the default parameters

args.xlbl = 'Date (time)';
args.ylbl = 'Y';
args.fontSize = 12;
args.smoothInd=1;

args = tsEasyParseNamedArgs(varargin, args);

xlbl = args.xlbl;
ylbl = args.ylbl;
fontSize = args.fontSize;
smoothInd = args.smoothInd;
% set some parameters
labelMark=(["(a)","(b)","(c)","(d)","(e)"]);

%read some parameters

methodology=copulaAnalysis.methodology;

%set figure and axes properties; keep axes handles as ana array

axxArray=[];
rt=1;
b0=27;
l0=27;
spMan = tsLcSubplotManager(b0, l0, 'CellXSize', round(27*rt), 'CellYSize', round(27*rt), 'gap', [0 0]);
spMan.initFigure;
h=[7,7,7,11,11];
b=[13,13,13,11,11];
h0=[7,16,25,11.5,24.5];
b0=[2,2,2,17.5,17.5];
for ij=1:length(h0)
    % format    [top left height width]   unlike the usual MATLAB
    % format of [left bottom width height]
    axx=spMan.createAxes(num2str(ij),h0(ij),b0(ij),h(ij),b(ij));
    axxArray=[axxArray,axx];
end

%%% plot time series, joint peaks and threshold curves

% read stationary copula information for plotting time series
if copulaAnalysis.timeVaryingCopula==1

    yMax=copulaAnalysis.yMax;
    tMax=copulaAnalysis.tMax;
else
    yMax=copulaAnalysis.jointExtremes;              %non-stationary peaks
    tMax=copulaAnalysis.jointExtremeTimeStamps;     %time of occurrences of non-stationary peaks
end

thresholdPotNS=copulaAnalysis.thresholdPotNS;   %curves representing non-stationary threshold

% stationary joint peaks were sorted in the sampling process; However,
% non-stationary peaks aren't necessarily sorted; since we use
% non-stationary peaks as a way of color-coding the data, we'll have them
% sorted at this stage;

[~,iyMax]=sort(mean(yMax,2),'descend');
yMax=yMax(iyMax,:);
tMax=tMax(iyMax,:);
%
% obtain original input series data based on data contained withing copulaAnalysis

marginalAnalysis=copulaAnalysis.marginalAnalysis;
nonStatSeriesCell=cellfun(@(x) x{2}.nonStatSeries,marginalAnalysis,'UniformOutput',0);
timeStampsCell=cellfun(@(x) x{2}.timeStamps,marginalAnalysis,'UniformOutput',0);
nonStatSeries=cell2mat(nonStatSeriesCell);
timeStamps=cell2mat(timeStampsCell);

% plotting
positionLabel=[];
yLabel=[];
for jx=1:size(nonStatSeries,2)

    plot(axxArray(jx),datetime(datevec(timeStamps(:,jx))),nonStatSeries(:,jx))
    set(axxArray(jx),'nextplot','add')
    if strcmpi(methodology,'gpd')
        plot(axxArray(jx),datetime(datevec(timeStamps(:,jx))),thresholdPotNS(:,jx),'LineWidth',2)
    end
    scatterPlot=scatter(axxArray(jx),datetime(datevec(tMax(:,jx))),yMax(:,jx),[],mean(yMax,2),'filled');
    yl=ylabel(axxArray(jx),ylbl{jx});
    positionLabel=[positionLabel;yl.Position];
    yLabel=[yLabel,yl];
    xlabel(axxArray(jx),xlbl)

    text(axxArray(jx),0.05,0.9,labelMark(jx),'Units','normalized')
    grid on
end
minLabel=min(positionLabel);
yLabel(1).Position(1)=minLabel(1);
yLabel(2).Position(1)=minLabel(1);


Cdata = scatterPlot.CData;
cmap = colormap("jet");
% make it into a index image.
cmin = min(Cdata(:));
cmax = max(Cdata(:));
m = length(cmap);
index = fix((Cdata-cmin)/(cmax-cmin)*m)+1; %A
% Then to RGB
RGB = squeeze(ind2rgb(index,cmap));

grid(axxArray(1),'on')
grid(axxArray(2),'on')
pval=cellfun(@(x) x{2}.pValueChange,copulaAnalysis.marginalAnalysis);
text(axxArray(1),0.5,0.1,['{\it p-value}_{MK}= ',sprintf('%0.3g',pval(1))],'units','normalized')
text(axxArray(2),0.5,0.1,['{\it p-value}_{MK}= ',sprintf('%0.3g',pval(2))],'units','normalized')
set(axxArray(1),'FontSize',fontSize)
set(axxArray(2),'FontSize',fontSize)
%%% plotting scatters of samples (joint peaks) co-plotted with Monte Carlo-generated samples using the
%%% chosen copula model; two scatters are plotted one at the begining and
%%% one at the ending of the series
%%
%read relevant input data
resampleLevel=copulaAnalysis.resampleLevel;
if copulaAnalysis.timeVaryingCopula==1

    yMaxLevel=copulaAnalysis.jointExtremes;
elseif copulaAnalysis.timeVaryingCopula==0
    yMaxLevel={copulaAnalysis.jointExtremes};

end
[~,Locb]=ismember(yMaxLevel{1},yMax,'rows');
[~,Locb2]=ismember(yMaxLevel{end},yMax,'rows');

couplingParam=copulaAnalysis.copulaParam.rho;
copulaFamily=copulaAnalysis.copulaParam.family;
%perform plotting
scatterMontCarl01=scatter(axxArray(4),resampleLevel{1}(:,1), resampleLevel{1}(:,2));
set(scatterMontCarl01,'LineWidth',1,'Marker','o','MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor',[0.65,0.65,0.65],...
    'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.6,'HandleVisibility', 'off')
hold(axxArray(4),'on')
scatter(axxArray(4),yMaxLevel{1}(:,1), yMaxLevel{1}(:,2),[],RGB(Locb,:),'filled','HandleVisibility', 'off');

scatterMontCarl02=scatter(axxArray(5),resampleLevel{end}(:,1), resampleLevel{end}(:,2));
set(scatterMontCarl02,'LineWidth',1,'Marker','o','MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor',[0.65,0.65,0.65],...
    'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.6,'HandleVisibility', 'off')
hold(axxArray(5),'on')
scatter(axxArray(5),yMaxLevel{end}(:,1), yMaxLevel{end}(:,2),[],RGB(Locb2,:),'filled','HandleVisibility', 'off');

if strcmpi(copulaFamily,'Gaussian')
    par01=round(couplingParam{1}(2)*100)/100;
    par02=round(couplingParam{end}(2)*100)/100;
elseif strcmpi(copulaFamily,'Frank') || strcmpi(copulaFamily,'Clayton') ||strcmpi(copulaFamily,'Gumbel')
    if isfield(copulaAnalysis.copulaParam,'rhoMean')
        couplingParamMean=copulaAnalysis.copulaParam.rhoMean;
        par01=round(couplingParamMean{1}(1)*100)/100;
        par02=round(couplingParamMean{end}(1)*100)/100;
    else
        par01=round(couplingParam{1}(1)*100)/100;
        par02=round(couplingParam{end}(1)*100)/100;

    end
end
if copulaAnalysis.timeVaryingCopula==1

    inputtimestampsWindowCell=copulaAnalysis.copulaParam.inputtimestampsWindowCell;

elseif copulaAnalysis.timeVaryingCopula==0
    inputtimestampsWindowCell=cellfun(@(x) x{2}.timeStamps,copulaAnalysis.marginalAnalysis,'UniformOutput',0);

end

t1x=datestr(inputtimestampsWindowCell{1}(1),'yyyy');
t2x=datestr(inputtimestampsWindowCell{1}(end),'yyyy');
t1x2=datestr(inputtimestampsWindowCell{end}(1),'yyyy');
t2x2=datestr(inputtimestampsWindowCell{end}(end),'yyyy');
if strcmpi(copulaFamily,'Gaussian')
    ht=title(axxArray(4),{[t1x,' - ',t2x];['Gaussian (\rho = ',num2str(par01),')']});
    ht2=title(axxArray(5),{[t1x2,' - ',t2x2];['Gaussian (\rho = ',num2str(par02),')']});
elseif strcmpi(copulaFamily,'Frank') || strcmpi(copulaFamily,'Clayton') ||strcmpi(copulaFamily,'Gumbel')
    ht=title(axxArray(4),{[t1x,' - ',t2x];[char(copulaFamily),' (\theta = ',num2str(par01),')']});
    ht2=title(axxArray(5),{[t1x2,' - ',t2x2];[char(copulaFamily),' (\theta = ',num2str(par02),')']});
end

ht.VerticalAlignment='top';
ht2.VerticalAlignment='top';

xlabel(axxArray(5),ylbl{1});
ylabel(axxArray(5),ylbl{2});

xlabel(axxArray(4),ylbl{1});
ylabel(axxArray(4),ylbl{2});

set(axxArray(5),'FontSize',fontSize)
set(axxArray(4),'FontSize',fontSize)
set(axxArray(5),'xgrid','on','ygrid','on')
set(axxArray(4),'xgrid','on','ygrid','on')
text(axxArray(4),0.05,0.9,labelMark(4),'Units','normalized')
text(axxArray(5),0.05,0.9,labelMark(5),'Units','normalized')

xlims=[get(axxArray(4),'xlim');get(axxArray(5),'xlim')];
xlimsnew=[min(xlims(:,1)),max(xlims(:,2))];
set(axxArray(4),'xlim',xlimsnew)
set(axxArray(5),'xlim',xlimsnew)

ylims=[get(axxArray(4),'ylim');get(axxArray(5),'ylim')];
ylimsnew=[min(ylims(:,1)),max(ylims(:,2))];
set(axxArray(4),'ylim',ylimsnew)
set(axxArray(5),'ylim',ylimsnew)
%%

ttRho=linspace(inputtimestampsWindowCell{1}(1),inputtimestampsWindowCell{end}(end),length(inputtimestampsWindowCell));
positionLabel2=[];
yLabel2=[];
corrSpearmanSamplex=gofStatistics.corrSpearmanSamplex;
corrSpearmanMontex=gofStatistics.corrSpearmanMontex;
if strcmpi(copulaFamily,'gaussian')
    if copulaAnalysis.timeVaryingCopula==0

        rho=cellfun(@(x) x(2),couplingParam);
        x11=datetime(datevec(ttRho));
        y11=[rho;rho];
        x22=[datetime(datevec(ttRho)),datetime(datevec(ttRho))];
        y22=[[corrSpearmanSamplex{:};corrSpearmanSamplex{:}],[corrSpearmanMontex{:};corrSpearmanMontex{:}]];
        [hAx,hLine1,hLine2]=plotyy(axxArray(3),x11,y11,x22,y22);
        set(hAx,{'ycolor'},{'k';'k'})
        hLine1.LineWidth=1;
        hLine2(1).LineWidth=1;
        hLine2(2).LineWidth=1;
        str=lower(char(copulaFamily));
        idx=regexp([' ' str],'(?<=\s+)\S','start')-1;
        str(idx)=upper(str(idx));
        yll=ylabel(hAx(1),sprintf('\\rho_{%s}', str));

        positionLabel2=[positionLabel2;yll.Position];
        yLabel2=[yLabel2,yll];
    elseif copulaAnalysis.timeVaryingCopula==1
        rho=cellfun(@(x) x(2),couplingParam);
        x11=datetime(datevec(ttRho));
        y11=smooth(rho,smoothInd);
        x22=[datetime(datevec(ttRho)),datetime(datevec(ttRho))];
        y22=[smooth([corrSpearmanSamplex{:}]',smoothInd),smooth([corrSpearmanMontex{:}]',smoothInd)];
        [hAx,hLine1,hLine2]=plotyy(axxArray(3),x11,y11,x22,y22);

        if isfield(copulaAnalysis.copulaParam,'rhoMean')
            set(hAx(1),'NextPlot','add')
            plot(hAx(1),datetime(datevec(ttRho)),cell2mat(copulaAnalysis.copulaParam.rhoMean),'--','LineWidth',1)
        end
        set(hAx,{'ycolor'},{'k';'k'})
        hLine1.LineWidth=1;
        hLine2(1).LineWidth=1;
        hLine2(2).LineWidth=1;
        [~,p_value]=tsMann_Kendall(rho,0.05);

        str=lower(char(copulaFamily));
        idx=regexp([' ' str],'(?<=\s+)\S','start')-1;
        str(idx)=upper(str(idx));
        yll=ylabel(hAx(1),sprintf('\\rho_{%s}', str));

        positionLabel2=[positionLabel2;yll.Position];
        yLabel2=[yLabel2,yll];
    end
elseif strcmpi(copulaFamily,'clayton') || strcmpi(copulaFamily,'gumbel') || strcmpi(copulaFamily,'frank')
    if copulaAnalysis.timeVaryingCopula==0

        x11=datetime(datevec(ttRho));
        y11=[cell2mat(couplingParam);cell2mat(couplingParam)];
        x22=[datetime(datevec(ttRho)),datetime(datevec(ttRho))];
        y22=[[corrSpearmanSamplex{:};corrSpearmanSamplex{:}],[corrSpearmanMontex{:};corrSpearmanMontex{:}]];
        [hAx,hLine1,hLine2]=plotyy(axxArray(3),x11,y11,x22,y22);
        set(hAx,{'ycolor'},{'k';'k'})
        hLine1.LineWidth=1;
        hLine2(1).LineWidth=1;
        hLine2(2).LineWidth=1;
        str=lower(char(copulaFamily));
        idx=regexp([' ' str],'(?<=\s+)\S','start')-1;
        str(idx)=upper(str(idx));
        yll=ylabel(hAx(1),sprintf('\\theta_{%s}', str));

        positionLabel2=[positionLabel2;yll.Position];
        yLabel2=[yLabel2,yll];

    elseif copulaAnalysis.timeVaryingCopula==1
        x11=datetime(datevec(ttRho));
        y11=smooth(cell2mat(couplingParam),smoothInd);
        x22=[datetime(datevec(ttRho)),datetime(datevec(ttRho))];
        y22=[smooth([corrSpearmanSamplex{:}]',smoothInd),smooth([corrSpearmanMontex{:}]',smoothInd)];
        [hAx,hLine1,hLine2]=plotyy(axxArray(3),x11,y11,x22,y22);
        if isfield(copulaAnalysis.copulaParam,'rhoMean')
            set(hAx(1),'NextPlot','add')
            plot(hAx(1),datetime(datevec(ttRho)),cell2mat(copulaAnalysis.copulaParam.rhoMean),'--','LineWidth',1)
        end
        set(hAx,{'ycolor'},{'k';'k'})
        hLine1.LineWidth=1;
        hLine2(1).LineWidth=1;
        hLine2(2).LineWidth=1;
        [~,p_value]=tsMann_Kendall(cell2mat(couplingParam),0.05);

        str=lower(char(copulaFamily));
        idx=regexp([' ' str],'(?<=\s+)\S','start')-1;
        str(idx)=upper(str(idx));
        yll=ylabel(hAx(1),sprintf('\\theta_{%s}', str));

        positionLabel2=[positionLabel2;yll.Position];
        yLabel2=[yLabel2,yll];
    end
end
grid(axxArray(3),'on')
if copulaAnalysis.timeVaryingCopula==1
    text(axxArray(3),0.65,0.1,['{\it p-value}_{MK}= ',sprintf('%0.3g',p_value)],'units','normalized','HorizontalAlignment','left')
end
xlabel(axxArray(3),'Date (time)')

snSample=gofStatistics.snSample;
corrKendallSampleDelta=gofStatistics.corrKendallSampleDelta;
corrSpearmanSampleDelta=gofStatistics.corrSpearmanSampleDelta;

formattedValue = sprintf('%.2g', corrSpearmanSampleDelta); %
formattedValue2 = sprintf('%.1g', snSample); %
formattedValue3 = sprintf('%.2g', corrKendallSampleDelta); %
if copulaAnalysis.timeVaryingCopula==1
    latexString = sprintf('$\\overline{\\Delta \\rho}_{\\mathrm{Spearman}} = %s$', formattedValue);
    latexString2 = sprintf('$\\overline{S_n} = %s$', formattedValue2);
    latexString3 = sprintf('$\\overline{\\Delta \\tau}_{\\mathrm{Kendall}} = %s$', formattedValue3);
elseif copulaAnalysis.timeVaryingCopula==0
    latexString = sprintf('${\\Delta \\rho}_{\\mathrm{Spearman}} = %s$', formattedValue);
    latexString2 = sprintf('${S_n} = %s$', formattedValue2);
    latexString3 = sprintf('${\\Delta \\tau}_{\\mathrm{Kendall}} = %s$', formattedValue3);
end
text(axxArray(3),0.65, 0.4, latexString,'units','normalized','HorizontalAlignment','left','Interpreter', 'latex', 'FontSize', 14);
text(axxArray(3),0.65, 0.3, latexString3,'units','normalized','HorizontalAlignment','left','Interpreter', 'latex','FontSize', 14);
text(axxArray(3),0.65, 0.2, latexString2 ,'units','normalized','HorizontalAlignment','left','Interpreter', 'latex','FontSize', 14);

set(axxArray(3),'FontSize',fontSize)

text(axxArray(3),0.05,0.9,labelMark(3),'Units','normalized')
minLabel=min([positionLabel;positionLabel2]);

legend1 = legend(axxArray(3),'show');
if isfield(copulaAnalysis.copulaParam,'rhoMean')
    set(legend1,...
        'Position',[0.0647530040053402 0.237293699277636 0.140854472630173 0.0767690253671562],...
        'EdgeColor','none',...
        'Color','none',...
        'AutoUpdate','off','String',{['$',yll.String,'$'],'$\overline{\theta}_{Gumbel}$',['$','\rho_{Spearman, S}','$'],['$','\rho_{Spearman, MC}','$']},'Interpreter','latex');
else
    set(legend1,...
        'Position',[0.0647530040053402 0.237293699277636 0.140854472630173 0.0767690253671562],...
        'EdgeColor','none',...
        'Color','none',...
        'AutoUpdate','off','String',{yll.String,"\rho_{Spearman, S}","\rho_{Spearman, MC}"});
end

end






