function [axxArray] = tsCopulaTimeVaryingPlot(copulaAnalysis,gofStatistics,varargin)

%tsCopulaTimeVaryingPlot plotting joint peaks and Monte-Carlo resampled
%values
% [handles] = tsCopulaTimeVaryingPlot(copulaAnalysis, varargin)
%                     two plots are generated based on the copulaAnalysis
% structure file. First plot shows time series of input series along with
% joint peaks marked with asterisk. In the second plot, a scatter plot is
% generated using the joint peaks and the Monte-Carlo resampled values.

% input:
%  copulaAnalysis          - a variable of type structure provided as the
%                            output of tsCopulaCompoundGPD or
%                            tsCopulaCompoundGPDMontecarlo functions
%  gofStatistics           - a variable of type structure provided as the
%                             output of tsCopulaUncertainty function


% output:
%  handles:                - a set of handles concerning the generated plots;
%                           first entry is with regard to the time series plot
%                           and second entry is with regard to the scatter plots
%
%
%
%

% M.H.Bahmanpour, 2024

% setting the default parameters

args.xlbl = 'Date (time)';
args.ylbl = 'Y';
args.fontSize = 12;

args.panels=[];
args.scale=1;

args = tsEasyParseNamedArgs(varargin, args);

xlbl = args.xlbl;
ylbl = args.ylbl;
fontSize = args.fontSize;

% set some parameters
labelMark=(["(a)","(b)","(c)","(d)","(e)"]);

%read some parameters

methodology=copulaAnalysis.methodology;

inputtimestampsWindowCell=copulaAnalysis.copulaParam.inputtimestampsWindowCell;


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

c=(copulaAnalysis.jointExtremes);
ct=(copulaAnalysis.jointExtremeTimeStamps);
ct2=vertcat(ct{:});
[yMax,iB,~] = unique(vertcat(c{:}),'stable','rows'); 
tMax=ct2(iB,:);
% yMax=copulaAnalysisStat.jointExtremes;              %non-stationary peaks
% tMax=copulaAnalysisStat.jointExtremeTimeStamps;     %time of occurrences of non-stationary peaks
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


%%% plotting scatters of samples (joint peaks) co-plotted with Monte Carlo-generated samples using the
%%% chosen copula model; Due to the time-varying nature of the coupling, plotting is performed both for
%%% the initial and ending of the time series using a pre-defined time window

%read relevant input data
resampleLevel=copulaAnalysis.resampleLevel;
yMaxLevel=copulaAnalysis.jointExtremes;
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
grid(axxArray(1),'on')
grid(axxArray(2),'on')
text(axxArray(4),0.05,0.9,labelMark(4),'Units','normalized')
text(axxArray(5),0.05,0.9,labelMark(5),'Units','normalized')

%annotate correlations
corrSp01M=corr(resampleLevel{1},'type','spearman');
corrSp01M=round(corrSp01M(2)*100)/100;
corrSp01S=corr(yMaxLevel{1},'type','spearman');
corrSp01S=round(corrSp01S(2)*100)/100;
corrSp02M=corr(resampleLevel{end},'type','spearman');
corrSp02M=round(corrSp02M(2)*100)/100;
corrSp02S=corr(yMaxLevel{end},'type','spearman');
corrSp02S=round(corrSp02S(2)*100)/100;

% text(axxArray(4),0.6,0.08,sprintf('\\rho_{Spearman, S} =  %2g\n\\rho_{Spearman, MC} =  %2g', corrSp01S,corrSp01M),'Units','normalized')
% text(axxArray(5),0.6,0.08,sprintf('\\rho_{Spearman, S} =  %2g\n\\rho_{Spearman, MC} =  %2g', corrSp02S,corrSp02M),'Units','normalized')

%%% calculation of p-value (using Mann-Kendall test) for time-varying copula parameter
% ttRho=cellfun(@(x) x(round(length(x)/2)),inputtimestampsWindowCell);

ttRho=linspace(inputtimestampsWindowCell{1}(1),inputtimestampsWindowCell{end}(end),length(inputtimestampsWindowCell));
positionLabel2=[];
yLabel2=[];
if strcmpi(copulaFamily,'gaussian')
    rho=cellfun(@(x) x(2),couplingParam);
    plot(axxArray(3),datetime(datevec(ttRho)),(rho),'LineWidth',1)
    [H,p_value]=tsMann_Kendall((rho),0.05);

    str=lower(char(copulaFamily));
    idx=regexp([' ' str],'(?<=\s+)\S','start')-1;
    str(idx)=upper(str(idx));
    yll=ylabel(axxArray(3),sprintf('\\rho_{%s}', str));
    % ylabel(axxArray(3),'\rho_{Gaussian}')
    positionLabel2=[positionLabel2;yll.Position];
    yLabel2=[yLabel2,yll];
elseif strcmpi(copulaFamily,'clayton') || strcmpi(copulaFamily,'gumbel') || strcmpi(copulaFamily,'frank')
    corrSpearmanSamplex=gofStatistics.corrSpearmanSamplex;
    corrSpearmanMontex=gofStatistics.corrSpearmanMontex;
    % [hAx,hLine1,hLine2] = plotyy(x,y1,[x',x'],[y2',y3']);
    % plot(axxArray(3),datetime(datevec(ttRho)),cell2mat(couplingParam),'LineWidth',1)
     [hAx,hLine1,hLine2]=plotyy(axxArray(3),datetime(datevec(ttRho)),cell2mat(couplingParam),[datetime(datevec(ttRho)),datetime(datevec(ttRho))],[[corrSpearmanSamplex{:}]',[corrSpearmanMontex{:}]']);
    if isfield(copulaAnalysis.copulaParam,'rhoMean')
        set(hAx(1),'NextPlot','add')
        plot(hAx(1),datetime(datevec(ttRho)),cell2mat(copulaAnalysis.copulaParam.rhoMean),'--','LineWidth',1)
    end
     set(hAx,{'ycolor'},{'k';'k'}) 
     hLine1.LineWidth=1;
     hLine2(1).LineWidth=1;
     hLine2(2).LineWidth=1;
    [H,p_value]=tsMann_Kendall(cell2mat(couplingParam),0.05);


    str=lower(char(copulaFamily));
    idx=regexp([' ' str],'(?<=\s+)\S','start')-1;
    str(idx)=upper(str(idx));
    yll=ylabel(hAx(1),sprintf('\\theta_{%s}', str));
    % yll2=ylabel(hAx(2),'\rho_{Spearman}');
    positionLabel2=[positionLabel2;yll.Position];
    yLabel2=[yLabel2,yll];
end
grid(axxArray(3),'on')
text(axxArray(3),0.65,0.1,['{\it p-value}_{MK}= ',sprintf('%0.3g',p_value)],'units','normalized','HorizontalAlignment','left')
xlabel(axxArray(3),'Date (time)')

pval=cellfun(@(x) x{2}.pValueChange,copulaAnalysis.marginalAnalysis);
text(axxArray(1),0.5,0.1,['{\it p-value}_{MK}= ',sprintf('%0.3g',pval(1))],'units','normalized')
text(axxArray(2),0.5,0.1,['{\it p-value}_{MK}= ',sprintf('%0.3g',pval(2))],'units','normalized')

indexToCopula=find(strcmp(copulaFamily,{gofStatistics.copulaFamily}));
snSample=gofStatistics(indexToCopula).snSample;
corrKendallSampleDelta=gofStatistics(indexToCopula).corrKendallSampleDelta;
corrSpearmanSampleDelta=gofStatistics(indexToCopula).corrSpearmanSampleDelta;

% text(axxArray(3),0.65,0.2,['Sn = ',sprintf('%0.1g',snSample)],'units','normalized','HorizontalAlignment','left')
% text(axxArray(3),0.65,0.3,['\Delta\tau_{Kendall} = ',sprintf('%0.1g',corrKendallSampleDelta)],'units','normalized','HorizontalAlignment','left')
% text(axxArray(3),0.65,0.4,['\Delta\rho_{Spearman} = ',sprintf('%0.1g',corrSpearmanSampleDelta)],'units','normalized','HorizontalAlignment','left')
formattedValue = sprintf('%.2g', corrSpearmanSampleDelta); % 
formattedValue2 = sprintf('%.1g', snSample); % 
formattedValue3 = sprintf('%.2g', corrKendallSampleDelta); % 

latexString = sprintf('$\\overline{\\Delta \\rho}_{\\mathrm{Spearman}} = %s$', formattedValue);
% latexString2 = sprintf('$\\overline{S_n} = %s$', formattedValue2);
latexString2 = sprintf('$\\overline{S_n} = %s$', formattedValue2);
latexString3 = sprintf('$\\overline{\\Delta \\tau}_{\\mathrm{Kendall}} = %s$', formattedValue3);
text(axxArray(3),0.65, 0.4, latexString,'units','normalized','HorizontalAlignment','left','Interpreter', 'latex', 'FontSize', 14);
text(axxArray(3),0.65, 0.3, latexString3,'units','normalized','HorizontalAlignment','left','Interpreter', 'latex','FontSize', 14);
text(axxArray(3),0.65, 0.2, latexString2 ,'units','normalized','HorizontalAlignment','left','Interpreter', 'latex','FontSize', 14);

set(axxArray(1),'FontSize',fontSize)
set(axxArray(2),'FontSize',fontSize)
set(axxArray(3),'FontSize',fontSize)

text(axxArray(3),0.05,0.9,labelMark(3),'Units','normalized')
minLabel=min([positionLabel;positionLabel2]);


% string(['\bar','{',yll.String,'}'])
legend1 = legend(axxArray(3),'show');
if isfield(copulaAnalysis.copulaParam,'rhoMean')
hl=set(legend1,...
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
xlims=[get(axxArray(4),'xlim');get(axxArray(5),'xlim')];
xlimsnew=[min(xlims(:,1)),max(xlims(:,2))];
set(axxArray(4),'xlim',xlimsnew)
set(axxArray(5),'xlim',xlimsnew)

ylims=[get(axxArray(4),'ylim');get(axxArray(5),'ylim')];
ylimsnew=[min(ylims(:,1)),max(ylims(:,2))];
set(axxArray(4),'ylim',ylimsnew)
set(axxArray(5),'ylim',ylimsnew)
handles=[];
end






