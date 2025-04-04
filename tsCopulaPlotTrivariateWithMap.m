function [axxArray] = tsCopulaPlotTrivariateWithMap(copulaAnalysis,gofStatistics,varargin)
% tsCopulaPlotTrivariate  plotting joint peaks and Monte-Carlo resampled
%values

% [axxArray] = tsCopulaPlotTrivariate(copulaAnalysis,gofStatistics,varargin)
% used for trivariate plotting of joint extremes and time-variation of the
% coupling parameter.
% used to plot figure 3 of Bahmanpour et al. 2025, creates a map in the
% topleft panel, needs the m_map toolbox.



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
args.locString=["Loc1","Loc2","Loc3"];
args.latlon=[];

args = tsEasyParseNamedArgs(varargin, args);

xlbl = args.xlbl;
ylbl = args.ylbl;
fontSize = args.fontSize;
locString=args.locString;
latlon=args.latlon;
% set some parameters
labelMark=(["(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)"]);

%read some parameters

methodology=copulaAnalysis.methodology;

%set figure and axes properties; keep axes handles as an array

axxArray=[];
rt=1;
hf=27;
wf=27;
spMan = tsLcSubplotManager(hf, wf, 'CellXSize', round(27*rt), 'CellYSize', round(27*rt), 'gap', [0 0]);
spMan.initFigure;
h=[repmat(5,1,12)];
b=[repmat(7.5,1,12)];%
h0=[repmat([5,12,19,26],1,3)];
b0=[2,2,2,2,11.5,11.5,11.5,11.5,21,21,21,21];
for ij=1:length(h0)
    % format    [top left height width]   unlike the usual MATLAB
    % format of [left bottom width height]
    axx=spMan.createAxes(num2str(ij),h0(ij),b0(ij),h(ij),b(ij));
    if ij==1

        m_proj('azimuthal equidistant','lon',[150], ...
            'lat',[-10],'rad',30);
        m_gshhs_f('save','coastLinePacific','patch',[.7 .9 .7])
        m_usercoast('coastLinePacific','color','k')
        m_grid('linestyle','none','linewidth',2,'tickdir','out',...
            'xaxisloc','bottom','yaxisloc','left','fontsize',6);

        for ii=1:size(latlon,2)
            m_line(latlon(1,ii),latlon(2,ii),'marker','square','markersize',4,'color','r');
        end
        m_ruler([.2 .5],.8,1,'fontsize',8);
        set(gca,'position',[-0.02    0.7959    0.2611    0.1741])

        axx_=spMan.createAxes(num2str(19),2.7,6.45,3,3.75);

        m_proj('albers equal-area','longitudes',[floor(min(latlon(1,:)))-2 ceil(max(latlon(1,:)))+2], ...
            'latitudes',[floor(min(latlon(2,:)))-1 ceil(max(latlon(2,:)))+1],'rect','on');
        m_gshhs_f('save','coastlineZoomed','patch',[.7 .9 .7]);
        m_usercoast('coastlineZoomed','color','k')
        m_grid('linestyle','none','linewidth',2,'tickdir','in',...
            'xaxisloc','top','yaxisloc','left','fontsize',6);
        for ii=1:size(latlon,2)
            m_line(latlon(1,ii),latlon(2,ii),'marker','square','markersize',4,'color','r');
            m_text(latlon(1,ii),latlon(2,ii),locString(ii),'color','k','vertical','top','fontweight','bold');
        end
        m_ruler([.2 .5],.8,1,'fontsize',8);
        xll=axx_.XLim;
        yll=axx_.YLim;
        [longx,latx]=m_xy2ll(xll,yll);
        m_proj('azimuthal equidistant','lon',[150], ...
            'lat',[-10],'rad',30);
        [xx,yy]=m_ll2xy(longx,latx);
        set(axx,'NextPlot','add')
        plot(axx,[xx(1) xx(2) xx(2) xx(1) xx(1)],[yy(1) yy(1) yy(2) yy(2) yy(1)],'LineStyle','-','LineWidth',1,'Color','k')

        annotation(gcf,'line',[0.152202937249666 0.218958611481976],...
            [0.945595460614152 0.974632843791722]);
        annotation(gcf,'line',[0.154873164218959 0.218958611481976],...
            [0.934579439252336 0.881174899866489]);
    end
    axxArray=[axxArray,axx];
end

set(axxArray(1),'FontSize',fontSize)

%%% plot time series, joint peaks and threshold curves

% read stationary copula information for plotting time series
if copulaAnalysis.timeVaryingCopula==1

    yMax=copulaAnalysis.yMax;
    tMax=copulaAnalysis.tMax;
else
    yMax=copulaAnalysis.jointExtremes;              %non-stationary peaks
    tMax=copulaAnalysis.jointExtremeTimeStamps;     %time of occurrences of non-stationary peaks
end


% stationary joint peaks were sorted in the sampling process; However,
% non-stationary peaks aren't necessarily sorted; since we use
% non-stationary peaks as a way of color-coding the data, we'll have them
% sorted at this stage;

[~,iyMax]=sort(mean(yMax,2),'descend');
yMax=yMax(iyMax,:);
tMax=tMax(iyMax,:);
%
sc3=scatter3(axxArray(5),yMax(:,1),yMax(:,2),yMax(:,3),[],yMax(:,1),'filled');
view(57.5,30)
hold(axxArray(5),'on')
ylblx2=sprintf('{\\it%s}_{%s}',locString(1),ylbl{1});
xlabel(axxArray(5),ylblx2);
ylblx2=sprintf('{\\it%s}_{%s}',locString(2),ylbl{2});
ylabel(axxArray(5),ylblx2);
ylblx3=sprintf('{\\it%s}_{%s}',locString(3),ylbl{3});
zlabel(axxArray(5),ylblx3);


set(axxArray(5),'FontSize',fontSize)
text(axxArray(5),0.05,0.9,labelMark(5),'Units','normalized')
% obtain original input series data based on data contained withing copulaAnalysis

marginalAnalysis=copulaAnalysis.marginalAnalysis;
nonStatSeriesCell=cellfun(@(x) x{2}.nonStatSeries,marginalAnalysis,'UniformOutput',0);
timeStampsCell=cellfun(@(x) x{2}.timeStamps,marginalAnalysis,'UniformOutput',0);
nonStatSeries=cell2mat(nonStatSeriesCell);
timeStamps=cell2mat(timeStampsCell);

snSample=gofStatistics.snSample;
corrKendallSampleDelta=gofStatistics.corrKendallSampleDelta;
corrSpearmanSampleDelta=gofStatistics.corrSpearmanSampleDelta;

formattedValue = sprintf('%.2g', corrSpearmanSampleDelta); %
formattedValue2 = sprintf('%.1g', snSample); %
formattedValue3 = sprintf('%.2g', corrKendallSampleDelta); %

latexString = sprintf('$\\overline{\\Delta \\rho}_{\\mathrm{Spearman}} = %s$', formattedValue);
latexString2 = sprintf('$\\overline{S_n} = %s$', formattedValue2);
latexString3 = sprintf('$\\overline{\\Delta \\tau}_{\\mathrm{Kendall}} = %s$', formattedValue3);

text(axxArray(5),0.65, 0.4, latexString,'units','normalized','HorizontalAlignment','left','Interpreter', 'latex', 'FontSize', 14);
text(axxArray(5),0.65, 0.3, latexString3,'units','normalized','HorizontalAlignment','left','Interpreter', 'latex','FontSize', 14);
text(axxArray(5),0.65, 0.2, latexString2 ,'units','normalized','HorizontalAlignment','left','Interpreter', 'latex','FontSize', 14);
%%
%%% plotting scatters of samples (joint peaks) co-plotted with Monte Carlo-generated samples using the
%%% chosen copula model; Due to the time-varying nature of the coupling, plotting is performed both for
%%% the initial and ending of the time series using a pre-defined time window

%read relevant input data

yMaxLevel=copulaAnalysis.jointExtremes;
resampleLevel=copulaAnalysis.resampleLevel;

couplingParam=copulaAnalysis.copulaParam.rho;
copulaFamily=copulaAnalysis.copulaParam.family;
%perform plotting

scatterMontCarl01=scatter(axxArray(3),resampleLevel{1}(:,1), resampleLevel{1}(:,2));
set(scatterMontCarl01,'LineWidth',1,'Marker','o','MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor',[0.65,0.65,0.65],...
    'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.6)
hold(axxArray(3),'on')
sca3=scatter(axxArray(3),yMaxLevel{1}(:,1), yMaxLevel{1}(:,2),[],yMaxLevel{1}(:,1),'filled');

scatterMontCarl02=scatter(axxArray(7),resampleLevel{1}(:,1), resampleLevel{1}(:,3));
set(scatterMontCarl02,'LineWidth',1,'Marker','o','MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor',[0.65,0.65,0.65],...
    'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.6)
hold(axxArray(7),'on')
sca7=scatter(axxArray(7),yMaxLevel{1}(:,1), yMaxLevel{1}(:,3),[],yMaxLevel{1}(:,1),'filled');

scatterMontCarl03=scatter(axxArray(11),resampleLevel{1}(:,2), resampleLevel{1}(:,3));
set(scatterMontCarl03,'LineWidth',1,'Marker','o','MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor',[0.65,0.65,0.65],...
    'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.6)
hold(axxArray(11),'on')
sca11=scatter(axxArray(11),yMaxLevel{1}(:,2), yMaxLevel{1}(:,3),[],yMaxLevel{1}(:,1),'filled');


scatterMontCarl01e=scatter(axxArray(4),resampleLevel{end}(:,1), resampleLevel{end}(:,2));
set(scatterMontCarl01e,'LineWidth',1,'Marker','o','MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor',[0.65,0.65,0.65],...
    'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.6)
hold(axxArray(4),'on')
sca4=scatter(axxArray(4),yMaxLevel{end}(:,1), yMaxLevel{end}(:,2),[],yMaxLevel{end}(:,1),'filled');

scatterMontCarl02e=scatter(axxArray(8),resampleLevel{end}(:,1), resampleLevel{end}(:,3));
set(scatterMontCarl02e,'LineWidth',1,'Marker','o','MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor',[0.65,0.65,0.65],...
    'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.6)
hold(axxArray(8),'on')
sca8=scatter(axxArray(8),yMaxLevel{end}(:,1), yMaxLevel{end}(:,3),[],yMaxLevel{end}(:,1),'filled');

scatterMontCarl03e=scatter(axxArray(12),resampleLevel{end}(:,2), resampleLevel{end}(:,3));
set(scatterMontCarl03e,'LineWidth',1,'Marker','o','MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor',[0.65,0.65,0.65],...
    'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.6)
hold(axxArray(12),'on')
sca12=scatter(axxArray(12),yMaxLevel{end}(:,2), yMaxLevel{end}(:,3),[],yMaxLevel{end}(:,1),'filled');

par01=couplingParam{1}(1,2);
par02=couplingParam{1}(1,3);
par03=couplingParam{1}(2,3);

par01e=couplingParam{end}(1,2);
par02e=couplingParam{end}(1,3);
par03e=couplingParam{end}(2,3);

ylblx=sprintf('{\\it%s}_{%s}',locString(1),ylbl{1});
xlabel(axxArray(3),ylblx);
xlabel(axxArray(4),ylblx);

ylblx=sprintf('{\\it%s}_{%s}',locString(2),ylbl{2});
ylabel(axxArray(3),ylblx);
ylabel(axxArray(4),ylblx);

ylblx=sprintf('{\\it%s}_{%s}',locString(1),ylbl{1});
xlabel(axxArray(7),ylblx);
xlabel(axxArray(8),ylblx);
ylblx=sprintf('{\\it%s}_{%s}',locString(3),ylbl{3});
ylabel(axxArray(7),ylblx);
ylabel(axxArray(8),ylblx);

ylblx=sprintf('{\\it%s}_{%s}',locString(2),ylbl{2});
xlabel(axxArray(11),ylblx);
xlabel(axxArray(12),ylblx);
ylblx=sprintf('{\\it%s}_{%s}',locString(3),ylbl{3});
ylabel(axxArray(11),ylblx);
ylabel(axxArray(12),ylblx);

grid(axxArray(3),'on')
grid(axxArray(7),'on')
grid(axxArray(11),'on')

grid(axxArray(4),'on')
grid(axxArray(8),'on')
grid(axxArray(12),'on')

text(axxArray(3),0.05,0.9,labelMark(3),'Units','normalized')
text(axxArray(7),0.05,0.9,labelMark(7),'Units','normalized')
text(axxArray(11),0.05,0.9,labelMark(11),'Units','normalized')

text(axxArray(4),0.05,0.9,labelMark(4),'Units','normalized')
text(axxArray(8),0.05,0.9,labelMark(8),'Units','normalized')
text(axxArray(12),0.05,0.9,labelMark(12),'Units','normalized')

%annotate correlations

set(axxArray(3),'FontSize',fontSize)
set(axxArray(4),'FontSize',fontSize)
set(axxArray(7),'FontSize',fontSize)
set(axxArray(8),'FontSize',fontSize)
set(axxArray(11),'FontSize',fontSize)
set(axxArray(12),'FontSize',fontSize)

axis(axxArray(7),'tight')
axis(axxArray(8),'tight')
axis(axxArray(3),'tight')
axis(axxArray(4),'tight')
axis(axxArray(11),'tight')
axis(axxArray(12),'tight')

%%
%
inputtimestampsWindowCell=copulaAnalysis.copulaParam.inputtimestampsWindowCell;

thresholdPotNS=copulaAnalysis.thresholdPotNS;   %curves representing non-stationary threshold
%
% % plotting
positionLabel=[];
yLabel=[];
jx2=[2,6,10];

for jx=1:size(nonStatSeries,2)

    plot(axxArray(jx2(jx)),datetime(datevec(timeStamps(:,jx))),nonStatSeries(:,jx))
    set(axxArray(jx2(jx)),'nextplot','add')
    if strcmpi(methodology,'gpd')
        plot(axxArray(jx2(jx)),datetime(datevec(timeStamps(:,jx))),thresholdPotNS(:,jx),'LineWidth',2)
    end
    scatterPlot=scatter(axxArray(jx2(jx)),datetime(datevec(tMax(:,jx))),yMax(:,jx),[],mean(yMax,2),'filled');
    ylblx=sprintf('{\\it%s}_{%s}',locString(jx),ylbl{jx});

    yl=ylabel(axxArray(jx2(jx)),ylblx);
    positionLabel=[positionLabel;yl.Position];
    yLabel=[yLabel,yl];
    xlabel(axxArray(jx2(jx)),xlbl)

    text(axxArray(jx2(jx)),0.05,0.9,labelMark(jx2(jx)),'Units','normalized')
    grid on
end
minLabel=min(positionLabel);
yLabel(1).Position(1)=minLabel(1);
yLabel(2).Position(1)=minLabel(1);

pval=cellfun(@(x) x{2}.pValueChange,copulaAnalysis.marginalAnalysis);
text(axxArray(2),0.5,0.1,['{\it p-value}_{MK}= ',sprintf('%0.3g',pval(1))],'units','normalized')
text(axxArray(6),0.5,0.1,['{\it p-value}_{MK}= ',sprintf('%0.3g',pval(2))],'units','normalized')
text(axxArray(10),0.5,0.1,['{\it p-value}_{MK}= ',sprintf('%0.3g',pval(3))],'units','normalized')

set(axxArray(2),'FontSize',fontSize)
set(axxArray(6),'FontSize',fontSize)
set(axxArray(10),'FontSize',fontSize)

grid(axxArray(2),'on')
grid(axxArray(6),'on')

grid(axxArray(10),'on')

Cdata = scatterPlot.CData;
cmap = colormap("jet");
% make it into a index image.
cmin = min(Cdata(:));
cmax = max(Cdata(:));
m = length(cmap);
index = fix((Cdata-cmin)/(cmax-cmin)*m)+1; %A
% Then to RGB
RGB = squeeze(ind2rgb(index,cmap));
%
%
% %%% plotting scatters of samples (joint peaks) co-plotted with Monte Carlo-generated samples using the
% %%% chosen copula model; Due to the time-varying nature of the coupling, plotting is performed both for
% %%% the initial and ending of the time series using a pre-defined time window
%
%read relevant input data
resampleLevel=copulaAnalysis.resampleLevel;
yMaxLevel=copulaAnalysis.jointExtremes;
[~,Locb]=ismember(yMaxLevel{1},yMax,'rows');
[~,Locb2]=ismember(yMaxLevel{end},yMax,'rows');

couplingParam=copulaAnalysis.copulaParam.rho;
couplingParamRaw=copulaAnalysis.copulaParam.rhoRaw; % used to compute the Mann-Kendall test
copulaFamily=copulaAnalysis.copulaParam.family;

sc3.CData=RGB;
sca3.CData=RGB(Locb,:);
sca7.CData=RGB(Locb,:);
sca11.CData=RGB(Locb,:);
sca4.CData=RGB(Locb2,:);
sca8.CData=RGB(Locb2,:);
sca12.CData=RGB(Locb2,:);

t1x=datestr(inputtimestampsWindowCell{1}(1),'yyyy');
t2x=datestr(inputtimestampsWindowCell{1}(end),'yyyy');
t1x2=datestr(inputtimestampsWindowCell{end}(1),'yyyy');
t2x2=datestr(inputtimestampsWindowCell{end}(end),'yyyy');
%
ht=title(axxArray(3),{[t1x,' - ',t2x];['Gaussian (\rho = ',num2str(par01, '%.2f'),')']});
ht2=title(axxArray(4),{[t1x2,' - ',t2x2];['Gaussian (\rho = ',num2str(par01e, '%.2f'),')']});
ht3=title(axxArray(7),{[t1x,' - ',t2x];['Gaussian (\rho = ',num2str(par02, '%.2f'),')']});
ht4=title(axxArray(8),{[t1x2,' - ',t2x2];['Gaussian (\rho = ',num2str(par02e, '%.2f'),')']});
ht5=title(axxArray(11),{[t1x,' - ',t2x];['Gaussian (\rho = ',num2str(par03, '%.2f'),')']});
ht6=title(axxArray(12),{[t1x2,' - ',t2x2];['Gaussian (\rho = ',num2str(par03e, '%.2f'),')']});

ht.VerticalAlignment='top';
ht2.VerticalAlignment='top';
ht3.VerticalAlignment='top';
ht4.VerticalAlignment='top';
ht5.VerticalAlignment='top';
ht6.VerticalAlignment='top';
%
%%
ttRho=linspace(inputtimestampsWindowCell{1}(1),inputtimestampsWindowCell{end}(end),length(inputtimestampsWindowCell));
positionLabel2=[];
yLabel2=[];

couplingParamMat=(cell2mat(cellfun(@(x) x(find(tril(x,-1))),couplingParam,'UniformOutput',0)))';
couplingParamMat = arrayfun(@(col) couplingParamMat(:, col), 1:size(couplingParamMat, 2), 'UniformOutput', false);
couplingParamMat = cell2mat(couplingParamMat);  % Convert cell array back to matrix
plot(axxArray(9),datetime(datevec(ttRho)),couplingParamMat);

couplingParamMat=(cell2mat(cellfun(@(x) x(find(tril(x,-1))),couplingParamRaw,'UniformOutput',0)))';
couplingParamMat = arrayfun(@(col) couplingParamMat(:, col), 1:size(couplingParamMat, 2), 'UniformOutput', false);
couplingParamMat = cell2mat(couplingParamMat);  % Convert cell array back to matrix
[~,p_value]=tsMann_Kendall(couplingParamMat(:,1),0.05);
[~,p_value2]=tsMann_Kendall(couplingParamMat(:,2),0.05);
[~,p_value3]=tsMann_Kendall(couplingParamMat(:,3),0.05);

str=lower(char(copulaFamily));
idx=regexp([' ' str],'(?<=\s+)\S','start')-1;
str(idx)=upper(str(idx));
yll=ylabel(axxArray(9),sprintf('\\rho_{%s}', str));

positionLabel2=[positionLabel2;yll.Position];
yLabel2=[yLabel2,yll];

grid(axxArray(9),'on')
text(axxArray(9),0.35,0.1,['{\it p-value}_{MK,1-2}= ',sprintf('%0.3g',p_value)],'units','normalized','HorizontalAlignment','left')
text(axxArray(9),0.35,0.3,['{\it p-value}_{MK,1-3}= ',sprintf('%0.3g',p_value2)],'units','normalized','HorizontalAlignment','left')
text(axxArray(9),0.35,0.5,['{\it p-value}_{MK,2-3}= ',sprintf('%0.3g',p_value3)],'units','normalized','HorizontalAlignment','left')

xlabel(axxArray(9),'Date (time)')

text(axxArray(9),0.05,0.9,labelMark(9),'Units','normalized')

legend1 = legend(axxArray(9),'show');
set(legend1,...
    'String',{'\rho_{Gaussian_{1,2}}','\rho_{Gaussian_{1,3}}','\rho_{Gaussian_{2,3}}'},'EdgeColor','none','Color','none');

end












