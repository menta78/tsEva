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
args.locString="";
args.panels=[];
args.scale=1;
args.latlon=[];
args = tsEasyParseNamedArgs(varargin, args);

xlbl = args.xlbl;
ylbl = args.ylbl;
fontSize = args.fontSize;
locString=args.locString;
locString=["Loc1","Loc2","Loc3"];
latlon=args.latlon;
% set some parameters
labelMark=(["(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)"]);

%read some parameters

methodology=copulaAnalysis.methodology;

% inputtimestampsWindowCell=copulaAnalysis.copulaParam.inputtimestampsWindowCell;


%set figure and axes properties; keep axes handles as ana array

axxArray=[];
rt=1;
hf=27;%23
wf=27;
spMan = tsLcSubplotManager(hf, wf, 'CellXSize', round(27*rt), 'CellYSize', round(27*rt), 'gap', [0 0]);
spMan.initFigure;
h=[repmat(5,1,12)];%6.5
b=[repmat(7.5,1,12)];%7.5
h0=[repmat([5,12,19,26],1,3)];
b0=[2,2,2,2,11.5,11.5,11.5,11.5,21,21,21,21];
for ij=1:length(h0)
    % format    [top left height width]   unlike the usual MATLAB
    % format of [left bottom width height]
    axx=spMan.createAxes(num2str(ij),h0(ij),b0(ij),h(ij),b(ij));
    if ij==1
        % path(path,'/Users/hadi/Downloads/m_map')

       %albers equal-area
         m_proj('azimuthal equidistant','lon',[150], ...
            'lat',[-10],'rad',30); 
         m_gshhs_f('save','coastLinePacific','patch',[.7 .9 .7])  
        m_usercoast('coastLinePacific','color','k')

        m_grid('linestyle','none','linewidth',2,'tickdir','out',...
            'xaxisloc','bottom','yaxisloc','left','fontsize',6);

        % m_proj('albers equal-area','longitudes',[floor(min(latlon(1,:))) ceil(max(latlon(1,:)))], ...
        %     'latitudes',[floor(min(latlon(2,:))) ceil(max(latlon(2,:)))],'rect','on'); 
        % % m_gshhs_f('save','coastLinePortugal','patch',[.7 .9 .7])  
        % m_usercoast('coastLineMarshall','color','k')
        % 
        % m_grid('linestyle','none','linewidth',2,'tickdir','out',...
        %     'xaxisloc','top','yaxisloc','right','fontsize',6);

       diffArray=[latlon(1,:)-latlon(1,1);latlon(1,:)-latlon(1,2);latlon(1,:)-latlon(1,3)];
        diffArray2=[latlon(2,:)-latlon(2,1);latlon(2,:)-latlon(2,2);latlon(2,:)-latlon(2,3)];
        distArray=(sqrt(diffArray.^2+diffArray2.^2));
        distArray(distArray==0)=inf;
        [im,~]=min(min(distArray));
        [ir,ic]=find(distArray==im,1);
        ext1=min([latlon(:,ir),latlon(:,ic)],[],2);
        ext2=max([latlon(:,ir),latlon(:,ic)],[],2);
        for ii=1:size(latlon,2)
            m_line(latlon(1,ii),latlon(2,ii),'marker','square','markersize',4,'color','r');
           
           % m_text(latlon(1,ii),latlon(2,ii),locString(ii),'color','k','vertical','top','fontweight','bold');

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

        % latlon2=[latlon(:,ir),latlon(:,ic)];
        % locString2=[locString(:,ir),locString(:,ic)];
        for ii=1:size(latlon,2)
            m_line(latlon(1,ii),latlon(2,ii),'marker','square','markersize',4,'color','r');
            % if ii==1
                m_text(latlon(1,ii),latlon(2,ii),locString(ii),'color','k','vertical','top','fontweight','bold');
            % else
                % m_text(latlon(1,ii),latlon(2,ii),locString(ii),'color','k','vertical','bottom','fontweight','bold');
            % end

        end
        m_ruler([.2 .5],.8,1,'fontsize',8);
        xll=axx_.XLim;
        yll=axx_.YLim;
        [longx,latx]=m_xy2ll(xll,yll);
        % m_proj('albers equal-area','longitudes',[-10 0], ...
        %     'latitudes',[36 44],'rect','on');
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


%%% plot time series, joint peaks and threshold curves

% read stationary copula information for plotting time series
c=(copulaAnalysis.jointExtremes);
ct=(copulaAnalysis.jointExtremeTimeStamps);
ct2=vertcat(ct{:});
[yMax,iB,~] = unique(vertcat(c{:}),'stable','rows'); 
tMax=ct2(iB,:);
% yMax=copulaAnalysisStat.jointExtremes;              %non-stationary peaks
% tMax=copulaAnalysisStat.jointExtremeTimeStamps;     %time of occurrences of non-stationary peaks
%thresholdPotNS=copulaAnalysis.thresholdPotNS;   %curves representing non-stationary threshold

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
% resampleLevel=copulaAnalysisStat.resampleLevel;
% scatterMontCarl03d=scatter3(axxArray(5),resampleLevel{1}(:,1),resampleLevel{1}(:,2),resampleLevel{1}(:,3));
% set(scatterMontCarl03d,'LineWidth',1,'Marker','o','MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor',[0.65,0.65,0.65],...
%     'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)

ylblx2=sprintf('{\\it%s}_{%s}',locString(1),ylbl{1});
xlabel(axxArray(5),ylblx2);
ylblx2=sprintf('{\\it%s}_{%s}',locString(2),ylbl{2});
ylabel(axxArray(5),ylblx2);
ylblx3=sprintf('{\\it%s}_{%s}',locString(3),ylbl{3});
zlabel(axxArray(5),ylblx3);
% obtain original input series data based on data contained withing copulaAnalysis

marginalAnalysis=copulaAnalysis.marginalAnalysis;
nonStatSeriesCell=cellfun(@(x) x{2}.nonStatSeries,marginalAnalysis,'UniformOutput',0);
timeStampsCell=cellfun(@(x) x{2}.timeStamps,marginalAnalysis,'UniformOutput',0);
nonStatSeries=cell2mat(nonStatSeriesCell);
timeStamps=cell2mat(timeStampsCell);


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
if strcmpi(copulaFamily,'Gaussian')
    par01=round(couplingParam{1}(1,2)*100)/100;
    par02=round(couplingParam{1}(1,3)*100)/100;
    par03=round(couplingParam{1}(2,3)*100)/100;

    par01e=round(couplingParam{end}(1,2)*100)/100;
    par02e=round(couplingParam{end}(1,3)*100)/100;
    par03e=round(couplingParam{end}(2,3)*100)/100;
elseif strcmpi(copulaFamily,'Frank') || strcmpi(copulaFamily,'Clayton') ||strcmpi(copulaFamily,'Gumbel')
    % par01=round(couplingParam{1}(1)*100)/100;

end

if strcmpi(copulaFamily,'Gaussian')
    ht=title(axxArray(3),{['Gaussian (\rho = ',num2str(par01),')']});
    ht2=title(axxArray(7),{['Gaussian (\rho = ',num2str(par02),')']});
    ht3=title(axxArray(11),{['Gaussian (\rho = ',num2str(par03),')']});

    hte=title(axxArray(4),{['Gaussian (\rho = ',num2str(par01e),')']});
    ht2e=title(axxArray(8),{['Gaussian (\rho = ',num2str(par02e),')']});
    ht3e=title(axxArray(12),{['Gaussian (\rho = ',num2str(par03e),')']});
elseif strcmpi(copulaFamily,'Frank') || strcmpi(copulaFamily,'Clayton') ||strcmpi(copulaFamily,'Gumbel')
    % ht=title(axxArray(2),{[char(copulaFamily),' (\theta = ',num2str(par01),')']});
    % ht2=title(axxArray(4),{[char(copulaFamily),' (\theta = ',num2str(par01),')']});
    % ht3=title(axxArray(6),{[char(copulaFamily),' (\theta = ',num2str(par01),')']});
end

ht.VerticalAlignment='top';
ht2.VerticalAlignment='top';
ht3.VerticalAlignment='top';

hte.VerticalAlignment='top';
ht2e.VerticalAlignment='top';
ht3e.VerticalAlignment='top';
% labelMarks2=["Loc1","Loc2","Loc3","Loc4","Loc5","Loc6","Loc7"];
% ylblx=cellfun(@(x,y) strtrim(x(1:y-1)),ylbl,regexp(ylbl,'('),'UniformOutput',0)
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
corrSp01M=corr(resampleLevel{1},'type','spearman');
corrSp01M=round(corrSp01M*100)/100;
corrSp01S=corr(yMaxLevel{1},'type','spearman');
corrSp01S=round(corrSp01S*100)/100;

corrSp01Me=corr(resampleLevel{end},'type','spearman');
corrSp01Me=round(corrSp01Me*100)/100;
corrSp01Se=corr(yMaxLevel{end},'type','spearman');
corrSp01Se=round(corrSp01Se*100)/100;


% text(axxArray(2),0.5,0.12,sprintf('\\rho_{Spearman, S} =  %2g\n\\rho_{Spearman, MC} =  %2g', corrSp01S(1,2),corrSp01M(1,2)),'Units','normalized')
% text(axxArray(5),0.5,0.12,sprintf('\\rho_{Spearman, S} =  %2g\n\\rho_{Spearman, MC} =  %2g', corrSp01S(1,3),corrSp01M(1,3)),'Units','normalized')
% text(axxArray(8),0.5,0.12,sprintf('\\rho_{Spearman, S} =  %2g\n\\rho_{Spearman, MC} =  %2g', corrSp01S(2,3),corrSp01M(2,3)),'Units','normalized')
% 
% text(axxArray(3),0.5,0.12,sprintf('\\rho_{Spearman, S} =  %2g\n\\rho_{Spearman, MC} =  %2g', corrSp01Se(1,2),corrSp01Me(1,2)),'Units','normalized')
% text(axxArray(6),0.5,0.12,sprintf('\\rho_{Spearman, S} =  %2g\n\\rho_{Spearman, MC} =  %2g', corrSp01Se(1,3),corrSp01Me(1,3)),'Units','normalized')
% text(axxArray(9),0.5,0.12,sprintf('\\rho_{Spearman, S} =  %2g\n\\rho_{Spearman, MC} =  %2g', corrSp01Se(2,3),corrSp01Me(2,3)),'Units','normalized')


indexToCopula=find(strcmp(copulaFamily,{gofStatistics.copulaFamily}));
snSample=gofStatistics(indexToCopula).snSample;
corrKendallSampleDelta=gofStatistics(indexToCopula).corrKendallSampleDelta;
corrSpearmanSampleDelta=gofStatistics(indexToCopula).corrSpearmanSampleDelta;

% text(axxArray(9),0.1,0.2,['Sn = ',sprintf('%0.3g',snSample)],'units','normalized','HorizontalAlignment','left','FontSize',14)
% text(axxArray(9),0.1,0.5,['\Delta\tau_{Kendall} = ',sprintf('%0.3g',corrKendallSampleDelta)],'units','normalized','HorizontalAlignment','left','FontSize',14)
% text(axxArray(9),0.1,0.8,['\Delta\rho_{Spearman} = ',sprintf('%0.3g',corrSpearmanSampleDelta)],'units','normalized','HorizontalAlignment','left','FontSize',14)
% set(axxArray(9),'xtick',[])
% set(axxArray(9),'xticklabel',[])
% set(axxArray(9),'ytick',[])
% set(axxArray(9),'yticklabel',[])

set(axxArray(3),'FontSize',fontSize)
set(axxArray(4),'FontSize',fontSize)
 set(axxArray(5),'FontSize',fontSize)
set(axxArray(7),'FontSize',fontSize)
set(axxArray(8),'FontSize',fontSize)
% set(axxArray(9),'FontSize',fontSize)
set(axxArray(11),'FontSize',fontSize)
set(axxArray(12),'FontSize',fontSize)

% text(axxArray(9),0.05,0.9,labelMark(7),'Units','normalized')
 text(axxArray(5),0.05,0.9,labelMark(5),'Units','normalized')

handles=[];





% 
% args.xlbl = 'Date (time)';
% args.ylbl = 'Y';
% args.fontSize = 12;
% 
% args.panels=[];
% args.scale=1;
% 
% args = tsEasyParseNamedArgs(varargin, args);
% 
% xlbl = args.xlbl;
% ylbl = args.ylbl;
% fontSize = args.fontSize;
% 
% % set some parameters
% labelMark=(["(a)","(b)","(c)","(d)","(e)"]);
% 
% %read some parameters
% 
% methodology=copulaAnalysis.methodology;
% 
 inputtimestampsWindowCell=copulaAnalysis.copulaParam.inputtimestampsWindowCell;
% 
% 
% %set figure and axes properties; keep axes handles as ana array
% 
% axxArray=[];
% rt=1;
% b0=27;
% l0=27;
% spMan = tsLcSubplotManager(b0, l0, 'CellXSize', round(27*rt), 'CellYSize', round(27*rt), 'gap', [0 0]);
% spMan.initFigure;
% h=[7,7,7,11,11];
% b=[13,13,13,11,11];
% h0=[7,16,25,11.5,24.5];
% b0=[2,2,2,17.5,17.5];
% for ij=1:length(h0)
%     % format    [top left height width]   unlike the usual MATLAB
%     % format of [left bottom width height]
%     axx=spMan.createAxes(num2str(ij),h0(ij),b0(ij),h(ij),b(ij));
%     axxArray=[axxArray,axx];
% end
% 
% 
% 
% %%% plot time series, joint peaks and threshold curves
% 
% % read stationary copula information for plotting time series
% 
%  yMax_=copulaAnalysisStat.jointExtremes;              %non-stationary peaks
% tMax=copulaAnalysisStat.jointExtremeTimeStamps;     %time of occurrences of non-stationary peaks
 thresholdPotNS=copulaAnalysis.thresholdPotNS;   %curves representing non-stationary threshold
% 
% % stationary joint peaks were sorted in the sampling process; However,
% % non-stationary peaks aren't necessarily sorted; since we use
% % non-stationary peaks as a way of color-coding the data, we'll have them
% % sorted at this stage;
% 
 [~,iyMax]=sort(mean(yMax,2),'descend');
yMax=yMax(iyMax,:);
 tMax=tMax(iyMax,:);
% %
% % obtain original input series data based on data contained withing copulaAnalysis
% 


% 
% % plotting
 positionLabel=[];
 yLabel=[];
 jx2=[2,6,10];
 mean(yMax,2)
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
copulaFamily=copulaAnalysis.copulaParam.family;
%perform plotting
% scatterMontCarl01=scatter(axxArray(4),resampleLevel{1}(:,1), resampleLevel{1}(:,2));
% set(scatterMontCarl01,'LineWidth',1,'Marker','o','MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor',[0.65,0.65,0.65],...
%     'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.6)
% hold(axxArray(4),'on')
sca3.CData=RGB(Locb,:);
sc3.CData=RGB;
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
if strcmpi(copulaFamily,'Gaussian')
    ht=title(axxArray(3),{[t1x,' - ',t2x];['Gaussian (\rho = ',num2str(par01),')']});
    ht2=title(axxArray(4),{[t1x2,' - ',t2x2];['Gaussian (\rho = ',num2str(par01e),')']});
    ht3=title(axxArray(7),{[t1x,' - ',t2x];['Gaussian (\rho = ',num2str(par02),')']});
    ht4=title(axxArray(8),{[t1x2,' - ',t2x2];['Gaussian (\rho = ',num2str(par02e),')']});
    ht5=title(axxArray(11),{[t1x,' - ',t2x];['Gaussian (\rho = ',num2str(par03),')']});
    ht6=title(axxArray(12),{[t1x2,' - ',t2x2];['Gaussian (\rho = ',num2str(par03e),')']});
    
elseif strcmpi(copulaFamily,'Frank') || strcmpi(copulaFamily,'Clayton') ||strcmpi(copulaFamily,'Gumbel')
    ht=title(axxArray(4),{[t1x,' - ',t2x];[char(copulaFamily),' (\theta = ',num2str(par01),')']});
    ht2=title(axxArray(5),{[t1x2,' - ',t2x2];[char(copulaFamily),' (\theta = ',num2str(par02),')']});
end
ht.VerticalAlignment='top';
ht2.VerticalAlignment='top';
h3.VerticalAlignment='top';
ht4.VerticalAlignment='top';
ht5.VerticalAlignment='top';
ht6.VerticalAlignment='top';
% 

% 





 
 ttRho=linspace(inputtimestampsWindowCell{1}(1),inputtimestampsWindowCell{end}(end),length(inputtimestampsWindowCell));
positionLabel2=[];
yLabel2=[];
if strcmpi(copulaFamily,'gaussian')

    corrSpearmanSamplex=gofStatistics.corrSpearmanSamplex;
    corrSpearmanMontex=gofStatistics.corrSpearmanMontex;
    % [hAx,hLine1,hLine2] = plotyy(x,y1,[x',x'],[y2',y3']);
    % plot(axxArray(3),datetime(datevec(ttRho)),cell2mat(couplingParam),'LineWidth',1)
     % 
     couplingParamMat=(cell2mat(cellfun(@(x) x(find(tril(x,-1))),couplingParam,'UniformOutput',0)))';
   
     corrSpearmanSamplexMat=[corrSpearmanSamplex{:}]';
     corrSpearmanMontexMat=[corrSpearmanMontex{:}]';
    [hPlot]=plot(axxArray(9),datetime(datevec(ttRho)),[couplingParamMat]);
     % set(hAx,{'ycolor'},{'k';'k'}) 
     % hLine1(1).LineWidth=1;
     % hLine1(2).LineWidth=1;
     % hLine1(3).LineWidth=1;
     % hLine2(1).LineWidth=1;
     % hLine2(2).LineWidth=1;
    [H,p_value]=tsMann_Kendall(couplingParamMat(:,1),0.05);
 [H2,p_value2]=tsMann_Kendall(couplingParamMat(:,2),0.05);
  [H3,p_value3]=tsMann_Kendall(couplingParamMat(:,3),0.05);

    str=lower(char(copulaFamily));
    idx=regexp([' ' str],'(?<=\s+)\S','start')-1;
    str(idx)=upper(str(idx));
    yll=ylabel(axxArray(9),sprintf('\\rho_{%s}', str));
   
    positionLabel2=[positionLabel2;yll.Position];
    yLabel2=[yLabel2,yll];

    % rho=cellfun(@(x) x(2),couplingParam);
    % plot(axxArray(9),datetime(datevec(ttRho)),(rho),'LineWidth',1)
    % [H,p_value]=tsMann_Kendall((rho),0.05);
    % 
    % str=lower(char(copulaFamily));
    % idx=regexp([' ' str],'(?<=\s+)\S','start')-1;
    % str(idx)=upper(str(idx));
    % yll=ylabel(axxArray(9),sprintf('\\rho_{%s}', str));
    % % ylabel(axxArray(3),'\rho_{Gaussian}')
    % positionLabel2=[positionLabel2;yll.Position];
    % yLabel2=[yLabel2,yll];
elseif strcmpi(copulaFamily,'clayton') || strcmpi(copulaFamily,'gumbel') || strcmpi(copulaFamily,'frank')
    corrSpearmanSamplex=gofStatistics.corrSpearmanSamplex;
    corrSpearmanMontex=gofStatistics.corrSpearmanMontex;
    % [hAx,hLine1,hLine2] = plotyy(x,y1,[x',x'],[y2',y3']);
    % plot(axxArray(3),datetime(datevec(ttRho)),cell2mat(couplingParam),'LineWidth',1)
     [hAx,hLine1,hLine2]=plotyy(axxArray(3),datetime(datevec(ttRho)),cell2mat(couplingParam),[datetime(datevec(ttRho)),datetime(datevec(ttRho))],[[corrSpearmanSamplex{:}]',[corrSpearmanMontex{:}]']);
     set(hAx,{'ycolor'},{'k';'k'}) 
     hLine1.LineWidth=1;
     hLine2(1).LineWidth=1;
     hLine2(2).LineWidth=1;
    [H,p_value]=tsMann_Kendall(cell2mat(couplingParam),0.05);


    str=lower(char(copulaFamily));
    idx=regexp([' ' str],'(?<=\s+)\S','start')-1;
    str(idx)=upper(str(idx));
    yll=ylabel(hAx(1),sprintf('\\theta_{%s}', str));
    yll2=ylabel(hAx(2),'\rho_{Spearman}');
    positionLabel2=[positionLabel2;yll.Position];
    yLabel2=[yLabel2,yll];
end
grid(axxArray(9),'on')
text(axxArray(9),0.35,0.1,['{\it p-value}_{MK,1-2}= ',sprintf('%0.3g',p_value)],'units','normalized','HorizontalAlignment','left')
text(axxArray(9),0.35,0.3,['{\it p-value}_{MK,1-3}= ',sprintf('%0.3g',p_value2)],'units','normalized','HorizontalAlignment','left')
text(axxArray(9),0.35,0.5,['{\it p-value}_{MK,2-3}= ',sprintf('%0.3g',p_value3)],'units','normalized','HorizontalAlignment','left')

xlabel(axxArray(9),'Date (time)')

pval=cellfun(@(x) x{2}.pValueChange,copulaAnalysis.marginalAnalysis);
text(axxArray(2),0.5,0.1,['{\it p-value}_{MK}= ',sprintf('%0.3g',pval(1))],'units','normalized')
text(axxArray(6),0.5,0.1,['{\it p-value}_{MK}= ',sprintf('%0.3g',pval(2))],'units','normalized')
text(axxArray(10),0.5,0.1,['{\it p-value}_{MK}= ',sprintf('%0.3g',pval(3))],'units','normalized')

indexToCopula=find(strcmp(copulaFamily,{gofStatistics.copulaFamily}));
snSample=gofStatistics(indexToCopula).snSample;
corrKendallSampleDelta=gofStatistics(indexToCopula).corrKendallSampleDelta;
corrSpearmanSampleDelta=gofStatistics(indexToCopula).corrSpearmanSampleDelta;

% text(axxArray(5),0.65,0.2,['Sn = ',sprintf('%0.1g',snSample)],'units','normalized','HorizontalAlignment','left')
% text(axxArray(5),0.65,0.3,['\Delta\tau_{Kendall} = ',sprintf('%0.2g',corrKendallSampleDelta)],'units','normalized','HorizontalAlignment','left')
% % text(axxArray(5),0.65,0.4,['\Delta\rho_{Spearman} = ',sprintf('%0.2g',corrSpearmanSampleDelta)],'units','normalized','HorizontalAlignment','left')
% text(axxArray(5),0.65,0.4,['\Delta\rho_{Spearman} = ',sprintf('%0.2g',corrSpearmanSampleDelta)],'units','normalized','HorizontalAlignment','left')

formattedValue = sprintf('%.2g', corrSpearmanSampleDelta); % 
formattedValue2 = sprintf('%.1g', snSample); % 
formattedValue3 = sprintf('%.2g', corrKendallSampleDelta); % 

latexString = sprintf('$\\overline{\\Delta \\rho}_{\\mathrm{Spearman}} = %s$', formattedValue);
% latexString2 = sprintf('$\\overline{S_n} = %s$', formattedValue2);
latexString2 = sprintf('$\\overline{S_n} = %s$', formattedValue2);
latexString3 = sprintf('$\\overline{\\Delta \\tau}_{\\mathrm{Kendall}} = %s$', formattedValue3);
text(axxArray(5),0.65, 0.4, latexString,'units','normalized','HorizontalAlignment','left','Interpreter', 'latex', 'FontSize', 14);
text(axxArray(5),0.65, 0.3, latexString3,'units','normalized','HorizontalAlignment','left','Interpreter', 'latex','FontSize', 14);
text(axxArray(5),0.65, 0.2, latexString2 ,'units','normalized','HorizontalAlignment','left','Interpreter', 'latex','FontSize', 14);

set(axxArray(1),'FontSize',fontSize)
set(axxArray(2),'FontSize',fontSize)
set(axxArray(3),'FontSize',fontSize)
set(axxArray(6),'FontSize',fontSize)
set(axxArray(10),'FontSize',fontSize)




legend1 = legend(axxArray(9),'show');
set(legend1,...
    'String',{'\rho_{Gaussian_{1,2}}','\rho_{Gaussian_{1,3}}','\rho_{Gaussian_{2,3}}'},'EdgeColor','none','Color','none');
grid(axxArray(2),'on')
grid(axxArray(6),'on')

grid(axxArray(10),'on')
axis(axxArray(7),'tight')
axis(axxArray(8),'tight')
axis(axxArray(3),'tight')
axis(axxArray(4),'tight')
axis(axxArray(11),'tight')
axis(axxArray(12),'tight')
 text(axxArray(9),0.05,0.9,labelMark(9),'Units','normalized')

end












