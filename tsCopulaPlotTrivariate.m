function [handles] = tsCopulaTimeVaryingPlot(copulaAnalysis,gofStatistics,varargin)

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
labelMark=(["(a)","(b)","(c)","(d)","(e)","f"]);

%read some parameters

methodology=copulaAnalysis.methodology;

% inputtimestampsWindowCell=copulaAnalysis.copulaParam.inputtimestampsWindowCell;


%set figure and axes properties; keep axes handles as ana array

axxArray=[];
rt=1;
b0=15;
l0=27;
spMan = tsLcSubplotManager(b0, l0, 'CellXSize', round(27*rt), 'CellYSize', round(27*rt), 'gap', [0 0]);
spMan.initFigure;
h=[6,6,6,6,6,6];
b=[7.5,7.5,7.5,7.5,7.5,7.5];
h0=[6.5,14,6.5,14,6.5,14];
b0=[2,2,11.5,11.5,21,21];
for ij=1:length(h0)
    % format    [top left height width]   unlike the usual MATLAB
    % format of [left bottom width height]
    axx=spMan.createAxes(num2str(ij),h0(ij),b0(ij),h(ij),b(ij));
    if ij==1
        % path(path,'/Users/hadi/Downloads/m_map')

        %mercator
        m_proj('albers equal-area','longitudes',[-10 0], ...
            'latitudes',[36 44],'rect','on');  % repeated here so cut-n-paste simplified
        % m_gshhs_f('save','coastLinePortugal','patch',[.7 .9 .7])
        m_usercoast('coastLinePortugal','color','k')
      %         M_GSHHS_F('save',FILENAME) saves the extracted coastline data
%         for the current projection in a file FILENAME. This allows 
%         speedier replotting using M_USERCOAST(FILENAME). 
        m_grid('linestyle','none','linewidth',2,'tickdir','out',...
            'xaxisloc','top','yaxisloc','right','fontsize',6);

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
            if ii~=ir && ii~=ic
                m_text(latlon(1,ii),latlon(2,ii),locString(ii),'color','k','vertical','top','fontweight','bold','horizontalalignment','right');
            end
            m_ruler([.5 .8],.2,3,'fontsize',8);
        end
        axx_=spMan.createAxes(num2str(7),4.7,4.45,3,3.75);

        m_proj('albers equal-area','longitudes',[-9.1 -8.35], ...
            'latitudes',[40.8 41.4],'rect','on');
        % m_gshhs_f('patch',[.7 .9 .7]);
          m_usercoast('coastLinePortugal','color','k')
        m_grid('linestyle','none','linewidth',2,'tickdir','out',...
            'xaxisloc','top','yaxisloc','right','fontsize',6);

        latlon2=[latlon(:,ir),latlon(:,ic)];
        locString2=[locString(:,ir),locString(:,ic)];
        for ii=1:size(latlon2,2)
            m_line(latlon2(1,ii),latlon2(2,ii),'marker','square','markersize',4,'color','r');
            if ii==1
                m_text(latlon2(1,ii),latlon2(2,ii),locString2(ii),'color','k','vertical','top','fontweight','bold');
            else
                m_text(latlon2(1,ii),latlon2(2,ii),locString2(ii),'color','k','vertical','bottom','fontweight','bold');
            end

        end
        xll=axx_.XLim;
        yll=axx_.YLim;
        [longx,latx]=m_xy2ll(xll,yll);
        m_proj('albers equal-area','longitudes',[-10 0], ...
            'latitudes',[36 44],'rect','on');

        [xx,yy]=m_ll2xy(longx,latx);
        set(axx,'NextPlot','add')
        plot(axx,[xx(1) xx(2) xx(2) xx(1) xx(1)],[yy(1) yy(1) yy(2) yy(2) yy(1)],'LineStyle','-','LineWidth',1,'Color','k')


        annotation(gcf,'line',[0.11214953271028 0.162883845126836],...
            [0.814117647058824 0.861176470588235]);
        annotation(gcf,'line',[0.113484646194927 0.164218958611482],...
            [0.788235294117647 0.675294117647059]);

    end
    axxArray=[axxArray,axx];
end


%%% plot time series, joint peaks and threshold curves

% read stationary copula information for plotting time series

yMax=copulaAnalysis.jointExtremes;              %non-stationary peaks
tMax=copulaAnalysis.jointExtremeTimeStamps;     %time of occurrences of non-stationary peaks
thresholdPotNS=copulaAnalysis.thresholdPotNS;   %curves representing non-stationary threshold

% stationary joint peaks were sorted in the sampling process; However,
% non-stationary peaks aren't necessarily sorted; since we use
% non-stationary peaks as a way of color-coding the data, we'll have them
% sorted at this stage;

[~,iyMax]=sort(mean(yMax,2),'descend');
yMax=yMax(iyMax,:);
tMax=tMax(iyMax,:);
%
scatter3(axxArray(3),yMax(:,1),yMax(:,2),yMax(:,3),[],yMax(:,1),'filled')
hold(axxArray(3),'on')
resampleLevel=copulaAnalysis.resampleLevel;
scatterMontCarl03d=scatter3(axxArray(3),resampleLevel{1}(:,1),resampleLevel{1}(:,2),resampleLevel{1}(:,3));
set(scatterMontCarl03d,'LineWidth',1,'Marker','o','MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor',[0.65,0.65,0.65],...
    'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)

ylblx2=sprintf('{\\it%s}_{%s}',locString(1),ylbl{1});
xlabel(axxArray(3),ylblx2);
ylblx2=sprintf('{\\it%s}_{%s}',locString(2),ylbl{2});
ylabel(axxArray(3),ylblx2);
ylblx3=sprintf('{\\it%s}_{%s}',locString(3),ylbl{3});
zlabel(axxArray(3),ylblx3);
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


couplingParam=copulaAnalysis.copulaParam.rho;
copulaFamily=copulaAnalysis.copulaParam.family;
%perform plotting
scatterMontCarl01=scatter(axxArray(2),resampleLevel{1}(:,1), resampleLevel{1}(:,2));
set(scatterMontCarl01,'LineWidth',1,'Marker','o','MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor',[0.65,0.65,0.65],...
    'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.6)
hold(axxArray(2),'on')
scatter(axxArray(2),yMaxLevel(:,1), yMaxLevel(:,2),[],yMaxLevel(:,1),'filled');

scatterMontCarl02=scatter(axxArray(4),resampleLevel{1}(:,1), resampleLevel{1}(:,3));
set(scatterMontCarl02,'LineWidth',1,'Marker','o','MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor',[0.65,0.65,0.65],...
    'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.6)
hold(axxArray(4),'on')
scatter(axxArray(4),yMaxLevel(:,1), yMaxLevel(:,3),[],yMaxLevel(:,1),'filled');

scatterMontCarl03=scatter(axxArray(6),resampleLevel{1}(:,2), resampleLevel{1}(:,3));
set(scatterMontCarl03,'LineWidth',1,'Marker','o','MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor',[0.65,0.65,0.65],...
    'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.6)
hold(axxArray(6),'on')
scatter(axxArray(6),yMaxLevel(:,2), yMaxLevel(:,3),[],yMaxLevel(:,1),'filled');

if strcmpi(copulaFamily,'Gaussian')
    par01=round(couplingParam{1}(1,2)*100)/100;
    par02=round(couplingParam{1}(1,3)*100)/100;
    par03=round(couplingParam{1}(2,3)*100)/100;
elseif strcmpi(copulaFamily,'Frank') || strcmpi(copulaFamily,'Clayton') ||strcmpi(copulaFamily,'Gumbel')
    par01=round(couplingParam{1}(1)*100)/100;

end

if strcmpi(copulaFamily,'Gaussian')
    ht=title(axxArray(2),{['Gaussian (\rho = ',num2str(par01),')']});
    ht2=title(axxArray(4),{['Gaussian (\rho = ',num2str(par02),')']});
    ht3=title(axxArray(6),{['Gaussian (\rho = ',num2str(par03),')']});
elseif strcmpi(copulaFamily,'Frank') || strcmpi(copulaFamily,'Clayton') ||strcmpi(copulaFamily,'Gumbel')
    ht=title(axxArray(2),{[char(copulaFamily),' (\theta = ',num2str(par01),')']});
    ht2=title(axxArray(4),{[char(copulaFamily),' (\theta = ',num2str(par01),')']});
    ht3=title(axxArray(6),{[char(copulaFamily),' (\theta = ',num2str(par01),')']});
end

ht.VerticalAlignment='top';
ht2.VerticalAlignment='top';
ht3.VerticalAlignment='top';

% labelMarks2=["Loc1","Loc2","Loc3","Loc4","Loc5","Loc6","Loc7"];
% ylblx=cellfun(@(x,y) strtrim(x(1:y-1)),ylbl,regexp(ylbl,'('),'UniformOutput',0)
ylblx=sprintf('{\\it%s}_{%s}',locString(1),ylbl{1});
xlabel(axxArray(2),ylblx);
ylblx=sprintf('{\\it%s}_{%s}',locString(2),ylbl{2});
ylabel(axxArray(2),ylblx);

ylblx=sprintf('{\\it%s}_{%s}',locString(1),ylbl{1});
xlabel(axxArray(4),ylblx);
ylblx=sprintf('{\\it%s}_{%s}',locString(3),ylbl{3});
ylabel(axxArray(4),ylblx);

ylblx=sprintf('{\\it%s}_{%s}',locString(2),ylbl{2});
xlabel(axxArray(6),ylblx);
ylblx=sprintf('{\\it%s}_{%s}',locString(3),ylbl{3});
ylabel(axxArray(6),ylblx);


grid(axxArray(2),'on')
grid(axxArray(4),'on')
grid(axxArray(6),'on')
text(axxArray(2),0.05,0.9,labelMark(2),'Units','normalized')
text(axxArray(4),0.05,0.9,labelMark(4),'Units','normalized')
text(axxArray(6),0.05,0.9,labelMark(6),'Units','normalized')
%annotate correlations
corrSp01M=corr(resampleLevel{1});
corrSp01M=round(corrSp01M*100)/100;
corrSp01S=corr(yMaxLevel);
corrSp01S=round(corrSp01S*100)/100;


text(axxArray(2),0.5,0.12,sprintf('\\rho_{Spearman, S} =  %2g\n\\rho_{Spearman, MC} =  %2g', corrSp01S(1,2),corrSp01M(1,2)),'Units','normalized')
text(axxArray(4),0.5,0.12,sprintf('\\rho_{Spearman, S} =  %2g\n\\rho_{Spearman, MC} =  %2g', corrSp01S(1,3),corrSp01M(1,3)),'Units','normalized')
text(axxArray(6),0.5,0.12,sprintf('\\rho_{Spearman, S} =  %2g\n\\rho_{Spearman, MC} =  %2g', corrSp01S(2,3),corrSp01M(2,3)),'Units','normalized')


indexToCopula=find(strcmp(copulaFamily,{gofStatistics.copulaFamily}));
snSample=gofStatistics(indexToCopula).snSample;
corrKendallSampleDelta=gofStatistics(indexToCopula).corrKendallSampleDelta;
corrSpearmanSampleDelta=gofStatistics(indexToCopula).corrSpearmanSampleDelta;

text(axxArray(5),0.1,0.2,['Sn = ',sprintf('%0.3g',snSample)],'units','normalized','HorizontalAlignment','left','FontSize',20)
text(axxArray(5),0.1,0.5,['\Delta\tau_{Kendall} = ',sprintf('%0.3g',corrKendallSampleDelta)],'units','normalized','HorizontalAlignment','left','FontSize',20)
text(axxArray(5),0.1,0.8,['\Delta\tau_{Spearman} = ',sprintf('%0.3g',corrSpearmanSampleDelta)],'units','normalized','HorizontalAlignment','left','FontSize',20)
set(axxArray(5),'xtick',[])
set(axxArray(5),'xticklabel',[])
set(axxArray(5),'ytick',[])
set(axxArray(5),'yticklabel',[])

set(axxArray(2),'FontSize',fontSize)
set(axxArray(3),'FontSize',fontSize)
set(axxArray(4),'FontSize',fontSize)
set(axxArray(5),'FontSize',fontSize)
set(axxArray(6),'FontSize',fontSize)

text(axxArray(5),0.05,0.9,labelMark(5),'Units','normalized')
text(axxArray(3),0.05,0.9,labelMark(3),'Units','normalized')

handles=[];
end






