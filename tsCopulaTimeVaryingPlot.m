function [handles] = tsCopulaTimeVaryingPlot(copulaAnalysis, gofStatistics,varargin)

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

args.xlbl = 'X';
args.ylbl = 'Y';
args.fontSize = 15;
args.numberofverticalpanels=2;
args.numberofhorizontalpanels=3;

args = tsEasyParseNamedArgs(varargin, args);

xlbl = args.xlbl;
ylbl = args.ylbl;
fontSize = args.fontSize;
numberofverticalpanels=args.numberofverticalpanels;
numberofhorizontalpanels=args.numberofhorizontalpanels;

%

% obtain original input series data based on data contained withing copulaAnalysis
marginalAnalysis=copulaAnalysis.marginalAnalysis;
nonStatSeriesCell=cellfun(@(x) x{2}.nonStatSeries,marginalAnalysis,'UniformOutput',0);
timeStampsCell=cellfun(@(x) x{2}.timeStamps,marginalAnalysis,'UniformOutput',0);
nonStatSeries=cell2mat(nonStatSeriesCell);
timeStamps=cell2mat(timeStampsCell);



if iscell(copulaAnalysis.jointExtremes) %time-varying copula
    %read information about joint peaks and thresholds
    yMax=copulaAnalysis.jointExtremes;
    tMax=copulaAnalysis.jointExtremeTimeStamps;
    thresholdPotNS=copulaAnalysis.thresholdPotNS;

    % obtain information about the peaks lumping together all peaks in each
    % time window when the calculated copula was of a time-varying type
    ytMax=unique([cell2mat(yMax'),cell2mat(tMax')],'rows','stable');
    numVar=size(ytMax,2)/2;
    yMax=ytMax(:,1:numVar);
    tMax=ytMax(:,numVar+1:end);

    %first plot is a simple plot of input time series and sampled joint
    %peaks
    for jx=1:size(nonStatSeries,2)
        subplot(size(nonStatSeries,2),1,jx)
        plot(datetime(datevec(timeStamps(:,jx))),nonStatSeries(:,jx))
        hold on
        plot(datetime(datevec(timeStamps(:,jx))),thresholdPotNS(:,jx))
        plot(datetime(datevec(tMax(:,jx))),yMax(:,jx),'.r')
    end
    hFig1=get(gcf,'Children');
    %determine nSubPlot and nFig based on number of time windows withing which
    %a copula has been estimated
    nSubPlot=size(copulaAnalysis.jointExtremes,2);
    if mod(nSubPlot,numberofverticalpanels*numberofhorizontalpanels)==0
        nFig=nSubPlot/(numberofverticalpanels*numberofhorizontalpanels);
    else
        nFig=fix(nSubPlot/(numberofverticalpanels*numberofhorizontalpanels))+1;
    end

    % next plot concerns a scatter plot of sampled joint peaks and resampled
    % Monte-Carlo values that were sampled in accordance with the chosen copula
    % model
    yMaxLevel=copulaAnalysis.jointExtremes;
    Rho=copulaAnalysis.copulaParam.rho;
    resampleLevel=copulaAnalysis.resampleLevel;
    inputtimestampsWindowCell=copulaAnalysis.copulaParam.inputtimestampsWindowCell;

    %extract gof information
    aicSample=gofStatistics.aicSample;
    bicSample=gofStatistics.bicSample;
    corrKendallSampleDelta=gofStatistics.corrKendallSampleDelta;
    corrPearsonSampleDelta=gofStatistics.corrPearsonSampleDelta;
    corrSpearmanSampleDelta=gofStatistics.corrSpearmanSampleDelta;
    llSample=gofStatistics.llSample;
    sncSample=gofStatistics.sncSample;

    %calculate limits of x and y axes
    minXYAxes=zeros(1,numVar);
    maxXYAxes=zeros(1,numVar);
    for iVarC=1:numVar

        [minXYAxes(iVarC)]=min(cellfun(@(x,y) min([min(x(:,iVarC)),min(y(:,iVarC))]),yMaxLevel,resampleLevel,'UniformOutput',1));
        [maxXYAxes(iVarC)]=max(cellfun(@(x,y) max([max(x(:,iVarC)),max(y(:,iVarC))]),yMaxLevel,resampleLevel,'UniformOutput',1));
    end
    minmaxXYAxes=[minXYAxes;maxXYAxes];
    minmaxXX=minmaxXYAxes(:,1);
    minmaxYY=minmaxXYAxes(:,2);
    %lopp through number of figures
    axn=1;
    for ij=1:nFig

        %set some parameters that determine sizing of the figure
        rt=1;
        b0=27;
        l0=27;
        spMan = tsLcSubplotManager(b0, l0, 'CellXSize', round(27*rt), 'CellYSize', round(27*rt), 'gap', [0 0]);
        spMan.initFigure;

        %determine dimensions of the axes and other parameters related with
        %the axes
        ratiovertical=6;
        ratiohoriz=3;
        bmargin=1.2;
        bmargin2=0.1;
        leftpanels=[-numberofverticalpanels -numberofverticalpanels+1;1 -ratiovertical]\[-0.97*(b0)+bmargin+bmargin2-1;0];
        nb=leftpanels(1);
        uddelta=leftpanels(2);
        lmargin=2.5;
        lmargin2=0.1;
        horizpanels=[-numberofhorizontalpanels -numberofhorizontalpanels+1;1 -ratiohoriz]\[-0.97*(l0)+lmargin+lmargin2-1;0];
        nl=horizpanels(1);
        lrdelta=horizpanels(2);

        for ik=1:numberofhorizontalpanels
            for ij=1:numberofverticalpanels

                if axn>nSubPlot
                    break
                end
                %create axis, top left corner is (0,0), positive x is
                %towards right and positive y is downward. The usual MATLAB
                %convention of defining origin of the axes and defining its
                %dimensions (width and height) is used
                eval(['b',num2str(ij),'=','(ij)*(nb)+(ij-1)*(uddelta)+bmargin2;'])
                eval(['l',num2str(ik),'=','(ik-1)*(nl)+(ik-1)*(lrdelta)+lmargin;'])
                axx=spMan.createAxes(num2str(axn),eval(['b',num2str(ij)]),eval(['l',num2str(ik)]),nb,nl);
                axn=axn+1;

                %plot a scatter plot
                rsmplSctr = scatter(axx,resampleLevel{1}(:,1), resampleLevel{1}(:,2));
                resampleLevel=resampleLevel(2:end);
                hold on;
                ymaxSctr = scatter(yMaxLevel{1}(:,1), yMaxLevel{1}(:,2), 'markerfacecolor', 'r');
                yMaxLevel=yMaxLevel(2:end);

             
                %annotate correlation parameter of the copula and date of each time window
                if strcmpi(copulaAnalysis.copulaParam.family,'Gaussian') ||strcmpi(copulaAnalysis.copulaParam.family,'t')
                    nr=round(Rho{1}(2)*1000)/1000;
                elseif strcmpi(copulaAnalysis.copulaParam.family,'Frank') || strcmpi(copulaAnalysis.copulaParam.family,'Clayton') ||strcmpi(copulaAnalysis.copulaParam.family,'Gumbel')
                    nr=round(Rho{1}(1)*100)/100;
                end
                t1x=datestr(inputtimestampsWindowCell{1}(1),'yyyy');
                t2x=datestr(inputtimestampsWindowCell{1}(end),'yyyy');
                Rho=Rho(2:end);
                inputtimestampsWindowCell=inputtimestampsWindowCell(2:end);

                if strcmpi(copulaAnalysis.copulaParam.family,'Gaussian') ||strcmpi(copulaAnalysis.copulaParam.family,'t')
                    ht=title({[t1x,' - ',t2x];['Rho = ',num2str(nr)]});
                elseif strcmpi(copulaAnalysis.copulaParam.family,'Frank') || strcmpi(copulaAnalysis.copulaParam.family,'Clayton') ||strcmpi(copulaAnalysis.copulaParam.family,'Gumbel')
                    ht=title({[t1x,' - ',t2x];['Alpha = ',num2str(nr)]});
                end

                ht.VerticalAlignment='top';

                xlabel(xlbl);
                ylabel(ylbl);
                grid on;
                ax = gca;
                ax.FontSize = fontSize;

                xlim(minmaxXX)
                ylim(minmaxYY)
                text(0.65,0.1,['AIC = ',num2str(round(aicSample(1)*10)/10)],'units','normalized')
                aicSample=aicSample(2:end);
                text(0.65,0.06,['BIC = ',num2str(round(bicSample(1)*10)/10)],'units','normalized')
                bicSample=bicSample(2:end);
                text(0.65,0.02,['ll = ',num2str(round(llSample(1)*10)/10)],'units','normalized')
                llSample=llSample(2:end);
                text(0.65,0.14,['SnC = ',num2str(round(sncSample(1)*100)/100)],'units','normalized')
                sncSample=sncSample(2:end);
                text(0.40,0.26,[' | \tau_{Mn} - \tau_{Sm} |_{K}= ',sprintf('%.1E',corrKendallSampleDelta(1))],'units','normalized','interpreter','tex')
                corrKendallSampleDelta=corrKendallSampleDelta(2:end);
                text(0.40,0.22,[' | \tau_{Mn} - \tau_{Sm} |_{P}= ',sprintf('%.1E',corrPearsonSampleDelta(1))],'units','normalized','interpreter','tex')
                corrPearsonSampleDelta=corrPearsonSampleDelta(2:end);
                text(0.40,0.18,[' | \tau_{Mn} - \tau_{Sm} |_{S}= ',sprintf('%.1E',corrSpearmanSampleDelta(1))],'units','normalized','interpreter','tex')
                corrSpearmanSampleDelta=corrSpearmanSampleDelta(2:end);
 
 
            end
        end

    end
    hFig2=get(gcf,'Children');
    handles=[hFig1;hFig2];
else %stationary copula
    %read information about the peaks and thresholds
    yMax=copulaAnalysis.jointExtremes; %non-stationary peaks
    tMax=copulaAnalysis.jointExtremeTimeStamps;
    thresholdPotNS=copulaAnalysis.thresholdPotNS; %non-stationary threshold parameters
    %generate first plot: plotting of input time series and thresholds
    %along with sampled joint peaks
    for jx=1:size(nonStatSeries,2)

        subplot(size(nonStatSeries,2),1,jx)
        plot(datetime(datevec(timeStamps(:,jx))),nonStatSeries(:,jx))
        hold on

        plot(datetime(datevec(timeStamps(:,jx))),thresholdPotNS(:,jx))
        plot(datetime(datevec(tMax(:,jx))),yMax(:,jx),'.r')
    end

    hFig1=get(gcf,'Children');

    % read joint peaks and other related information
    yMaxLevel=copulaAnalysis.jointExtremes;
    Rho=copulaAnalysis.copulaParam.rho;

    aicSample=gofStatistics.aicSample;
    bicSample=gofStatistics.bicSample;
    corrKendallSampleDelta=gofStatistics.corrKendallSampleDelta;
    corrPearsonSampleDelta=gofStatistics.corrPearsonSampleDelta;
    corrSpearmanSampleDelta=gofStatistics.corrSpearmanSampleDelta;
    llSample=gofStatistics.llSample;
    sncSample=gofStatistics.sncSample;
    % set-up figure dimensions
    b0=20;
    l0=20;
    spMan = tsLcSubplotManager(b0, l0, 'CellXSize', l0, 'CellYSize', b0, 'gap', [0.00 0.00]);
    spMan.initFigure;

    %set-up axes dimensions
    bmargin=2;
    nb=b0-bmargin;

    lmargin=3;
    nl=l0-lmargin;

    ij=1;
    eval(['b',num2str(1),'=','b0-bmargin;'])
    eval(['l',num2str(1),'=','lmargin;'])
    axx=spMan.createAxes(num2str(1),eval(['b',num2str(ij)]),eval(['l',num2str(ij)]),nb,nl);%measured from top left corner positive downward and to the right

    %plot scatter plots of joint peaks and resampled values that come from the
    %Monte-Carlo simulation
    scatter(axx,copulaAnalysis.resampleLevel(:,1), copulaAnalysis.resampleLevel(:,2));
    hold on
    scatter(yMaxLevel(:,1), yMaxLevel(:,2), 'markerfacecolor', 'r');
    %annotate correlation parameter of the copula and the time duration of
    %input series
    if strcmpi(copulaAnalysis.copulaParam.family,'Gaussian') ||strcmpi(copulaAnalysis.copulaParam.family,'t')
        nr=round(Rho(2)*1000)/1000;
    elseif strcmpi(copulaAnalysis.copulaParam.family,'Frank') || strcmpi(copulaAnalysis.copulaParam.family,'Clayton') ||strcmpi(copulaAnalysis.copulaParam.family,'Gumbel')
        nr=round(Rho(1)*100)/100;
    end
    t1x=datestr(timeStamps(1,1),'yyyy');
    t2x=datestr(timeStamps(end,1),'yyyy');

    if strcmpi(copulaAnalysis.copulaParam.family,'Gaussian') ||strcmpi(copulaAnalysis.copulaParam.family,'t')
        ht=title({[t1x,' - ',t2x];['Rho = ',num2str(nr)]});
    elseif strcmpi(copulaAnalysis.copulaParam.family,'Frank') || strcmpi(copulaAnalysis.copulaParam.family,'Clayton') ||strcmpi(copulaAnalysis.copulaParam.family,'Gumbel')
        ht=title({[t1x,' - ',t2x];['Alpha = ',num2str(nr)]});
    end

    ht.VerticalAlignment='top';

    xlabel(xlbl);
    ylabel(ylbl);
    grid on;
    ax = gca;
    ax.FontSize = fontSize;


    text(0.65,0.1,['AIC = ',num2str(round(aicSample(1)*10)/10)],'units','normalized')
    aicSample=aicSample(2:end);
    text(0.65,0.06,['BIC = ',num2str(round(bicSample(1)*10)/10)],'units','normalized')
    bicSample=bicSample(2:end);
    text(0.65,0.02,['ll = ',num2str(round(llSample(1)*10)/10)],'units','normalized')
    llSample=llSample(2:end);
    text(0.65,0.14,['SnC = ',num2str(round(sncSample(1)*100)/100)],'units','normalized')
    sncSample=sncSample(2:end);
    text(0.65,0.26,[' | \tau_{Mn} - \tau_{Sm} |_{K}= ',sprintf('%.1E',corrKendallSampleDelta(1))],'units','normalized','interpreter','tex')
    corrKendallSampleDelta=corrKendallSampleDelta(2:end);
    text(0.65,0.22,[' | \tau_{Mn} - \tau_{Sm} |_{P}= ',sprintf('%.1E',corrPearsonSampleDelta(1))],'units','normalized','interpreter','tex')
    corrPearsonSampleDelta=corrPearsonSampleDelta(2:end);
    text(0.65,0.18,[' | \tau_{Mn} - \tau_{Sm} |_{S}= ',sprintf('%.1E',corrSpearmanSampleDelta(1))],'units','normalized','interpreter','tex')
    corrSpearmanSampleDelta=corrSpearmanSampleDelta(2:end);

    hFig2=get(gcf,'Children');

    handles=[hFig1;hFig2];

end


end




 

