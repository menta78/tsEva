function [handles] = tsCopulaTimeVaryingPlot(resampleLevel, copulaAnalysis, varargin)

args.xlbl = 'X';
args.ylbl = 'Y';
args.figPosition = [100, 100, 800, 800];
args.fontSize = 15;
args.numberofverticalpanels=4;
args.numberofhorizontalpanels=4;
args = tsEasyParseNamedArgs(varargin, args);
xlbl = args.xlbl;
ylbl = args.ylbl;
figPosition = args.figPosition;
fontSize = args.fontSize;
numberofverticalpanels=args.numberofverticalpanels;
numberofhorizontalpanels=args.numberofhorizontalpanels;
if iscell(resampleLevel)
    nSubPlot=size(resampleLevel,2);
else
    nSubPlot=1;
end
if mod(nSubPlot,numberofverticalpanels*numberofhorizontalpanels)==0
    nFig=nSubPlot/(numberofverticalpanels*numberofhorizontalpanels);
else
    nFig=fix(nSubPlot/(numberofverticalpanels*numberofhorizontalpanels))+1;
end
yMaxLevel=copulaAnalysis.jointExtremes;
Rho=copulaAnalysis.copulaParam.rho;
if iscell(resampleLevel)
    inputtimestampsWindowCell=copulaAnalysis.copulaParam.inputtimestampsWindowCell;
else
    inputtimestampsWindowCell=copulaAnalysis.marginalAnalysis{1}{2}.timeStamps;
end
axn=1;

for ij=1:nFig

    rt=1;
    b0=27;
    l0=21;
    spMan = tsLcSubplotManager(b0, l0, 'CellXSize', round(21*rt), 'CellYSize', round(27*rt), 'gap', [0.03 0.01]);

    spMan.initFigure;

    ratiovertical=8;
    ratiohoriz=8;
    bmargin=0.2;
    leftpanels=[-numberofverticalpanels -numberofverticalpanels+1;1 -ratiovertical]\[-b0+bmargin;0];
    nb=leftpanels(1);
    uddelta=leftpanels(2);
    lmargin=1;

    horizpanels=[-numberofhorizontalpanels -numberofhorizontalpanels+1;1 -ratiohoriz]\[-l0+lmargin;0];
    nl=horizpanels(1);
    lrdelta=horizpanels(2);

    for ik=1:numberofhorizontalpanels
        for ij=1:numberofverticalpanels

            if axn>nSubPlot
                break
            end
            eval(['b',num2str(ij),'=','b0-(numberofverticalpanels-ij)*(nb)-(numberofverticalpanels-ij+1)*(uddelta);'])
            eval(['l',num2str(ik),'=','l0-(numberofhorizontalpanels-ik+1)*(nl)-(numberofhorizontalpanels-ik)*(lrdelta)+lmargin;'])

            axx=spMan.createAxes(num2str(axn),eval(['b',num2str(ij)]),eval(['l',num2str(ik)]),nb,nl);

            axn=axn+1;
            if iscell(resampleLevel)
                rsmplSctr = scatter(axx,resampleLevel{1}(:,1), resampleLevel{1}(:,2));
                resampleLevel=resampleLevel(2:end);
                hold on;
                ymaxSctr = scatter(yMaxLevel{1}(:,1), yMaxLevel{1}(:,2), 'markerfacecolor', 'r');
                yMaxLevel=yMaxLevel(2:end);
                if strcmp(copulaAnalysis.copulaParam.family,'Gaussian') ||strcmp(copulaAnalysis.copulaParam.family,'t')
                    nr=round(Rho{1}(2)*1000)/1000;
                elseif strcmp(copulaAnalysis.copulaParam.family,'Frank') || strcmp(copulaAnalysis.copulaParam.family,'Clayton') ||strcmp(copulaAnalysis.copulaParam.family,'Gumbel')
                    nr=round(Rho{1}(1));
                end
                t1x=datestr(inputtimestampsWindowCell{1}(1),'yyyy');
                t2x=datestr(inputtimestampsWindowCell{1}(end),'yyyy');
                Rho=Rho(2:end);
                inputtimestampsWindowCell=inputtimestampsWindowCell(2:end);
            else
                rsmplSctr = scatter(axx,resampleLevel(:,1), resampleLevel(:,2));
                hold on
                ymaxSctr = scatter(yMaxLevel(:,1), yMaxLevel(:,2), 'markerfacecolor', 'r');
                if strcmp(copulaAnalysis.copulaParam.family,'Gaussian') ||strcmp(copulaAnalysis.copulaParam.family,'t')
                    nr=round(Rho(2)*1000)/1000;
                elseif strcmp(copulaAnalysis.copulaParam.family,'Frank') || strcmp(copulaAnalysis.copulaParam.family,'Clayton') ||strcmp(copulaAnalysis.copulaParam.family,'Gumbel')
                    nr=round(Rho(1));
                end
                t1x=datestr(inputtimestampsWindowCell(1),'yyyy');
                t2x=datestr(inputtimestampsWindowCell(end),'yyyy');

            end

            xlabel(xlbl);
            ylabel(ylbl);
            grid on;
            ax = gca;
            ax.FontSize = fontSize;

            ht=title({[t1x,' - ',t2x];['Rho = ',num2str(nr)]});
            ht.VerticalAlignment='top';

        end
    end

end

handles=[];

