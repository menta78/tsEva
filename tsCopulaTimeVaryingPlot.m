function [handles] = tsCopulaTimeVaryingPlot(resampleLevel, copulaAnalysis, varargin)

args.xlbl = 'X';
args.ylbl = 'Y';
args.figPosition = [100, 100, 800, 800];
args.fontSize = 15;

args = tsEasyParseNamedArgs(varargin, args);
xlbl = args.xlbl;
ylbl = args.ylbl;
figPosition = args.figPosition;
fontSize = args.fontSize;

nSubPlot=size(resampleLevel,2);
nFig=fix(nSubPlot/16)+1;

yMaxLevel=copulaAnalysis.jointExtremes;
Rho=copulaAnalysis.copulaParam.rho;
inputtimestampsWindowCell=copulaAnalysis.copulaParam.inputtimestampsWindowCell;
for ij=1:nFig

rt=1;
b0=27;
l0=21;
spMan = tsLcSubplotManager(b0, l0, 'CellXSize', round(21*rt), 'CellYSize', round(27*rt), 'gap', [0.03 0.01]);


spMan.initFigure;

numberofverticalpanels=4;
numberofhorizontalpanels=4;
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


axn=1;
for ik=1:numberofhorizontalpanels
for ij=1:numberofverticalpanels

eval(['b',num2str(ij),'=','b0-(numberofverticalpanels-ij)*(nb)-(numberofverticalpanels-ij+1)*(uddelta)'])
eval(['l',num2str(ik),'=','l0-(numberofhorizontalpanels-ik+1)*(nl)-(numberofhorizontalpanels-ik)*(lrdelta)+lmargin'])

axx=spMan.createAxes(num2str(axn),eval(['b',num2str(ij)]),eval(['l',num2str(ik)]),nb,nl);

if axn>nSubPlot
break
end
axn=axn+1;


rsmplSctr = scatter(axx,resampleLevel{1}(:,1), resampleLevel{1}(:,2));
resampleLevel=resampleLevel(2:end);
% rsmplSctrObj = findall(rsmplSctr(1), 'type', 'line');
% rsmplDist1Ax = rsmplSctr(2);
% rsmplDist2Ax = rsmplSctr(3);
hold on;
ymaxSctr = scatter(yMaxLevel{1}(:,1), yMaxLevel{1}(:,2), 'markerfacecolor', 'r');
yMaxLevel=yMaxLevel(2:end);
xlabel(xlbl);
ylabel(ylbl);
grid on;
ax = gca;
ax.FontSize = fontSize;
nr=round(Rho{1}(2)*1000)/1000;
t1x=datestr(inputtimestampsWindowCell{1}(1),'yyyy');
t2x=datestr(inputtimestampsWindowCell{1}(end),'yyyy');
Rho=Rho(2:end);
inputtimestampsWindowCell=inputtimestampsWindowCell(2:end);
ht=title({[t1x,' - ',t2x];['Rho = ',num2str(nr)]});
ht.VerticalAlignment='top';

% lgnd = legend([ymaxSctr, rsmplSctr], {'Peaks', ['Joint Distribution' newline 'Montecarlo']}, 'fontsize', fontSize);
% lgndPos = lgnd.Position;
% newHigh = lgndPos(4)*2;
% lgnd.Position = [rsmplDist2Ax.Position(1), rsmplDist1Ax.Position(2) + rsmplDist1Ax.Position(4) - newHigh, lgndPos(3), newHigh];
% set(fig, 'paperpositionmode', 'auto');
end
end




end



handles=[];
% handles = [fig, rsmplSctrObj, ymaxSctr, lgnd, rsmplSctr(1), rsmplSctr(2), rsmplSctr(3)];
