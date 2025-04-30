function  [copulaAnalysis,hFig1]=tsCopulaJointReturnPeriod(copulaAnalysis,axxArray,varargin)
%tsCopulaPlotJointReturnPeriod plotting of multivariate return periods

% tsCopulaPlotJointReturnPeriod(copulaAnalysis,varargin)
%       % plots the multivariate return period according to AND/OR Scenarios
%https://doi.org/10.1002/2015WR017225
% Parts of the code were reworked from MvCAT toolbox https://doi.org/10.1002/2016WR020242
% input:
%  copulaAnalysis          - a variable of type structure provided as the
%                            output of tsCopulaCompoundGPD or
%                            tsCopulaCompoundGPDMontecarlo functions



% output:
%  :                -
%

% M.H.Bahmanpour, 2024
colorPool='rgbcmyk';
args.xlbl = 'X';
args.ylbl = 'Y';
args.RL=[10,20,50];
args.copulaAnalysisForCount=[];
args.timeWindowNonStat=[];
args = tsEasyParseNamedArgs(varargin, args);

copulaAnalysisForCount=args.copulaAnalysisForCount;
timeWindowNonStat=args.timeWindowNonStat;
xlbl = args.xlbl;
ylbl = args.ylbl;
RL=args.RL;
cRL=colorPool(1:length(RL));
margDist = copulaAnalysis.methodology;

%append RL to copulaAnalysis for the output
copulaAnalysis.RL=RL;

% read some information from the copulaAnalysis structure file
jointExtremes=copulaAnalysis.jointExtremes;
Family=copulaAnalysis.copulaParam.family;
if isfield(copulaAnalysis.copulaParam,"rhoMean")
PAR=copulaAnalysis.copulaParam.rhoMean;
else
PAR=copulaAnalysis.copulaParam.rho;
end
% PAR2=copulaAnalysis.copulaParam.nu;

timeIndexArray=copulaAnalysis.timeIndexArray;
timeIndexArray=num2cell(timeIndexArray);
if strcmpi(margDist,'gpd')
    margDist='gp';
    eps =cellfun(@(x)  x{1}(2).parameters.epsilon,copulaAnalysis.marginalAnalysis,'UniformOutput',1);
    sig=cellfun(@(x)  x{1}(2).parameters.sigma,copulaAnalysis.marginalAnalysis,'UniformOutput',0);
    thr=cellfun(@(x)  x{1}(2).parameters.threshold,copulaAnalysis.marginalAnalysis,'UniformOutput',0);

    % assuming that the scaling of the univariates (number of peaks/number
    % of years) remain almost constant with time
    nYear =cellfun(@(x)  (x{1}(2).parameters.timeHorizonEnd-x{1}(2).parameters.timeHorizonStart)/365,copulaAnalysis.marginalAnalysis,'UniformOutput',1);
    nPeak=cellfun(@(x)  length(x{1}(2).objs.peakIndexes),copulaAnalysis.marginalAnalysis,'UniformOutput',1);
    %scaling to account for average number of peaks per year- needed only when GPD return levels are to be calculated
    Scl=nPeak./nYear;

elseif strcmpi(margDist,'gev')
    eps=cellfun(@(x)  x{1}(1).parameters.epsilon,copulaAnalysis.marginalAnalysis,'UniformOutput',1);
    sig=cellfun(@(x)  x{1}(1).parameters.sigma,copulaAnalysis.marginalAnalysis,'UniformOutput',0);
    thr=cellfun(@(x)  x{1}(1).parameters.mu,copulaAnalysis.marginalAnalysis,'UniformOutput',0);
    Scl=1;
end

%define a uniform square (or hypercube) over which analytical copulaCDF would be calculated
% Since we are only concerned with bivariate return period, this is always a 2-by-2 matrix

x = [linspace(1e-5,.95,400) linspace(0.95+1e-5,1-1e-5,500)];

[xx,yy] = meshgrid(x,x);
U = [xx(:),yy(:)];


% extend the uniform square to the length of windows
uCell=repmat({U},1,size(PAR,2));
familyCell=repmat({Family},1,size(PAR,2));

% calculate copula cdf across each window

yCell=cellfun(@(x,y,z) copulacdf(z,x,y),uCell,PAR,familyCell,'UniformOutput',0);

% % extract window information
%
% timeStampsByTimeWindow=copulaAnalysis.copulaParam.timeStampsByTimeWindow;


% calculate the OR return period
if strcmpi(margDist,'gev') 
rPCellOr=cellfun(@(x) 1./(1-x),yCell,'UniformOutput',0);
else
    numJointOrPeaks=cellfun(@(x) length(x),copulaAnalysisForCount.jointExtremes,'UniformOutput',0);
numYearWindow=repmat({timeWindowNonStat},1,size(numJointOrPeaks,2));
scaling=cellfun(@(x,y) x/y,numYearWindow,numJointOrPeaks,'UniformOutput',0);
rPCellOr=cellfun(@(x,y) y./(1-x),yCell,scaling,'UniformOutput',0);

end

% Sort joint return periods
[RpSortCell, IdRpSrtCell]=cellfun(@(x) sort(x),rPCellOr,'UniformOutput',0);

% Sort associated uniform marginal probabilities accordingly
U1cell=cellfun(@(x,y) x(y,1),uCell,IdRpSrtCell,'UniformOutput',0);
U2cell=cellfun(@(x,y) x(y,2),uCell,IdRpSrtCell,'UniformOutput',0);

% define return level contours and their acceptable lower and upper bounds
RpLb = RL - 0.005*RL;
RpUb = RL + 0.005*RL;

% convert marginal distribution parameters to cell format
epsCell=repmat({eps},1,size(PAR,2));
sigCell=repmat({sig},1,size(PAR,2));
thrCell=repmat({thr},1,size(PAR,2));
SclCell=repmat({Scl},1,size(PAR,2));
margDistCell=repmat({margDist},1,size(PAR,2));

%to store output information
cellRP=cell(size(RL));

%loop through each return period contour

% Find indices associated with each probability contour


for ij=1:size(RL,2)
    RpLbCell=repmat({RpLb(ij)},1,size(RpSortCell,2));
    RpUbCell=repmat({RpUb(ij)},1,size(RpSortCell,2));
    IdRpCell=cellfun(@(x,y,z) find(x>=y & x<=z),RpSortCell,RpLbCell,RpUbCell,'UniformOutput',0);
    % Sort first uniform marginal
    [UUcell, IdUCell]=cellfun(@(x,y) sort(x(y)),U1cell,IdRpCell,'UniformOutput',0);


    % Associated second uniform marginal
    VVcell=cellfun(@(x,y) x(y),U2cell,IdRpCell,'UniformOutput',0);

    VVcell=cellfun(@(x,y) x(y),VVcell,IdUCell,'UniformOutput',0);

    % Compute inverse of cdf for each varaible (return levels)
    if strcmpi(margDist,'gp')
    IUUcell=cellfun(@(x,x1,x2,x3,x4,x5,x6) icdf(x5,1-1./(x4(1).*(1./(1-x))),x1(1),x2{1}(x6),x3{1}(x6)),UUcell,epsCell,sigCell,thrCell,SclCell,margDistCell,timeIndexArray,'UniformOutput',0);
    IUUcellend=cellfun(@(x,x1,x2,x3,x4,x5,x6) icdf(x5,1-1./(x4(1).*(1./(1-x))),x1(1),x2{1}(x6),x3{1}(x6)),UUcell,epsCell,sigCell,thrCell,SclCell,margDistCell,timeIndexArray,'UniformOutput',0);

    IVVCell=cellfun(@(x,x1,x2,x3,x4,x5,x6) icdf(x5,1-1./(x4(2).*(1./(1-x))),x1(2),x2{2}(x6),x3{2}(x6)),VVcell,epsCell,sigCell,thrCell,SclCell,margDistCell,timeIndexArray,'UniformOutput',0);
    IVVCellend=cellfun(@(x,x1,x2,x3,x4,x5,x6) icdf(x5,1-1./(x4(2).*(1./(1-x))),x1(2),x2{2}(x6),x3{2}(x6)),VVcell,epsCell,sigCell,thrCell,SclCell,margDistCell,timeIndexArray,'UniformOutput',0);
    else
IUUcell=cellfun(@(x,x1,x2,x3,x4,x5,x6) icdf(x5,1-1./(x4(1).*(1./(1-x))),x1(1),x2{1}(x6),x3{1}(x6)),UUcell,epsCell,sigCell,thrCell,SclCell,margDistCell,timeIndexArray,'UniformOutput',0);
    IUUcellend=cellfun(@(x,x1,x2,x3,x4,x5,x6) icdf(x5,1-1./(x4(1).*(1./(1-x))),x1(1),x2{1}(x6),x3{1}(x6)),UUcell,epsCell,sigCell,thrCell,SclCell,margDistCell,timeIndexArray,'UniformOutput',0);

    IVVCell=cellfun(@(x,x1,x2,x3,x4,x5,x6) icdf(x5,1-1./(x4(1).*(1./(1-x))),x1(2),x2{2}(x6),x3{2}(x6)),VVcell,epsCell,sigCell,thrCell,SclCell,margDistCell,timeIndexArray,'UniformOutput',0);
    IVVCellend=cellfun(@(x,x1,x2,x3,x4,x5,x6) icdf(x5,1-1./(x4(1).*(1./(1-x))),x1(2),x2{2}(x6),x3{2}(x6)),VVcell,epsCell,sigCell,thrCell,SclCell,margDistCell,timeIndexArray,'UniformOutput',0);

    end



    % find non-nan and non-inf values
    IDGoodcell=cellfun(@(x,y) find(~isinf(x)&~isinf(y)&~isnan(x)&~isnan(y)),IUUcell,IVVCell,'UniformOutput',0);
    IDGoodcellend=cellfun(@(x,y) find(~isinf(x)&~isinf(y)&~isnan(x)&~isnan(y)),IUUcellend,IVVCellend,'UniformOutput',0);
    IUUcell=cellfun(@(x,y) x(y),IUUcell,IDGoodcell,'UniformOutput',0);
    IVVCell=cellfun(@(x,y) x(y),IVVCell,IDGoodcell,'UniformOutput',0);
    IUUcellend=cellfun(@(x,y) x(y),IUUcellend,IDGoodcellend,'UniformOutput',0);
    IVVCellend=cellfun(@(x,y) x(y),IVVCellend,IDGoodcellend,'UniformOutput',0);


    % Calculate densities along the probability isoline

    Denscell=cellfun(@(x,y,z,r) copulapdf(r,[x y],z),UUcell,VVcell,PAR,familyCell,'UniformOutput',0);


    % Normalize Densities to the highest density value
    Denscell=cellfun(@(x) x./max(x),Denscell,'UniformOutput',0);
    Denscell=cellfun(@(x,y) x(y),Denscell,IDGoodcell,'UniformOutput',0);

    % color code probability isoline with density level
    colcell= Denscell;
    %assign zero to the third variable (now just works for
    %bivariate)
    ZZcell=cellfun(@(x) zeros(size(x)),IUUcell,'UniformOutput',0);

    % construct x, y, z, and c, usable in surface function
    X=cellfun(@(x1) [x1';x1'],IUUcell,'UniformOutput',0);
    Y=cellfun(@(x1) [x1';x1'],IVVCell,'UniformOutput',0);
    Xend=cellfun(@(x1) [x1';x1'],IUUcellend,'UniformOutput',0);
    Yend=cellfun(@(x1) [x1';x1'],IVVCellend,'UniformOutput',0);
    Z=cellfun(@(x1) [x1';x1'],ZZcell,'UniformOutput',0);
    C=cellfun(@(x1) [x1';x1'],colcell,'UniformOutput',0);


    C=cellfun(@(x1) [x1';x1'],colcell,'UniformOutput',0);
    cellRP{1}=[X,Y,Z,C];
    cellRPend{1}=[Xend,Yend,Z,C];

jk2=[1,length(RpSortCell)];
    for jk=jk2
        if jk==1
           
            [xd,ixd]=unique(X{jk}(:));
            yd=Y{jk}(:);
            yd=yd(ixd);
       [yd,ixd]=unique(yd);       
            xd=xd(ixd);
       

% Define the new area (xrange and yrange) for extrapolation
x_range = get(axxArray(4),'xlim'); % New x-limits for extrapolation
x_extended=linspace(x_range(1),x_range(end),100);
       y_extended=interp1(xd,yd,x_extended,'linear','extrap');   


            plot(axxArray(4),x_extended,y_extended,cRL(ij),'LineWidth',2,'DisplayName',[num2str(RL(ij)),[' - year R.P.']])
            set(axxArray(4),'NextPlot','add')
            % colormap(axxArray(1),'jet')
            % cb = colorbar(axxArray(1),'location','eastoutside');
            % jointExtremes=jointExtremes';
            % text(axxArray(1),mean(IUUcell{jk}),mean(IVVCell{jk}),['Initial Or-JRP'],'Color','k','FontSize',10,'fontweight','bold','HorizontalAlignment','center')
            
            
            % jointExtremes2=unique(cell2mat(jointExtremes),'rows');
            % if jj==length(RL)
            % scatter(axxArray(1),jointExtremes{1}(:,1),jointExtremes{1}(:,2), 'markerfacecolor', 'r');
        elseif jk==length(RpSortCell)
            % surface(axxArray(5),Xend{jk},Yend{jk},Z{jk},C{jk},...
            %     'facecol','no',...
            %     'edgecol','flat',...
            %     'linew',2);
            [xd,ixd]=unique(Xend{jk}(:));
            yd=Yend{jk}(:);
            yd=yd(ixd);
[yd,ixd]=unique(yd);       
            xd=xd(ixd);
            x_range = get(axxArray(5),'xlim'); % New x-limits for extrapolation
x_extended=linspace(x_range(1),x_range(end),100);
       y_extended=interp1(xd,yd,x_extended,'linear','extrap'); 
                        plot(axxArray(5),x_extended,y_extended,cRL(ij),'LineWidth',2,'DisplayName',[num2str(RL(ij)),[' - year R.P.']])
            set(axxArray(5),'NextPlot','add')

         
            % text(axxArray(1),mean(IUUcellend{jk}),mean(IVVCellend{jk}),['Ending Or-JRP'],'Color','k','FontSize',10,'fontweight','bold','HorizontalAlignment','center')
            %scatter(axxArray(1),jointExtremes{end}(:,1),jointExtremes{end}(:,2), 'markerfacecolor', 'k');
           

        end





    end
end



legend(axxArray(4),'show')
legend(axxArray(5),'show')

%prepare output
copulaAnalysis.cellRP=cellRP;
hFig1=get(gcf,'Children');


