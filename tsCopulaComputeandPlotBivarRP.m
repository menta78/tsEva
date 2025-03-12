function  [rpAnalysis,hFig]=tsCopulaComputeandPlotBivarRP(copulaAnalysis,varargin)
%tsCopulaComputeandPlotBivarRP computing and plotting of bivariate return
%period of type "AND"

% [rpAnalysis,hFig]=tsCopulaComputeandPlotBivarRP(copulaAnalysis,axxArray,varargin)
% plots the bivariate joint return period based on "AND" scenario

% input:
% copulaAnalysis          - a variable of type structure provided by
%                           calling tsCopulaExtremes and
%                           tsCopulaCompoundGPDMontecarlo functions first





% output:
% rpAnalysis              - a variable of type structure containing
%                           output of return period analysis

% hFig                    - handle to the drawn figure


%REFERENCES

% [1] Bahmanpour, M.H., Mentaschi, L., Tilloy, A., Vousdoukas, M.,
%     Federico, I., Coppini, G., and Feyen, L., 2025,
%     Transformed-Stationary EVA 2.0: A Generalized Framework for
%     Non-stationary Joint Extreme Analysis (submitted to Hydrology and
%     Earth System Sciences; Feb 2025)

%[2]  Salvadori, G., F. Durante, C. De Michele, M. Bernardi, and L. Petrella
%    (2016), A multivariate copula-based framework for dealing with hazard
%     scenarios and failure probabilities, Water Resour. Res., 52, 3701–3721,
%     doi:10.1002/2015WR017225.

% M.H.Bahmanpour, 2025

colorPool='rgbcmyk';
args.axxArray=gca;

args.RL=[10,50];
args.timeWindowNonStat=[];

args = tsEasyParseNamedArgs(varargin, args);

timeWindowNonStat=(copulaAnalysis.timeWindow)/365.25;

axxArray=args.axxArray;
RL=args.RL;

cRL=colorPool(1:length(RL));

%append RL to copulaAnalysis for the output
rpAnalysis.returnLevels=RL;

% read some information from the copulaAnalysis structure file
jointExtremes=copulaAnalysis.jointExtremes;
margDist = copulaAnalysis.methodology;
Family=copulaAnalysis.copulaParam.family;

if isfield(copulaAnalysis.copulaParam,"rhoMean")
    PAR=copulaAnalysis.copulaParam.rhoMean;
else
    PAR=copulaAnalysis.copulaParam.rho;
end

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
    Scl=ones(1,2);
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

% calculate the AND return period since this type of return period is a
% better representation of "hazard" in the context of hazard mapping

if strcmpi(margDist,'gev')
    rpCellAND=cellfun(@(x,x1) 1./(1+x-x1(:,1)-x1(:,2)),yCell,uCell,'UniformOutput',0);
elseif strcmpi(margDist,'gp')
    if copulaAnalysis.timeVaryingCopula==1
        numJointOrPeaks=cellfun(@(x) length(x),copulaAnalysis.jointExtremes,'UniformOutput',0);
    elseif copulaAnalysis.timeVaryingCopula==0
        numJointOrPeaks={length(copulaAnalysis.jointExtremes)};
    end
    numYearWindow=repmat({timeWindowNonStat},1,size(numJointOrPeaks,2));
    scaling=cellfun(@(x,y) x/y,numYearWindow,numJointOrPeaks,'UniformOutput',0);
    rpCellAND=cellfun(@(x,y,x1) y./(1+x-x1(:,1)-x1(:,2)),yCell,scaling,uCell,'UniformOutput',0);
end
% Sort joint return periods
[rpSortCell, idRPSortCell]=cellfun(@(x) sort(x),rpCellAND,'UniformOutput',0);

% Sort associated uniform marginal probabilities accordingly
[u1Cell,u2Cell]=cellfun(@(x,y) deal(x(y,1), x(y,2)),uCell,idRPSortCell,'UniformOutput',0);

% define return level contours and their acceptable lower and upper bounds
rpLB = RL - 0.005*RL;
rpUB = RL + 0.005*RL;

% convert marginal distribution parameters to cell format
epsCell=repmat({eps},1,size(PAR,2));
sigCell=repmat({sig},1,size(PAR,2));
thrCell=repmat({thr},1,size(PAR,2));
sclCell=repmat({Scl},1,size(PAR,2));
margDistCell=repmat({margDist},1,size(PAR,2));

%to store output information
cellRP=cell(size(RL));

%loop through each return period contour

% Find indices associated with each probability contour
curveBegin={};curveEnd={};
for ij=1:size(RL,2)
    rpLBCell=repmat({rpLB(ij)},1,size(rpSortCell,2));
    rpUBCell=repmat({rpUB(ij)},1,size(rpSortCell,2));
    %find retrun periods in a particular range
    idTargetRPCell=cellfun(@(x,y,z) find(x>=y & x<=z),rpSortCell,rpLBCell,rpUBCell,'UniformOutput',0);

    % Sort first uniform marginal
    [u1TargetRPCell, idU1TargetRPCell]=cellfun(@(x,y) sort(x(y)),u1Cell,idTargetRPCell,'UniformOutput',0);

    % Associated second uniform marginal
    u2TargetRPCell=cellfun(@(x,y) x(y),u2Cell,idTargetRPCell,'UniformOutput',0);

    u2TargetRPCell=cellfun(@(x,y) x(y),u2TargetRPCell,idU1TargetRPCell,'UniformOutput',0);

    % Compute inverse of cdf for each varaible (return levels)

    IUUcell=cellfun(@(x,x1,x2,x3,x4,x5,x6) icdf(x5,1-1./(x4(1).*(1./(1-x))),x1(1),x2{1}(x6),x3{1}(x6)),u1TargetRPCell,epsCell,sigCell,thrCell,sclCell,margDistCell,timeIndexArray,'UniformOutput',0);

    IVVCell=cellfun(@(x,x1,x2,x3,x4,x5,x6) icdf(x5,1-1./(x4(2).*(1./(1-x))),x1(2),x2{2}(x6),x3{2}(x6)),u2TargetRPCell,epsCell,sigCell,thrCell,sclCell,margDistCell,timeIndexArray,'UniformOutput',0);

    % find non-nan and non-inf values
    IDGoodcell=cellfun(@(x,y) find(~isinf(x)&~isinf(y)&~isnan(x)&~isnan(y)),IUUcell,IVVCell,'UniformOutput',0);
    IUUcell=cellfun(@(x,y) x(y),IUUcell,IDGoodcell,'UniformOutput',0);
    IVVCell=cellfun(@(x,y) x(y),IVVCell,IDGoodcell,'UniformOutput',0);

    % construct x, y
    X=cellfun(@(x1) [x1';x1'],IUUcell,'UniformOutput',0);
    Y=cellfun(@(x1) [x1';x1'],IVVCell,'UniformOutput',0);
    
    for jk=1:length(rpSortCell)
        [xd,ixd]=unique(X{jk}(:));
        yd=Y{jk}(:);
        yd=yd(ixd);
        [yd,ixd]=unique(yd);
        xd=xd(ixd);
        if copulaAnalysis.timeVaryingCopula==1
            if jk==1
                % Define the new area (xrange and yrange) for extrapolation
                x_range = get(axxArray(4),'xlim'); % New x-limits for extrapolation
                x_extended=linspace(x_range(1),x_range(end),100);
                y_extended=interp1(xd,yd,x_extended,'linear','extrap');

                plot(axxArray(4),x_extended,y_extended,cRL(ij),'LineWidth',2,'DisplayName',[num2str(RL(ij)),[' - year R.P.']])
                set(axxArray(4),'NextPlot','add')
                curveBegin=[curveBegin,{[x_extended;y_extended]}];
            elseif jk==length(rpSortCell)

                x_range = get(axxArray(5),'xlim'); % New x-limits for extrapolation
                x_extended=linspace(x_range(1),x_range(end),100);
                y_extended=interp1(xd,yd,x_extended,'linear','extrap');
                plot(axxArray(5),x_extended,y_extended,cRL(ij),'LineWidth',2,'DisplayName',[num2str(RL(ij)),[' - year R.P.']])
                set(axxArray(5),'NextPlot','add')
                curveEnd=[curveEnd,{[x_extended;y_extended]}];

            end
        elseif copulaAnalysis.timeVaryingCopula==0
            
                x_range = get(axxArray(4),'xlim'); % New x-limits for extrapolation
                x_extended=linspace(x_range(1),x_range(end),100);
                y_extended=interp1(xd,yd,x_extended,'linear','extrap');

                plot(axxArray(4),x_extended,y_extended,cRL(ij),'LineWidth',2,'DisplayName',[num2str(RL(ij)),[' - year R.P.']])
                set(axxArray(4),'NextPlot','add')
                curveBegin=[curveBegin,{[x_extended;y_extended]}];           

                x_range = get(axxArray(5),'xlim'); % New x-limits for extrapolation
                x_extended=linspace(x_range(1),x_range(end),100);
                y_extended=interp1(xd,yd,x_extended,'linear','extrap');
                plot(axxArray(5),x_extended,y_extended,cRL(ij),'LineWidth',2,'DisplayName',[num2str(RL(ij)),[' - year R.P.']])
                set(axxArray(5),'NextPlot','add')
                curveEnd=[curveEnd,{[x_extended;y_extended]}];
            
        end
    end
end

legend(axxArray(4),'show')
legend(axxArray(5),'show')

%prepare output
rpAnalysis.curveBegin=curveBegin;
rpAnalysis.curveEnd=curveEnd;
hFig=get(gcf,'Children');


