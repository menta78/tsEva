function  [copulaAnalysis,hFig1]=tsCopulaJointReturnPeriod(copulaAnalysis,varargin)
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

args.xlbl = 'X';
args.ylbl = 'Y';
args.marginalDistributions = "gp";
args.plotType='AND';
args.RL=[2,20,50,100,200];
args.numberofverticalpanels=2;
args.numberofhorizontalpanels=3;

args = tsEasyParseNamedArgs(varargin, args);

xlbl = args.xlbl;
ylbl = args.ylbl;
margDist = args.marginalDistributions;
plotType= args.plotType;
RL=args.RL;
numberofverticalpanels=args.numberofverticalpanels;
numberofhorizontalpanels=args.numberofhorizontalpanels;

%append RL to copulaAnalysis for the output
copulaAnalysis.RL=RL;

% read some information from the copulaAnalysis structure file
jointExtremes=copulaAnalysis.jointExtremes;
Family=copulaAnalysis.copulaParam.family;
PAR=copulaAnalysis.copulaParam.rho;
if strcmpi(Family,'t')
    PAR2=copulaAnalysis.copulaParam.nu;
end
switch margDist
    case 'gp'
        eps =cellfun(@(x)  x{1}(2).parameters.epsilon,copulaAnalysis.marginalAnalysis,'UniformOutput',1);
        sig=cellfun(@(x)  x{1}(2).parameters.sigma(1),copulaAnalysis.marginalAnalysis,'UniformOutput',1);
        thr=cellfun(@(x)  x{1}(2).parameters.threshold(1),copulaAnalysis.marginalAnalysis,'UniformOutput',1);
        nYear =cellfun(@(x)  (x{1}(2).parameters.timeHorizonEnd-x{1}(2).parameters.timeHorizonStart)/365,copulaAnalysis.marginalAnalysis,'UniformOutput',1);
        nPeak=cellfun(@(x)  length(x{1}(2).objs.peakIndexes),copulaAnalysis.marginalAnalysis,'UniformOutput',1);
        %scaling to account for average number of peaks per year- needed only when GPD return levels are to be calculated
        Scl=nPeak./nYear;

    case 'gev'
        eps=cellfun(@(x)  x{1}(1).parameters.epsilon,copulaAnalysis.marginalAnalysis,'UniformOutput',1);
        sig=cellfun(@(x)  x{1}(1).parameters.sigma(1),copulaAnalysis.marginalAnalysis,'UniformOutput',1);
        thr=cellfun(@(x)  x{1}(1).parameters.mu(1),copulaAnalysis.marginalAnalysis,'UniformOutput',1);
        Scl=1;
end

%define a uniform square (or hypercube) over which analytical copulaCDF would be calculated
% for bivariate Archimedean copula, then u must be an n-by-2 matrix

x = [linspace(1e-5,.95,400) linspace(0.95+1e-5,1-1e-5,500)];
[xx,yy] = meshgrid(x,x);
U = [xx(:),yy(:)];

if ~iscell(PAR)
    %concerning a stationary copula

    if strcmpi(Family,'gaussian') || strcmpi(Family,'clayton') || strcmpi(Family,'frank') || strcmpi(Family,'gumbel')
        y = copulacdf(Family,U,PAR);
    elseif strcmpi(Family,'t')
        y = copulacdf(Family,U,PAR,PAR2);
    end

    switch plotType
        % by definition, joint return periods are normally calculated using
        % an "or" or "and" scenario with different formulas
        case 'OR'
            cBar = 1 - y;
            RP=1./cBar;
        case 'AND'
            cBar = 1 - U(:,1) - U(:,2) + y; %survival copula
            RP=1./cBar;
    end

    % Sort joint return periods
    [RpSort, IdRpSrt] = sort(RP);

    % sort associated uniform marginal probabilities accordingly
    U1  = U(IdRpSrt, 1);
    U2  = U(IdRpSrt, 2);

    % define return level contours and their acceptable lower and upper bounds
    RpLb = RL - 0.005*RL;
    RpUb = RL + 0.005*RL;

    %%define figure properties
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

    cellRP=cell(size(RL));
    % Loop through probability (or return period) contours
    for jj = 1:length(RL)
        % Find indices associated with each probability contour
        IdRp = find( RpSort(:) >= RpLb(jj) & RpSort(:) <= RpUb(jj) );

        % Select first uniform marginals that fall withing selected
        % return level and sort it accordingly
        [UU, IdU] = sort( U1(IdRp) );
        % Associated second uniform marginal
        VV = U2(IdRp);
        VV = VV(IdU);
        % calculate inverse cdf for each varaible using marginal information
        IUU = icdf( margDist, 1-1./(Scl(1)*(1./(1-UU))),eps(1),sig(1),thr(1));

        IVV = icdf( margDist, 1-1./(Scl(2)*(1./(1-VV))),eps(2),sig(2),thr(2));

        % remove infinity or nan if there is any
        IDNAN = find( isinf(IUU) | isinf(IVV) | isnan(IUU) | isnan(IVV) );
        IUU(IDNAN) = []; IVV(IDNAN) = [];

        % Calculate densities along the probability (or return level) isoline

        if strcmpi(Family,'Gaussian') || strcmpi(Family,'clayton') || strcmpi(Family,'frank') || strcmpi(Family,'gumbel')
            Dens = copulapdf(Family,[UU,VV],PAR);
        elseif strcmpi(Family,'t')
            Dens = copulapdf('t',[UU,VV],PAR,PAR2);
        end

        % Normalize Densities to the highest density value
        Dens = Dens / max(Dens);
        Dens(IDNAN) = [];

        % Color code probability (or return level) isoline with density level
        col = Dens;

        % ZZ could be used if three variates are included; in case of bivariate, ZZ
        % must be set to zero
        ZZ = zeros(size(IUU));

        %store information for the output
        cellRP{jj}=[IUU,IVV,ZZ,col];
        % generate surface plot
        h(jj)=surface([IUU';IUU'],[IVV';IVV'],[ZZ';ZZ'],[col';col'],...
            'facecol','no',...
            'edgecol','flat',...
            'linew',2);
        colormap(jet)
        %generate color bar and its label; also write title of the
        %plot
        if jj == 1
            cb = colorbar('location','south');
            set(cb, 'xlim', [0 1]);
            set(cb,'position',[cb.Position(1) cb.Position(2) 0.6*cb.Position(3) 0.5*cb.Position(4)])
            cb.Label.String='Copula pdf';

            t = title(strcat('\color{red}',{[upper(Family),' Copula'];['Joint Return Period - ', plotType,' Scenario']}));
            set(t, 'units', 'normalized', 'horizontalAlignment', 'center','fontname','times','fontweight','bold','fontsize',12);
        end
        % annotate each return level contour level
        text(max(IUU),min(IVV),num2str(RL(jj)),'Color','red','FontSize',10,'fontweight','bold')

    end
    box on;
    hold on
    h(end+1)=plot(jointExtremes(:,1),jointExtremes(:,2),'r*');
    xlabel(xlbl)
    ylabel(ylbl)
    ltt=legend(h(end),'Joint peaks');
    set(ltt,'Position',[0.4*ltt.Position(1) ltt.Position(2) ltt.Position(3) ltt.Position(4)])
    %prepare output
    copulaAnalysis.cellRP=cellRP;
    hFig1=get(gcf,'Children');


elseif iscell(PAR) % time-varying copula intended

    uCell=repmat({U},1,size(PAR,2));
    familyCell=repmat({Family},1,size(PAR,2));
    if strcmpi(Family,'Gaussian') || strcmpi(Family,'Frank') || strcmpi(Family,'Gumbel') || strcmpi(Family,'Clayton')
        yCell=cellfun(@(x,y,z) copulacdf(z,x,y),uCell,PAR,familyCell,'UniformOutput',0);
    elseif strcmpi(Family,'t')
        yCell=cellfun(@(x,y,z) copulacdf('t',x,y,z),uCell,PAR,PAR2,'UniformOutput',0);
    end

    inputtimestampsWindowCell=copulaAnalysis.copulaParam.inputtimestampsWindowCell;

    switch plotType

        case 'OR'

            rPCell=cellfun(@(x) 1./(1-x),yCell,'UniformOutput',0);

        case 'AND'

            rPCell=cellfun(@(x,y) 1./(1-y(:,1)-y(:,2)+x),yCell,uCell,'UniformOutput',0);
    end

    axxcell={}; %cell array containing all axes

    nSubPlot=size(copulaAnalysis.jointExtremes,2);
    if mod(nSubPlot,numberofverticalpanels*numberofhorizontalpanels)==0
        nFig=nSubPlot/(numberofverticalpanels*numberofhorizontalpanels);
    else
        nFig=fix(nSubPlot/(numberofverticalpanels*numberofhorizontalpanels))+1;
    end

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
                axxcell=[axxcell,axx];

                axn=axn+1;

            end
        end

    end

    % Sort joint return periods
    [RpSortCell, IdRpSrtCell]=cellfun(@(x) sort(x),rPCell,'UniformOutput',0);

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
    for jj = 1:length(RL)
        % Find indices associated with each probability contour

        RpLbCell=repmat({RpLb(jj)},1,size(RpSortCell,2));
        RpUbCell=repmat({RpUb(jj)},1,size(RpSortCell,2));
        IdRpCell=cellfun(@(x,y,z) find(x>=y & x<=z),RpSortCell,RpLbCell,RpUbCell,'UniformOutput',0);

        % Sort first uniform marginal
        [UUcell, IdUCell]=cellfun(@(x,y) sort(x(y)),U1cell,IdRpCell,'UniformOutput',0);


        % Associated second uniform marginal
        VVcell=cellfun(@(x,y) x(y),U2cell,IdRpCell,'UniformOutput',0);

        VVcell=cellfun(@(x,y) x(y),VVcell,IdUCell,'UniformOutput',0);

        % Compute inverse of cdf for each varaible (return levels)
        IUUcell=cellfun(@(x,x1,x2,x3,x4,x5) icdf(x5,1-1./(x4(1).*(1./(1-x))),x1(1),x2(1),x3(1)),UUcell,epsCell,sigCell,thrCell,SclCell,margDistCell,'UniformOutput',0);
        IVVCell=cellfun(@(x,x1,x2,x3,x4,x5) icdf(x5,1-1./(x4(2).*(1./(1-x))),x1(2),x2(2),x3(2)),VVcell,epsCell,sigCell,thrCell,SclCell,margDistCell,'UniformOutput',0);


        % find non-nan and non-inf values
        IDGoodcell=cellfun(@(x,y) find(~isinf(x)&~isinf(y)&~isnan(x)&~isnan(y)),IUUcell,IVVCell,'UniformOutput',0);
        IUUcell=cellfun(@(x,y) x(y),IUUcell,IDGoodcell,'UniformOutput',0);
        IVVCell=cellfun(@(x,y) x(y),IVVCell,IDGoodcell,'UniformOutput',0);

        % Calculate densities along the probability isoline
        if strcmpi(Family,'Gaussian') || strcmpi(Family,'Frank') || strcmpi(Family,'Gumbel') || strcmpi(Family,'Clayton')
            Denscell=cellfun(@(x,y,z,r) copulapdf(r,[x y],z),UUcell,VVcell,PAR,familyCell,'UniformOutput',0);
        elseif strcmpi(Family,'t')
            Denscell=cellfun(@(x,y,z,z1) copulapdf('t',[x y],z,z1),UUcell,VVcell,PAR,PAR2,'UniformOutput',0);
        end

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
        Z=cellfun(@(x1) [x1';x1'],ZZcell,'UniformOutput',0);
        C=cellfun(@(x1) [x1';x1'],colcell,'UniformOutput',0);
        cellRP{jj}=[X,Y,Z,C];
        for jk=1:size(RpSortCell,2)
            surface(axxcell(jk),X{jk},Y{jk},Z{jk},C{jk},...
                'facecol','no',...
                'edgecol','flat',...
                'linew',2);
            if jj==1
                set(axxcell(jk),'NextPlot','add')
                colormap(axxcell(jk),'jet')
                cb = colorbar(axxcell(jk),'location','south');
                if strcmpi(Family,'Gaussian')
                    nr=round(PAR{jk}(2)*100)/100;
                elseif strcmpi(Family,'Frank') || strcmpi(Family,'Gumbel') || strcmpi(Family,'Clayton')
                    nr=round(PAR{jk}*100)/100;
                elseif strcmpi(family,'t')
                    nr=round(PAR{jk}(2)*100)/100;
                    nr2=round(PAR2{jk}(2)*100)/100;
                end
                t1x=datestr(inputtimestampsWindowCell{jk}(1),'yyyy');
                t2x=datestr(inputtimestampsWindowCell{jk}(end),'yyyy');
                if strcmpi(Family,'Gaussian') || strcmpi(Family,'Frank') || strcmpi(Family,'Gumbel') || strcmpi(Family,'Clayton')
                    ht=title(axxcell(jk),{[t1x,' - ',t2x];['Rho = ',num2str(nr)]});
                elseif strcmpi(Family,'t')
                    ht=title(axxcell(jk),{[t1x,' - ',t2x];['Rho = ',num2str(nr)];['Nu = ',num2str(nr2)]});
                end
                ht.VerticalAlignment='top';
            end
            if jj==length(RL)
                scatter(axxcell(jk),jointExtremes{jk}(:,1),jointExtremes{jk}(:,2), 'markerfacecolor', 'r');
            end
            if mod(jj,2)==0
                text(axxcell(jk),max(IUUcell{jk}),min(IVVCell{jk}),num2str(RL(jj)),'Color','red','FontSize',10,'fontweight','bold','HorizontalAlignment','center')
            end
        end
    end
    %prepare output
    copulaAnalysis.cellRP=cellRP;
    hFig1=get(gcf,'Children');
end
end
