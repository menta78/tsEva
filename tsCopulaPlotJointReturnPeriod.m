function  tsCopulaPlotJointReturnPeriod(copulaAnalysis,varargin)

% plots the multivariate return period according to AND/OR Scenarios
% see https://doi.org/10.1002/2015WR017225
% Parts of the code were reworked from MvCAT toolbox https://doi.org/10.1002/2016WR020242
%
%Bahmanpour, M.H., 2023

args.marginalDistributions = "gp";
args.plotType='AND';
args.PL=[2,5,10,25,50,100,200];

args = tsEasyParseNamedArgs(varargin, args);

margDist = args.marginalDistributions;
plotType= args.plotType;
PL=args.PL;

jMax=copulaAnalysis.jointExtremes; %extreme compound events sampled previously
Family=copulaAnalysis.copulaParam.family;

switch margDist
    case 'gp'
        eps =cellfun(@(x)  x{1}(2).parameters.epsilon,copulaAnalysis.marginalAnalysis,'UniformOutput',1);
        sig=cellfun(@(x)  x{1}(2).parameters.sigma(1),copulaAnalysis.marginalAnalysis,'UniformOutput',1);
        thr=cellfun(@(x)  x{1}(2).parameters.threshold(1),copulaAnalysis.marginalAnalysis,'UniformOutput',1);
        nYear =cellfun(@(x)  (x{1}(2).parameters.timeHorizonEnd-x{1}(2).parameters.timeHorizonStart)/365,copulaAnalysis.marginalAnalysis,'UniformOutput',1);
        nPeak=cellfun(@(x)  length(x{1}(2).objs.peakIndexes),copulaAnalysis.marginalAnalysis,'UniformOutput',1);
        Scl2=nPeak./nYear; %scaling to account for average number of peaks per year- needed when GPD return levels are to be calculated

    case 'gev'
        eps=cellfun(@(x)  x{1}(1).parameters.epsilon,copulaAnalysis.marginalAnalysis,'UniformOutput',1);
        sig=cellfun(@(x)  x{1}(1).parameters.sigma(1),copulaAnalysis.marginalAnalysis,'UniformOutput',1);
        thr=cellfun(@(x)  x{1}(1).parameters.mu(1),copulaAnalysis.marginalAnalysis,'UniformOutput',1);
        Scl2=1; %annual maxima
end

PAR=copulaAnalysis.copulaParam.rho;


x = [linspace(1e-5,.95,400) linspace(0.95+1e-5,1-1e-5,500)];
[xx,yy] = meshgrid(x,x);
S = [xx(:),yy(:)];

if ~iscell(PAR)  %stationary copula is intended
    switch Family

        case 'Gaussian'
            P = copulacdf('Gaussian',S,PAR);
        case 't'
            PAR2=copulaAnalysis.copulaParam.nu;
            P = copulacdf('t',S,PAR,PAR2);
        case 'Clayton'
            P = copulacdf('Clayton',S,PAR);
        case 'Frank'
            P = copulacdf('Frank',S,PAR);
        case 'Gumbel'
            P = copulacdf('Gumbel',S,PAR);

    end
elseif iscell(PAR) % time-varying copula intended

    sCell=repmat({S},1,size(PAR,2));
    switch Family

        case 'Gaussian'
            % P = copulacdf('Gaussian',S,PAR);
            pCell=cellfun(@(x,y) copulacdf('Gaussian',x,y),sCell,PAR,'UniformOutput',0);
        case 't'
            %  P = copulacdf('t',S,PAR(1),PAR(2));
            PAR2=copulaAnalysis.copulaParam.nu;
           pCell=cellfun(@(x,y,z) copulacdf('t',x,y,z),sCell,PAR,PAR2,'UniformOutput',0);
            
        case 'Clayton'
            pCell=cellfun(@(x,y) copulacdf('Clayton',x,y),sCell,PAR,'UniformOutput',0);
        case 'Frank'
            pCell=cellfun(@(x,y) copulacdf('Frank',x,y),sCell,PAR,'UniformOutput',0);
        case 'Gumbel'
            pCell=cellfun(@(x,y) copulacdf('Gumbel',x,y),sCell,PAR,'UniformOutput',0);

    end
end

if ~iscell(PAR)
    switch plotType

        case 'OR'
            cBar = 1 - P; % P is copulaCDF
            RP=1./cBar; %OR scenario

            figure
            jointReturnPeriodPlot(S,RP,Scl2,PL,margDist,eps,sig,thr,Family,PAR,jMax,plotType,copulaAnalysis);
        case 'AND'
            EP1=S(:,1);
            EP2=S(:,2);
            cBar = 1 - EP1-EP2+ P; %survival copula
            RP=1./cBar; %AND scenario
            figure
            jointReturnPeriodPlot(S,RP,Scl2,PL,margDist,eps,sig,thr,Family,PAR,jMax,plotType,copulaAnalysis);

    end
elseif iscell(PAR)
inputtimestampsWindowCell=copulaAnalysis.copulaParam.inputtimestampsWindowCell;

    switch plotType

        case 'OR'
          
            rPCell=cellfun(@(x) 1./(1-x),pCell,'UniformOutput',0);
           
            jointReturnPeriodPlot(sCell,rPCell,Scl2,PL,margDist,eps,sig,thr,Family,PAR,jMax,plotType,copulaAnalysis,'timeStampWindowCell',inputtimestampsWindowCell);
        case 'AND'
           
            rPCell=cellfun(@(x,y) 1./(1-y(:,1)-y(:,2)+x),pCell,sCell,'UniformOutput',0);
            
            jointReturnPeriodPlot(sCell,rPCell,Scl2,PL,margDist,eps,sig,thr,Family,PAR,jMax,plotType,copulaAnalysis,'timeStampWindowCell',inputtimestampsWindowCell);

    end

end
end


function jointReturnPeriodPlot(EP,EBVP,Scl2,PL,margDist,eps,sig,thr,Family,PAR,jMax,plotType,copulaAnalysis,varargin)

args.timeStampWindowCell = {};

args = tsEasyParseNamedArgs(varargin, args);

inputtimestampsWindowCell = args.timeStampWindowCell;

if ~iscell(EBVP)
    hold on;
    % Sort data based on joint probability
    [P_Sort, ID_Pr] = sort(EBVP);
    % Associated uniform marginal probabilities
    U1  = EP(ID_Pr, 1);
    U2  = EP(ID_Pr, 2);

    % Probability contours and their acceptable lower and upper bounds

    P_LB = PL - 0.005*PL; P_UB = PL + 0.005*PL;

    % Define size of text
    SIZE = 12;

    % Loop through probability (or return period) contours
    for j = 1:length(PL)
        % Find indices associated with each probability contour
        ID_Contour = find( P_Sort(:) >= P_LB(j) & P_Sort(:) <= P_UB(j) );

        % Sort first uniform marginal
        [UU, ID_U] = sort( U1(ID_Contour) );
        % Associated second uniform marginal
        VV = U2(ID_Contour);
        VVV = VV(ID_U);

        IUU = icdf( margDist, 1-1./(Scl2(1)*(1./(1-UU))),eps(1),sig(1),thr(1));

        IVVV = icdf( margDist, 1-1./(Scl2(2)*(1./(1-VVV))),eps(2),sig(2),thr(2));
        % Compute inverse of cdf for each varaible (with appropriate
        % scaling for GPD return levels to account for average number of peaks per year)

        % Trick: remove infinity if there is any
        IDNAN = find( isinf(IUU) | isinf(IVVV) | isnan(IUU) | isnan(IVVV) );
        IUU(IDNAN) = []; IVVV(IDNAN) = [];

        % Calculate densities along the probability isoline
        switch Family

            case 'Gaussian'
                Dens = copulapdf('Gaussian',[UU,VV(ID_U)],PAR);
            case 't'
                PAR2=copulaAnalysis.copulaParam.nu;
                Dens = copulapdf('t',[UU,VV(ID_U)],PAR,PAR2);
            case 'Clayton'
                Dens = copulapdf('Clayton',[UU,VV(ID_U)],PAR);
            case 'Frank'
                Dens = copulapdf('Frank',[UU,VV(ID_U)],PAR);
            case 'Gumbel'
                Dens = copulapdf('Gumbel',[UU,VV(ID_U)],PAR);

        end
        % Normalize Densities to the highest density value
        Dens = Dens / max(Dens);
        Dens(IDNAN) = [];

        % Color code probability isoline with density level
        ZZ = zeros(size(IUU));
        col = Dens;  % This is the color, vary with x in this case.
        h(j)=surface([IUU';IUU'],[IVVV';IVVV'],[ZZ';ZZ'],[col';col'],...
            'facecol','no',...
            'edgecol','flat',...
            'linew',2);
        colormap(jet)
        if j == 1
            cb = colorbar('location','south');
            set(cb, 'xlim', [0 1]);
            set(cb,'position',[cb.Position(1) cb.Position(2) 0.6*cb.Position(3) 0.5*cb.Position(4)])
            cb.Label.String='Copula pdf'
        end

        text(max(IUU),min(IVVV),num2str(PL(j)),'Color','red','FontSize',10,'fontweight','bold')

        if j == 1 % Write the title once
            t = title(strcat('\color{red}',{[upper(Family),' Copula'];['Joint Return Period - ', plotType,' Scenario']}));
            set(t, 'units', 'normalized', 'horizontalAlignment', 'center','fontname','times','fontweight','bold','fontsize',SIZE);

        end

    end
    box on; %axis square
    hold on
    h(end+1)=plot(jMax(:,1),jMax(:,2),'r*');
    xlabel('Variable 1')
    ylabel('Variable 2')
    legend(h(end),'Joint peaks')

elseif iscell(EBVP) %time-varying joint return periods

    axxcell={}; %cell array containing all axes

    nSubPlot=size(EP,2);
    if mod(nSubPlot,8)==0
        nFig=nSubPlot/8;
    else
        nFig=fix(nSubPlot/8)+1;
    end

    axn=1;

    for ij=1:nFig

        rt=1;
        b0=27;
        l0=21;
        spMan = tsLcSubplotManager(b0, l0, 'CellXSize', round(21*rt), 'CellYSize', round(27*rt), 'gap', [0.03 0.01]);

        spMan.initFigure;

        numberofverticalpanels=4;
        numberofhorizontalpanels=2;

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
                axxcell=[axxcell,axx];

                axn=axn+1;

            end
        end

    end

    % Sort data based on joint probability
    
    [P_Sortcell, ID_Prcell]=cellfun(@(x) sort(x),EBVP,'UniformOutput',0);

    % Associated uniform marginal probabilities
    
    U1cell=cellfun(@(x,y) x(y,1),EP,ID_Prcell,'UniformOutput',0);
    U2cell=cellfun(@(x,y) x(y,2),EP,ID_Prcell,'UniformOutput',0);

    % Probability contours (or return periods) and their acceptable lower and upper bounds

    P_LB = PL - 0.005*PL; P_UB = PL + 0.005*PL;

    % Define size of text
    SIZE = 12;

    epscell=repmat({eps},1,size(PAR,2));
    sigcell=repmat({sig},1,size(PAR,2));
    thrcell=repmat({thr},1,size(PAR,2));
    Scl2cell=repmat({Scl2},1,size(PAR,2));
    margDistcell=repmat({margDist},1,size(PAR,2));

    for jk=1:size(P_Sortcell,2)
       
        for j = 1:length(PL)
            % Find indices associated with each probability contour
           
            P_LBcell=repmat({P_LB(j)},1,size(P_Sortcell,2));
            P_UBcell=repmat({P_UB(j)},1,size(P_Sortcell,2));
            ID_Contourcell=cellfun(@(x,y,z) find(x>=y & x<=z),P_Sortcell,P_LBcell,P_UBcell,'UniformOutput',0);

            % Sort first uniform marginal
            [UUcell, ID_Ucell]=cellfun(@(x,y) sort(x(y)),U1cell,ID_Contourcell,'UniformOutput',0);

         
            % Associated second uniform marginal
            VVcell=cellfun(@(x,y) x(y),U2cell,ID_Contourcell,'UniformOutput',0);
          
            VVVcell=cellfun(@(x,y) x(y),VVcell,ID_Ucell,'UniformOutput',0);

            % Compute inverse of cdf for each varaible (return levels)
            
            IUUcell=cellfun(@(x,x1,x2,x3,x4,x5) icdf(x5,1-1./(x4(1).*(1./(1-x))),x1(1),x2(1),x3(1)),UUcell,epscell,sigcell,thrcell,Scl2cell,margDistcell,'UniformOutput',0);
            IVVVcell=cellfun(@(x,x1,x2,x3,x4,x5) icdf(x5,1-1./(x4(2).*(1./(1-x))),x1(2),x2(2),x3(2)),VVVcell,epscell,sigcell,thrcell,Scl2cell,margDistcell,'UniformOutput',0);

          
            % Trick: find non-nan and non-inf values
            IDGoodcell=cellfun(@(x,y) find(~isinf(x)&~isinf(y)&~isnan(x)&~isnan(y)),IUUcell,IVVVcell,'UniformOutput',0);
            IUUcell=cellfun(@(x,y) x(y),IUUcell,IDGoodcell,'UniformOutput',0);
            IVVVcell=cellfun(@(x,y) x(y),IVVVcell,IDGoodcell,'UniformOutput',0);   

            % Calculate densities along the probability isoline
            switch Family

                case 'Gaussian'
                    % Dens = copulapdf('Gaussian',[UU,VV(ID_U)],PAR);
                    Denscell=cellfun(@(x,y,z) copulapdf('Gaussian',[x y],z),UUcell,VVVcell,PAR,'UniformOutput',0);
                case 't'
                    
                    PAR2=copulaAnalysis.copulaParam.nu;
                    Denscell=cellfun(@(x,y,z,z1) copulapdf('t',[x y],z,z1),UUcell,VVVcell,PAR,PAR2,'UniformOutput',0);
                case 'Clayton'
                    Denscell=cellfun(@(x,y,z) copulapdf('Clayton',[x y],z),UUcell,VVVcell,PAR,'UniformOutput',0);
                case 'Frank'
                    Denscell=cellfun(@(x,y,z) copulapdf('Frank',[x y],z),UUcell,VVVcell,PAR,'UniformOutput',0);
                case 'Gumbel'
                    Denscell=cellfun(@(x,y,z) copulapdf('Gumbel',[x y],z),UUcell,VVVcell,PAR,'UniformOutput',0);

            end

            % Normalize Densities to the highest density value
            Denscell=cellfun(@(x) x./max(x),Denscell,'UniformOutput',0);
      
            Denscell=cellfun(@(x,y) x(y),Denscell,IDGoodcell,'UniformOutput',0);
            
            ZZcell=cellfun(@(x) zeros(size(x)),IUUcell,'UniformOutput',0);

            % color code probability isoline with density level
         
            colcell= Denscell;
           
            % construct x, y, z, and c, usable in surface function
            X=cellfun(@(x1) [x1';x1'],IUUcell,'UniformOutput',0);
            Y=cellfun(@(x1) [x1';x1'],IVVVcell,'UniformOutput',0);
            Z=cellfun(@(x1) [x1';x1'],ZZcell,'UniformOutput',0);
            C=cellfun(@(x1) [x1';x1'],colcell,'UniformOutput',0);
           
            surface(axxcell(jk),X{jk},Y{jk},Z{jk},C{jk},...
                'facecol','no',...
                'edgecol','flat',...
                'linew',2);

            set(axxcell(jk),'NextPlot','add')
            if mod(j,2)==0
            text(axxcell(jk),max(IUUcell{jk}),min(IVVVcell{jk}),num2str(PL(j)),'Color','red','FontSize',10,'fontweight','bold','HorizontalAlignment','center')
            end
            if j==1
                colormap(axxcell(jk),'jet')
                cb = colorbar(axxcell(jk),'location','south');
                if strcmp(copulaAnalysis.copulaParam.family,'Gaussian') ||strcmp(copulaAnalysis.copulaParam.family,'t')
                    nr=round(PAR{jk}(2)*1000)/1000;
                elseif strcmp(copulaAnalysis.copulaParam.family,'Frank') || strcmp(copulaAnalysis.copulaParam.family,'Clayton') ||strcmp(copulaAnalysis.copulaParam.family,'Gumbel')
                    nr=round(PAR{jk}(1));
                end
              
                t1x=datestr(inputtimestampsWindowCell{jk}(1),'yyyy');
                t2x=datestr(inputtimestampsWindowCell{jk}(end),'yyyy');
                ht=title(axxcell(jk),{[t1x,' - ',t2x];['Rho = ',num2str(nr)]});
                ht.VerticalAlignment='top';
            end
        end
        scatter(axxcell(jk),jMax{jk}(:,1),jMax{jk}(:,2), 'markerfacecolor', 'r');
      
    end
   
end
end