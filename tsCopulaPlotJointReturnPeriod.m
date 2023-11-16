function  tsCopulaPlotJointReturnPeriod(copulaAnalysis,varargin)

% computes the multivariate return period according to AND, OR Scenarios
% other scenarios exist like Kendall and survival kendall but not
% implemented
% see https://doi.org/10.1002/2015WR017225
% Works for Both GPD and GEV (for the GEV, not yet fully implemented as
% the code for selection of joint peaks need to be adjusted accordingly)
% In case of GPD, a scaling factor (npeaks/nyears) needs to be used
% Parts of the code were reworked from MvCAT toolbox https://doi.org/10.1002/2016WR020242
%Bahmanpour, M.H., 2023

args.marginalDistributions = "gp";
args.plotType='AND';
args.PL=[2,5,10,25,50,100,200];

args = tsEasyParseNamedArgs(varargin, args);

margDist = args.marginalDistributions;
plotType= args.plotType;
PL=args.PL;
jMax=copulaAnalysis.jointExtremes; %non-stationary peaks

Family=copulaAnalysis.copulaParam.family;
switch margDist
    case 'gp'
        eps =cellfun(@(x)  x{1}(2).parameters.epsilon,copulaAnalysis.marginalAnalysis,'UniformOutput',1);
        sig=cellfun(@(x)  x{1}(2).parameters.sigma(1),copulaAnalysis.marginalAnalysis,'UniformOutput',1);
        thr=cellfun(@(x)  x{1}(2).parameters.threshold(1),copulaAnalysis.marginalAnalysis,'UniformOutput',1);
        nYear =cellfun(@(x)  (x{1}(2).parameters.timeHorizonEnd-x{1}(2).parameters.timeHorizonStart)/365,copulaAnalysis.marginalAnalysis,'UniformOutput',1);
        nPeak=size(copulaAnalysis.jointExtremes,1);
        Scl2=nPeak/nYear(1); %scaling to account for average number of peaks per year

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


switch Family

    case 'gaussian'
        P = copulacdf('Gaussian',S,PAR);
    case 't'
        P = copulacdf('t',S,PAR(1),PAR(2));
    case 'Clayton'
        P = copulacdf('Clayton',S,PAR);
    case 'Frank'
        P = copulacdf('Frank',S,PAR);
    case 'Gumbel'
        P = copulacdf('Gumbel',S,PAR);

end


switch plotType


    case 'OR'
        cBar = 1 - P; % P is copulaCDF
        RP=1./cBar; %Or scenario

        figure
        jointReturnPeriodPlot(S,RP,Scl2,PL,margDist,eps,sig,thr,Family,PAR,jMax,plotType);
    case 'AND'
        EP1=S(:,1);
        EP2=S(:,2);
        cBar = 1 - EP1-EP2+ P; %survival copula
        RP=1./cBar; %And scenario
        figure
        jointReturnPeriodPlot(S,RP,Scl2,PL,margDist,eps,sig,thr,Family,PAR,jMax,plotType);

end

end


% Compute empirical bivariate probability distribution

function jointReturnPeriodPlot(EP,EBVP,Scl2,PL,margDist,eps,sig,thr,Family,PAR,jMax,plotType)
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

% Loop through probability contours
for j = 1:length(PL)
    % Find indices associated with each probability contour
    ID_Contour = find( P_Sort(:) >= P_LB(j) & P_Sort(:) <= P_UB(j) );

    % Sort first uniform marginal
    [UU, ID_U] = sort( U1(ID_Contour) );
    % Associated second uniform marginal
    VV = U2(ID_Contour);
    VVV = VV(ID_U);

    IUU = icdf( margDist, 1-1./(Scl2*(1./(1-UU))),eps(1),sig(1),thr(1));

    IVVV = icdf( margDist, 1-1./(Scl2*(1./(1-VVV))),eps(2),sig(2),thr(2));
    % Compute inverse of cdf for each varaible


    % Trick: remove infinity if there is any
    IDNAN = find( isinf(IUU) | isinf(IVVV) | isnan(IUU) | isnan(IVVV) );
    IUU(IDNAN) = []; IVVV(IDNAN) = [];

    % Calculate densities along the probability isoline
    switch Family

        case 'gaussian'
            Dens = copulapdf('Gaussian',[UU,VV(ID_U)],PAR);
        case 't'
            Dens = copulapdf('t',[UU,VV(ID_U)],PAR(1,1),PAR(1,2));
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

    % if analytical density does exist, color code probability isoline with density level
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

end