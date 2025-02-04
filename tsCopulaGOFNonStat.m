function [gofStatistics] = tsCopulaGOF(copulaAnalysis,varargin)
%tsCopulaUncertainty estimation of copula goodness-of-fit
% [gofStatistics] = tsCopulaUncertainty(copulaAnalysis,varargin)
%                     returns a variable of type structure containing various parameters
%                     related to the goodness-of-fit of the copula

%goodness-of-fit is always assessed based on a stationary model of copula
%(otherwise, we have accepted the possibility of having different copula
%models in each time window, which is not desirable). 

% A battery of goodnes-of-fit parameters are provided by tsCopulaGOF function.
% These include different assessment of the correlation parameters
% (kendall, pearson and spearman), Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC), 
% and Log-likelihood; moreover, a measure of departure between the fitted and empirical copula is presented based on
% Cramer-Von mises statistic (CvM), including its two variants SnC (based
% on Rosenblatt transformation) and Sn (based on psudo-observations) (see
% ref [2])

% input:
%  copulaAnalysis                           - a variable of type structure provided as the output of tsCopulaCompoundGPD or
%                                              tsCopulaCompoundGPDMontecarlo functions



% output:
%  gofStatistics:                           - A variable of type structure containing:
%                                               xxx                         --
%
%

% M.H.Bahmanpour, 2024

%   References:
%       [1] Mentaschi, L., Vousdoukas, M., Voukouvalas, E., Sartini, L., Feyen, L., Besio, G., and Alfieri, L.:
%           The transformed-stationary approach: a generic and simplified methodology for non-stationary extreme value analysis,
%           Hydrol. Earth Syst. Sci., 20, 3527–3547, https://doi.org/10.5194/hess-20-3527-2016, 2016
%       [2] Genest, C., Rémillard, B., Beaudoin, D., Goodness-of-fit tests for copulas:
%           A review and a power study (Open Access),(2009) Insurance:
%           Mathematics and Economics, 44 (2), pp. 199-213, doi: 10.1016/j.insmatheco.2007.10.005

%       [3] Hofert, M., Kojadinovic, I., Mächler, M. & Yan, J. Elements of Copula Modeling with R (Springer, New York, 2018).

%%%%%%%%%%%%%%%%%%%%%%

% Some parts of the code are reworked from https://github.com/mscavnicky/copula-matlab


% setting the default parameters



args.pValSn=0;
args = tsEasyParseNamedArgs(varargin, args);

pValSn=args.pValSn;
copulaFamily = copulaAnalysis.copulaParam.family;
copulaParam=copulaAnalysis.copulaParam;

jointExtremesResampled=copulaAnalysis.resampleLevel; %from Monte-Carlo simulations

%non-stationary joint extremes
jointExtremes=copulaAnalysis.jointExtremes;
% calculate psuedo-observations
% uSample=tsPseudoObservations(jointExtremes);
uSample= cellfun(@(x) tsPseudoObservations(x),jointExtremes,'UniformOutput',0);
for iFamily=1:length(copulaFamily)
    rhoC=cell(1,length(copulaFamily));
    nuC=cell(1,length(copulaFamily));
if strcmpi(copulaFamily{iFamily},'Gaussian')

    %obtain an estimation of rho based on psudo-observations 
    % rho= copulafit('gaussian', uSample);
  rho=  cellfun(@(x) copulafit('gaussian',x),uSample,'UniformOutput',0);
    rhoC{iFamily}=rho;
    
elseif strcmpi(copulaFamily{iFamily}, 't')

    % [rho,nu]=copulafit('t', uSample);
    % 
    % rhoC{iFamily}=rho;
    % nuC{iFamily}=nu;
elseif strcmpi(copulaFamily{iFamily}, 'Gumbel') || strcmpi(copulaFamily{iFamily}, 'Clayton') || strcmpi(copulaFamily{iFamily}, 'Frank')

    % alpha=copulafit(copulaFamily{iFamily}, uSample);
alpha=  cellfun(@(x) copulafit(copulaFamily{iFamily},x),uSample,'UniformOutput',0);
    
    rhoC{iFamily}=alpha;
end
copulaParam.rho=rhoC;
copulaParam.nu=nuC;


%compute SnC statistic (a variant of CvM Cramer-Von Mises test)
% if strcmpi(copulaParam.family{iFamily},'gaussian') || strcmpi(copulaParam.family{iFamily},'gumbel') || strcmpi(copulaParam.family{iFamily},'clayton')  || strcmpi(copulaParam.family{iFamily},'frank') 
% 
% eSample = tsRosenblattTransform(uSample,'rho',copulaParam.rho{iFamily},'family',copulaParam.family{iFamily});
% elseif strcmpi(copulaParam.family{iFamily},'t')
% eSample = tsRosenblattTransform(uSample,'rho',copulaParam.rho{iFamily},'nu',copulaParam.nu{iFamily},'family',copulaParam.family{iFamily});
% 
% end
% 
% sncSample = sum((tsEmpirical(eSample) - prod(eSample, 2)) .^ 2);

%compute Sn statistic (a variant of CvM Cramer-Von Mises test)

if strcmpi(copulaParam.family{iFamily}, 'gaussian')
    % Y = copulacdf('gaussian', uSample, copulaParam.rho{iFamily});
    Y=cellfun(@(x,y) copulacdf('gaussian',x,y),uSample,copulaParam.rho{iFamily},'UniformOutput',0);
elseif  strcmpi(copulaParam.family{iFamily}, 't')
    % Y = copulacdf('t', uSample, copulaParam.rho{iFamily}, copulaParam.nu{iFamily});
elseif strcmpi(copulaParam.family{iFamily}, 'clayton') || strcmpi(copulaParam.family{iFamily}, 'frank') || strcmpi(copulaParam.family{iFamily}, 'gumbel')
    % Y = copulacdf(copulaParam.family{iFamily}, uSample, copulaParam.rho{iFamily});
       Y=cellfun(@(x,y) copulacdf(copulaParam.family{iFamily},x,y),uSample,copulaParam.rho{iFamily},'UniformOutput',0);
end
% snSample= sum((tsEmpirical(uSample) - Y) .^ 2);
snSample=cellfun(@(x,y) sum((tsEmpirical(x) - y) .^ 2),uSample,Y);
%calculate approximate p-value of sn estimate (details in Hofart book)
% if pValSn
%     N=1000;
%     Pval=tsApproxP(N,copulaParam.family{iFamily},copulaParam.rho{iFamily}, copulaParam.nu{iFamily},snSample,copulaParam.nSeries);
%     gofStatistics(iFamily).Pval=Pval;
% end
gofStatistics(iFamily).snSample=mean(round(snSample*1000)/1000);
gofStatistics(iFamily).copulaFamily=copulaFamily{iFamily};
% corrKendallSample=corr(jointExtremes,'type','Kendall');
corrKendallSample=cellfun(@(x) corr(x,'type','Kendall'),jointExtremes,'UniformOutput',0);
% corrSpearmanSample=corr(jointExtremes,'type','Spearman');
corrSpearmanSample=cellfun(@(x) corr(x,'type','Spearman'),jointExtremes,'UniformOutput',0);

% corrKendallMonte=corr(jointExtremesResampled{iFamily},'type','Kendall');
corrKendallMonte=cellfun(@(x) corr(x,'type','Kendall'),jointExtremesResampled,'UniformOutput',0);

% corrSpearmanMonte=corr(jointExtremesResampled{iFamily},'type','Spearman');
corrSpearmanMonte=cellfun(@(x) corr(x,'type','Spearman'),jointExtremesResampled,'UniformOutput',0);




% kendallDelta=abs(corrKendallSample-corrKendallMonte);
kendallDelta=cellfun(@(x,y) abs(x-y),corrKendallSample,corrKendallMonte,'UniformOutput',0);
% spearmanDelta=abs(corrSpearmanSample-corrSpearmanMonte);
spearmanDelta=cellfun(@(x,y) abs(x-y),corrSpearmanSample,corrSpearmanMonte,'UniformOutput',0);


nVar=copulaAnalysis.copulaParam.nSeries;
idNonDiag=eye(nVar,nVar);
idNonDiag=repmat({idNonDiag},size(kendallDelta));
% [kendallDelta,ia,~]=unique(kendallDelta(~idNonDiag),'stable');
[kendallDelta,ia,~]=cellfun(@(x,y) unique(x(~y),'stable'),kendallDelta,idNonDiag,'UniformOutput',0);

% gofStatistics(iFamily).corrKendallSampleDelta=round(mean(kendallDelta)*1000)/1000;
gofStatistics(iFamily).corrKendallSampleDelta=mean(cellfun(@(x) round(mean(x)*1000)/1000,kendallDelta,'UniformOutput',1))
% spearmanDelta=spearmanDelta(~idNonDiag);
spearmanDelta=cellfun(@(x,y) x(~y),spearmanDelta,idNonDiag,'UniformOutput',0);
% spearmanDelta=spearmanDelta(ia);
spearmanDelta=cellfun(@(x,y) x(y),spearmanDelta,ia,'UniformOutput',0);
% gofStatistics(iFamily).corrSpearmanSampleDelta=round(mean(spearmanDelta)*1000)/1000;
gofStatistics(iFamily).corrSpearmanSampleDelta=mean(cellfun(@(x) round(mean(x)*1000)/1000,spearmanDelta))


corrSpearmanSamplex=cellfun(@(x,y) x(~y),corrSpearmanSample,idNonDiag,'UniformOutput',0);

corrSpearmanSamplex=cellfun(@(x,y) x(y),corrSpearmanSamplex,ia,'UniformOutput',0);

corrSpearmanMontex=cellfun(@(x,y) x(~y),corrSpearmanMonte,idNonDiag,'UniformOutput',0);

corrSpearmanMontex=cellfun(@(x,y) x(y),corrSpearmanMontex,ia,'UniformOutput',0);
gofStatistics(iFamily).corrSpearmanSamplex=corrSpearmanSamplex;
gofStatistics(iFamily).corrSpearmanMontex=corrSpearmanMontex;

end


fields = fieldnames(gofStatistics);
igd=find(cellfun(@(x) isempty(regexp(x,'copulaFamily')),fields));
imnst=[];
for ifields=igd'
    if  strcmpi(fields{ifields},'snsample') || strcmpi(fields{ifields},'aicsample') || strcmpi(fields{ifields},'bicsample') || strcmpi(fields{ifields},'sncsample')

        ser0=extractfield(gofStatistics,fields{ifields});
        [mser0]=min(ser0);
        imns=find(ser0==mser0);
        imnst=[imnst,imns];
    elseif strcmpi(fields{ifields},'llsample')
        ser0=extractfield(gofStatistics,fields{ifields});
        [maser0]=max(ser0);
        imns=find(ser0==maser0);
        imnst=[imnst,imns];
    elseif strcmpi(fields{ifields},'corrkendallsampledelta')  || strcmpi(fields{ifields},'corrspearmansampledelta') || strcmpi(fields{ifields},'corrpearsonsampledelta')
        ser0=extractfield(gofStatistics,fields{ifields});
        if copulaAnalysis.copulaParam.nSeries==3
        ser00=ser0;
        elseif copulaAnalysis.copulaParam.nSeries==2
            ser00=ser0;
        end
        [mser02]=min((ser00));
         imns=find(ser00==mser02);
        imnst=[imnst,imns];
    end
end


[n,bin] = hist(imnst,unique(imnst));
bin=bin(find(n));
n=n(find(n));
[~,idx] = sort(-n);
 
res=cellfun(@(x,y) [x,' with ',num2str(y),' test passed'],copulaFamily(bin(idx)),num2cell(n(idx)),'UniformOutput',0);

disp(['GOF results:'])
disp([char(res')])



end





