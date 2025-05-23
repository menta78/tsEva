function [gofStatistics] = tsCopulaGOFNonStat(copulaAnalysis, monteCarloAnalysis, varargin)
%tsCopulaGOFNonStat estimation of copula goodness-of-fit and other battery
%of statistics

% [gofStatistics] = tsCopulaGOF(copulaAnalysis,varargin)
%                     returns a variable of type structure containing various parameters
%                     related to the goodness-of-fit of the copula


% To evaluate the goodness-of-fit (GOF) of the copula model, a
% multi-parameter approach was employed by analyzing a set of statistics
% that quantify the similarity between the fitted distribution and the
% empirical data. Specifically, the following statistics were considered:
% 	Cramér-von Mises statistic. The statistic Sn serves as a proxy for the
% 	distance
% between the empirical and theoretical distributions in probability space.
% In this study, we applied the rank-based version of the Cramér-von Mises
% statistic, where the ranks of Cn and Cθ are compared. For non-stationary
% distributions, Sn is estimated separately over different time windows,
% and the results are averaged, to provide the mean Cramér-von Mises
% statistic.
%
% 	For each bivariate sub-distribution, the goodness-of-fit was evaluated
% 	by comparing the correlation structure of the fitted copula to that of
% 	the original data. Specifically, the differences in Spearman’s rank
% 	correlation coefficient (Δρ_Spearmann) and Kendall’s tau (Δτ_Kendall)
% 	were computed between a Monte Carlo simulation of the fitted copula
% 	distribution and the empirical values derived from the original sample.
% 	This provides a measure of how well the fitted model captures the
% 	dependency structure of the data. For multivariate non-stationary
% 	copulas, this analysis was extended by computing the average
% 	differences over all the bivariate
% 	sub-distributions and the considered time windows.



% input:
%  copulaAnalysis                           - a variable of type structure provided as the output of tsCopulaExtremes or
%                                              tsCopulaCompoundGPDMontecarlo functions



% output:
%  gofStatistics:                           - A variable of type structure containing:
%                                               xxx                         --
%
%

% M.H.Bahmanpour, 2025

%REFERENCES

% [1] Bahmanpour, M.H., Mentaschi, L., Tilloy, A., Vousdoukas, M.,
%     Federico, I., Coppini, G., and Feyen, L., 2025,
%     Transformed-Stationary EVA 2.0: A Generalized Framework for
%     Non-stationary Joint Extreme Analysis (submitted to Hydrology and
%     Earth System Sciences; Feb 2025)
% [2] Mentaschi, L., Vousdoukas, M. I., Voukouvalas, E., Sartini, L.,
%     Feyen, L., Besio, G., & Alfieri, L. (2016). The
%     transformed-stationary approach: a generic and simplified methodology
%     for non-stationary extreme value analysis. Hydrology and Earth System
%     Sciences, 20(9), 3527–3547. https://doi.org/10.5194/hess-20-3527-2016
% [3] Genest, C., Rémillard, B., Beaudoin, D., Goodness-of-fit tests
%      for copulas: A review and a power study (Open Access),(2009) Insurance:
%      Mathematics and Economics, 44 (2), pp. 199-213, doi:
%      10.1016/j.insmatheco.2007.10.005
% [4] Hofert, M., Kojadinovic, I., Mächler, M. & Yan, J. Elements of
%      Copula Modeling with R (Springer, New York, 2018).

%%%%%%%%%%%%%%%%%%%%%%


% setting the default parameters

smoothInd = copulaAnalysis.copulaParam.smoothInd;
copulaFamily = copulaAnalysis.copulaParam.family;
copulaParam = copulaAnalysis.copulaParam;

%read non-stationary joint extremes
jointExtremes=copulaAnalysis.jointExtremes;

% calculate psuedo-observations (needed for Cramer-von Mises statistic)

if iscell(jointExtremes)
    uSample= cellfun(@(x) tsPseudoObservations(x),jointExtremes,'UniformOutput',0);
else
    uSample= cellfun(@(x) tsPseudoObservations(x),{jointExtremes},'UniformOutput',0);
end


    

    if strcmpi(copulaFamily,'Gaussian')

        %obtain an estimation of rho based on psudo-observations

        rhoC=  cellfun(@(x) copulafit('gaussian',x),uSample,'UniformOutput',0);
       

    elseif strcmpi(copulaFamily, 'Gumbel') || strcmpi(copulaFamily, 'Clayton') || strcmpi(copulaFamily, 'Frank')

        rhoC=  cellfun(@(x) copulafit(copulaFamily,x),uSample,'UniformOutput',0);
        
    end
    copulaParam.rho=rhoC;

    if strcmpi(copulaParam.family, 'gaussian')

        Y=cellfun(@(x,y) copulacdf('gaussian',x,y),uSample,copulaParam.rho,'UniformOutput',0);

    elseif strcmpi(copulaParam.family, 'clayton') || strcmpi(copulaParam.family, 'frank') || strcmpi(copulaParam.family, 'gumbel')

        Y=cellfun(@(x,y) copulacdf(copulaParam.family,x,y),uSample,copulaParam.rho,'UniformOutput',0);
    end

    snSample=cellfun(@(x,y) sum((tsEmpirical(x) - y) .^ 2),uSample,Y);

    gofStatistics.snSample=mean((snSample));
    gofStatistics.copulaFamily=copulaFamily;

    % calculate correlations in probability space
    jointExtremeMonovariateProbNS=copulaAnalysis.jointExtremeMonovariateProbNS;
   
    jointExtremesResampled = monteCarloAnalysis.resampleProb; %from Monte-Carlo simulations
    
    if iscell(jointExtremeMonovariateProbNS)
        corrKendallSample=cellfun(@(x) nonzeros(triu(corr(x,'type','Kendall'),1)),jointExtremeMonovariateProbNS,'UniformOutput',0);

        corrSpearmanSample=cellfun(@(x) nonzeros(triu(corr(x,'type','Spearman'),1)),jointExtremeMonovariateProbNS,'UniformOutput',0);
    else
        corrKendallSample=cellfun(@(x) nonzeros(triu(corr(x,'type','Kendall'),1)),{jointExtremeMonovariateProbNS},'UniformOutput',0);

        corrSpearmanSample=cellfun(@(x) nonzeros(triu(corr(x,'type','Spearman'),1)),{jointExtremeMonovariateProbNS},'UniformOutput',0);
    end

    corrKendallMonte=cellfun(@(x) nonzeros(triu(corr(x,'type','Kendall'),1)),jointExtremesResampled,'UniformOutput',0);
    corrSpearmanMonte=cellfun(@(x) nonzeros(triu(corr(x,'type','Spearman'),1)),jointExtremesResampled,'UniformOutput',0);

    if copulaParam.nSeries==2
        corrKendallSample=num2cell(smoothdata(cell2mat(corrKendallSample),'movmean',smoothInd));
        corrSpearmanSample=num2cell(smoothdata(cell2mat(corrSpearmanSample),'movmean',smoothInd));
        corrKendallMonte=num2cell(smoothdata(cell2mat(corrKendallMonte),'movmean',smoothInd));
        corrSpearmanMonte=num2cell(smoothdata(cell2mat(corrSpearmanMonte),'movmean',smoothInd));
    elseif copulaParam.nSeries==3
        for kk=1:4
            if kk==1
                XX=corrKendallSample;
            elseif kk==2
                XX=corrSpearmanSample;
            elseif kk==3
                XX=corrKendallMonte;
            elseif kk==4
                XX=corrSpearmanMonte;
            end
            N = length(corrKendallSample); % Number of 3x3 cell arrays

            % Preallocate arrays to store extracted values
            comp_12 = cell(1, N);
            comp_13 = cell(1, N);
            comp_23 = cell(1, N);

            % Extract the required components
            for ij = 1:N
                comp_12{ij} = XX{ij}(1); % Extract (1,2) component
                comp_13{ij} = XX{ij}(2); % Extract (1,3) component
                comp_23{ij} = XX{ij}(3); % Extract (2,3) component
            end
            comp_12=num2cell(smoothdata(cell2mat(comp_12),'movmean',smoothInd));
            comp_13=num2cell(smoothdata(cell2mat(comp_13),'movmean',smoothInd));
            comp_23=num2cell(smoothdata(cell2mat(comp_23),'movmean',smoothInd));
            for ij = 1:N
                XX{ij}(1)=comp_12{ij} ; % Extract (1,2) component
                XX{ij}(2)=comp_13{ij} ; % Extract (1,3) component
                XX{ij}(3)=comp_23{ij} ; % Extract (2,3) component
            end

            
            if kk==1
                corrKendallSample=XX;

            elseif kk==2
                corrSpearmanSample=XX;
            elseif kk==3
                corrKendallMonte=XX;
            elseif kk==4
                corrSpearmanMonte=XX;
            end
        end
    end
    kendallDelta=cellfun(@(x,y) abs(x-y),corrKendallSample,corrKendallMonte,'UniformOutput',0);
    spearmanDelta=cellfun(@(x,y) abs(x-y),corrSpearmanSample,corrSpearmanMonte,'UniformOutput',0);

    gofStatistics.corrKendallSampleDelta=mean(cellfun(@(x) (mean(x)),kendallDelta,'UniformOutput',1));
    gofStatistics.corrSpearmanSampleDelta=mean(cellfun(@(x) (mean(x)),spearmanDelta));

    gofStatistics.corrSpearmanSamplex=corrSpearmanSample;
    gofStatistics.corrSpearmanMontex=corrSpearmanMonte;
    gofStatistics.corrKendallSamplex=corrKendallSample;
    gofStatistics.corrKendallMontex=corrKendallMonte;

end





