function [resampleLevel, resampleProb] = tsCopulaCompoundGPDMontecarlo(copulaAnalysis,...
    nResample, ...
    varargin)


args.timeIndex = 1; %default time index to assess non-stationary return levels

args = tsEasyParseNamedArgs(varargin, args);

timeIndex = args.timeIndex;
copulaParam=copulaAnalysis.copulaParam;
marginalAnalysis=copulaAnalysis.marginalAnalysis;
copulaFamily = copulaParam.family;
timeVaryingCopula=copulaAnalysis.timeVaryingCopula;
switch timeVaryingCopula
    case false
        if strcmpi(copulaFamily, 'Gaussian')

            resampleProb = copularnd(copulaFamily, copulaParam.rho, nResample);
        elseif strcmpi(copulaFamily, 't')
            resampleProb = copularnd(copulaFamily, copulaParam.rho, copulaParam.nu, nResample);
        elseif strcmpi(copulaFamily, 'Gumbel') || strcmpi(copulaFamily, 'Clayton') || strcmpi(copulaFamily, 'Frank')
            resampleProb = copularnd(copulaFamily, copulaParam.rho, nResample);
        else
            error(['copulaFamily not supported: ' copulaFamily]);
        end

    case true
        resampleProb={};
        if strcmpi(copulaFamily, 'Gaussian') || strcmpi(copulaFamily, 'Gumbel') || strcmpi(copulaFamily, 'Clayton') || strcmpi(copulaFamily, 'Frank')
            rhoCell=copulaParam.rho;
            for ij=1:size(rhoCell,2)
                resampleProbT = copularnd(copulaFamily, rhoCell{ij}, nResample);
                resampleProb=[resampleProb,resampleProbT];
            end
        elseif strcmpi(copulaFamily, 't')
             rhoCell=copulaParam.rho;
             nuCell=copulaParam.nu;
             for ij=1:size(rhoCell,2)
                resampleProbT = copularnd(copulaFamily, rhoCell{ij}, nuCell{ij},nResample);
                resampleProb=[resampleProb,resampleProbT];
             end
            
        else
            error(['copulaFamily not supported: ' copulaFamily]);
        end
end


% on the basis of the timeIndex, find the non-stationary values of the
% thresold and scale parameter
nSeries = length(marginalAnalysis);
switch timeVaryingCopula
    case false

        for ivar = 1:nSeries

            nonStatEvaParams = marginalAnalysis{ivar}{1};
            thrshld = nonStatEvaParams(2).parameters.threshold(timeIndex);
            scaleParam = nonStatEvaParams(2).parameters.sigma(timeIndex);
            shapeParam = nonStatEvaParams(2).parameters.epsilon;
          
         
% transfrom probabilities to data scale using inverse sampling law
% no scaling is needed since thrshld parameter already transforms data with
% lowest probability corresponding with thrshld value
       resampleLevel(:,ivar)=gpinv(resampleProb(:,ivar),shapeParam, scaleParam, thrshld);
       
        end
    case true
        resampleLevelCell={};
        
        for ij=1:size(rhoCell,2)
            resampleLevel=[];
            
            resampleProbTemp=resampleProb{ij};
            for ivar = 1:nSeries

                nonStatEvaParams = marginalAnalysis{ivar}{1};
                thrshld = nonStatEvaParams(2).parameters.threshold(timeIndex);
                scaleParam = nonStatEvaParams(2).parameters.sigma(timeIndex);
                shapeParam = nonStatEvaParams(2).parameters.epsilon;

                resampleLevel(:,ivar)=gpinv(resampleProbTemp(:,ivar),shapeParam, scaleParam, thrshld);
    
            end
            resampleLevelCell=[resampleLevelCell,resampleLevel];
       
        end
        resampleLevel=resampleLevelCell;
       
end

