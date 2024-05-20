
function [copulaAnalysis] = tsCopulaCompoundGPDMontecarlo(copulaAnalysis,varargin)
%tsCopulaCompoundGPDMontecarlo Monte-carlo simulation (resampling) from a
%pre-determined copula function
% copulaAnalysis] = tsCopulaCompoundGPDMontecarlo(copulaAnalysis,nResample,varargin)
%                    appends to the variable (of type structure)
%                    copulaAnalysis results of monte-carlo simulation
%                    (i.e., resampleLevel, and resampleProb)



% input:

%  copulaAnalysis                              - a variable of type structure containing various parameters
%                                                of the fitted copula. Need to be the output of tsCopulaCompoundGPD function

% other (optional) inputs:

%  nResample                                   - 1d scalar indicating size of the Monte-Carlo
%                                                simulation to be performed. The default value is 1000
%  timeIndex                                   - a scalar parameter for indexing non-stationary parameters


% output:

%  CopulaAnalysis:                           - A variable of type structure same as the input with two additional appended variables:

%                                               resampleLevel        -- Resampled return levels. in case of a time-varying copula, a 1d cell array of length
%                                                                       matching with the number of time windows adopted for copula calculation;
%                                                                       in case of a stationary copula, a 2d array of size [nResample x nVar], where 
%                                                                       nVar indicates number of variables (e.g, 3 for trivariate case)
%                                               resampleProb         -- Resampled return probabilities. in case of a time-varying copula, a 1d cell array 
%                                                                       of length matching with the number of time windows adopted for copula calculation;
%                                                                       in case of a stationary copula, a 2d array of size [nResample x nVar], where 
%                                                                       nVar indicates number of variables (e.g, 3 for trivariate case)
%

% M.H.Bahmanpour, 2024


%%%%%%%%%%%%%%%%%%%%%%

% setting the default parameters

args.timeIndex = 1; %to activate reading of "timeindex" parameter if set by the user
args.nResample=1000;
args = tsEasyParseNamedArgs(varargin, args);

timeIndex = args.timeIndex;
nResample=args.nResample;

%obtain some details from the copulaAnalysis input variable
copulaParam=copulaAnalysis.copulaParam;
marginalAnalysis=copulaAnalysis.marginalAnalysis;
copulaFamily = copulaParam.family;
timeVaryingCopula=copulaAnalysis.timeVaryingCopula;


%differentiate between the way time-varying and time-invariant copula need
%to be dealt with
%resampling from the copula function using the copularnd function
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
        %for the case of a time-varying copula
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
            %if no timeindex is set by the user use time-index to assess
            % non-stationarity parameters at half the length of the time series
            nonStatEvaParams = marginalAnalysis{ivar}{1};

            if ~any(strcmpi(varargin,'timeindex'))
                timeIndex=round(length(nonStatEvaParams(2).parameters.threshold)/2);
            end

            thrshld = nonStatEvaParams(2).parameters.threshold(timeIndex);
            % thrshld = thresholdValues(ivar);
            scaleParam = nonStatEvaParams(2).parameters.sigma(timeIndex);
            shapeParam = nonStatEvaParams(2).parameters.epsilon;


            % transfrom probabilities to data scale using inverse sampling law
            % no scaling is needed since thrshld parameter already transforms data with
            % lowest probability corresponding with thrshld value
            resampleLevel(:,ivar)=gpinv(resampleProb(:,ivar),shapeParam, scaleParam, thrshld);

        end
    case true
        %in case of a time-varying copula
        resampleLevelCell={};
        jointExtremeTimeStampsCell=copulaAnalysis.jointExtremeTimeStamps;
        jointExtremeTimeStampsCellMin=cellfun(@(x) (min(x)),jointExtremeTimeStampsCell,'UniformOutput',0);
        jointExtremeTimeStampsCellMax=cellfun(@(x) (max(x)),jointExtremeTimeStampsCell,'UniformOutput',0);

        for ij=1:size(rhoCell,2)
            resampleLevel=[];

            resampleProbTemp=resampleProb{ij};
            minTime=jointExtremeTimeStampsCellMin{ij};
            maxTime=jointExtremeTimeStampsCellMax{ij};
            for ivar = 1:nSeries
                nonStatEvaParams = marginalAnalysis{ivar}{1};

                if ~any(strcmpi(varargin,'timeindex'))
                    minTimex=minTime(ivar);
                    maxTimex=maxTime(ivar);
                    timeStamps = marginalAnalysis{ivar}{2}.timeStamps;
                    iix=find(timeStamps>=minTimex&timeStamps<=maxTimex);
                    timeIndex=iix(round(length(iix)/2));
                end
                thrshld = nonStatEvaParams(2).parameters.threshold(timeIndex);
                scaleParam = nonStatEvaParams(2).parameters.sigma(timeIndex);
                shapeParam = nonStatEvaParams(2).parameters.epsilon;

                resampleLevel(:,ivar)=gpinv(resampleProbTemp(:,ivar),shapeParam, scaleParam, thrshld);

            end
            resampleLevelCell=[resampleLevelCell,resampleLevel];

        end
        resampleLevel=resampleLevelCell;

end
%append resampleLevel and resampleProb to the copulaAnalysis file of type
%structure
copulaAnalysis.resampleLevel=resampleLevel;

copulaAnalysis.resampleProb=resampleProb;
