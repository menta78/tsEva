
function [monteCarloAnalysis] = tsCopulaCompoundGPDMontecarlo(copulaAnalysis, varargin)
%tsCopulaCompoundGPDMontecarlo pefrom Monte-Carlo simulation (resampling) from a
%pre-determined copula function

% [copulaAnalysis] = tsCopulaCompoundGPDMontecarlo(copulaAnalysis,varargin)
%                    returns results of Monte-Carlo simulation including
%                    resampled data in probability and data space
%                    (i.e., monteCarloRsmpl, and resampleProb)



% input:

%  copulaAnalysis                              - a variable of type structure containing various parameters
%                                                of the fitted copula. Need to be the output of tsCopulaCompoundGPD function

% other (optional) inputs:

%  nResample                                   - 1d scalar indicating size of the Monte-Carlo
%                                                simulation to be performed. The default value is 1000
%  timeIndex                                   - a scalar parameter for indexing non-stationary parameters


% output:

%  CopulaAnalysis:                           - A variable of type structure same as the input with two additional appended variables:

%                                               monteCarloRsmpl        -- Resampled return levels. in case of a time-varying copula, a 1d cell array of length
%                                                                       matching with the number of time windows adopted for copula calculation;
%                                                                       in case of a stationary copula, a 2d array of size [nResample x nVar], where 
%                                                                       nVar indicates number of variables (e.g, 3 for trivariate case)
%                                               resampleProb         -- Resampled return probabilities. in case of a time-varying copula, a 1d cell array 
%                                                                       of length matching with the number of time windows adopted for copula calculation;
%                                                                       in case of a stationary copula, a 2d array of size [nResample x nVar], where 
%                                                                       nVar indicates number of variables (e.g, 3 for trivariate case)
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

%%%%%%%%%%%%%%%%%%%%%%

% setting the default parameters

args.timeIndex = 'middle'; 
args.nResample=1000;
args.nonStationarity = 'marginsandcoupling'; %two switches: "margins" ; "marginsandcoupling" 

args = tsEasyParseNamedArgs(varargin, args);

nonStationarity=args.nonStationarity;
timeIndex = args.timeIndex;
nResample=args.nResample;

%read input data
methodology=copulaAnalysis.methodology;
copulaParam=copulaAnalysis.copulaParam;
nSeries = copulaParam.nSeries;
copulaFamily = copulaParam.family;
marginalAnalysis=copulaAnalysis.marginalAnalysis;
timeVaryingCopula=copulaAnalysis.timeVaryingCopula;

%differentiate between the way time-varying and time-invariant copula need
%to be dealt with
%resampling from the copula function using the copularnd function
switch timeVaryingCopula
    case false
        
        
            if strcmpi(copulaFamily, 'Gaussian')

                resampleProb = tsCopulaRnd('gaussian', copulaParam.rho, nResample, []);

            elseif strcmpi(copulaFamily, 'Gumbel') || strcmpi(copulaFamily, 'Clayton') || strcmpi(copulaFamily, 'Frank')

                resampleProb = tsCopulaRnd(copulaFamily, copulaParam.rho, nResample,...
                    copulaAnalysis.jointExtremeMonovariateProbNS);

            else
                error(['copulaFamily not supported: ' copulaFamily]);
            end
        
    case true
        %for the case of a time-varying copula
        
        rhoCell=copulaParam.rho;
        nResampleC=repmat({nResample},1,size(rhoCell,2));
        resampleProb=cell(size(rhoCell));
            if strcmpi(nonStationarity,'margins') & strcmpi(copulaFamily, 'Gaussian')
                rhoCell=repmat({mean(cellfun(@(x) x(triu(true(size(x)),1)),rhoCell,'UniformOutput',1))},1,size(resampleProb,2));
            elseif strcmpi(nonStationarity,'margins') 
                rhoCell=repmat({mean([rhoCell{:}])},1,size(resampleProb,2));
            end
                   
                    resampleProb=cellfun(@(x,y,uSmpl) tsCopulaRnd(copulaFamily,x,y,uSmpl),...
                        rhoCell, nResampleC, copulaAnalysis.jointExtremeMonovariateProbNS, ...
                        'UniformOutput', 0);     
              
                if strcmpi(nonStationarity,'margins') 
                    monteCarloAnalysis.copulaParam=copulaParam;
                    monteCarloAnalysis.copulaParam.rhoMean=rhoCell;
                end


        
end


% on the basis of the timeIndex, find the non-stationary values of the
% thresold and scale parameter

switch timeVaryingCopula
    case false
        monteCarloRsmpl=cell(1,length(copulaFamily));
       
            for ivar = 1:nSeries
                %if no timeindex is set by the user use time-index to assess
                % non-stationarity parameters at half the length of the time series
                nonStatEvaParams = marginalAnalysis{ivar}{1};
                % statTransData = marginalAnalysis{ivar}{2};

                if strcmpi(timeIndex,'first') & ivar==1
                    timeIndex=1;
                    fprintf(['conversion of Monte-Carlo probabilities to data space is \n',...
                        'based on non-stationary values evaluated at\n' ...
                        'the first timeindex'])
                elseif strcmpi(timeIndex,'last') & ivar==1
                    fprintf(['conversion of Monte-Carlo probabilities to data space is \n',...
                        'based on non-stationary values evaluated at\n' ...
                        'the last timeindex'])
                    timeIndex=(length(nonStatEvaParams(2).parameters.threshold));
                elseif strcmpi(timeIndex,'middle') & ivar==1
                    fprintf(['conversion of Monte-Carlo probabilities to data space is \n',...
                        'based on non-stationary values evaluated at\n' ...
                        'the middle timeindex'])
                    timeIndex=ceil(length(nonStatEvaParams(2).parameters.threshold)/2);
                elseif isnumeric(timeIndex) & ivar==1
                    if timeIndex<1 || timeIndex>length(nonStatEvaParams(2).parameters.threshold)
                        error('timeIndex parameter must be chosen from {"first","last","middle"} or a valid index')
                    end
                end
                timeIndexArray=timeIndex;
                % transfrom probabilities to data scale using inverse sampling law
                % no scaling is needed since thrshld parameter already transforms data with
                % lowest probability corresponding with thrshld value
                monteCarloRsmpl{1}(:,ivar) = computeResampledLevels(resampleProb{1}(:,ivar), nonStatEvaParams, timeIndex,methodology);

            end
        
    case true

        %in case of a time-varying copula
        monteCarloRsmplCell=cell(size(resampleProb));
        timeStampsByTimeWindow=copulaParam.timeStampsByTimeWindow;

        timeStamps = marginalAnalysis{1}{2}.timeStamps;
        timeStampsCell=repmat({timeStamps},1,size(rhoCell,2));
        iixCell=cellfun(@(x,y) find(x>=min(y)&x<=max(y)),timeStampsCell,timeStampsByTimeWindow,'UniformOutput',0);
        timeIndexArray=cellfun(@(x) x(round(length(x)/2)),iixCell);
        if ~any(strcmpi(varargin,'timeindex'))
            disp('no timeindex set - middle timeindex (for each time-window) selected automatically')

        elseif isnumeric(timeIndex)

            fprintf(['numeric timeindex not accepted in case of \n',...
                'a time-varying copula use "first", "last" or "middle" instead\n' ...
                'middle timeindex for each time-window selected automatically'])
        elseif any(strcmpi(varargin,'first'))
            fprintf(['conversion of Monte-Carlo probabilities to data space is \n',...
                'based on non-stationary values evaluated at\n' ...
                'the first timeindex for each time-window'])
            timeIndexArray=cellfun(@(x) x(1),iixCell);
        elseif any(strcmpi(varargin,'last'))
            fprintf(['conversion of Monte-Carlo probabilities to data space is \n',...
                'based on non-stationary values evaluated at\n' ...
                'the last timeindex for each time-window'])
            timeIndexArray=cellfun(@(x) x(end),iixCell);
        elseif any(strcmpi(varargin,'middle'))
            fprintf(['conversion of Monte-Carlo probabilities to data space is \n',...
                'based on non-stationary values evaluated at\n' ...
                'the middle timeindex for each time-window'])

        end
           for ik=1:size(rhoCell,1)
               for ij=1:size(rhoCell,2)
                   
                   resampleProbTemp=resampleProb{ik,ij};

                   for ivar = 1:nSeries
                       nonStatEvaParams = marginalAnalysis{ivar}{1};
                       monteCarloRsmplCell{ik,ij}(:,ivar) = computeResampledLevels(resampleProbTemp(:,ivar), nonStatEvaParams, timeIndexArray(ij),methodology);

                   end
                  
               end
           end
        monteCarloRsmpl=monteCarloRsmplCell;

end
%append monteCarloRsmpl and resampleProb to the copulaAnalysis file of type
%structure
monteCarloAnalysis.monteCarloRsmpl=monteCarloRsmpl;

monteCarloAnalysis.resampleProb=resampleProb;

monteCarloAnalysis.timeIndexArray=timeIndexArray;

end


function monteCarloRsmpls = computeResampledLevels(resampleProb, nonStatEvaParams, timeIndex,methodology)

if strcmpi(methodology,'gpd')
    thrshld = nonStatEvaParams(2).parameters.threshold(timeIndex);
    scaleParam = nonStatEvaParams(2).parameters.sigma(timeIndex);
    shapeParam = nonStatEvaParams(2).parameters.epsilon;
    monteCarloRsmpls = gpinv(resampleProb, shapeParam, scaleParam, thrshld);
elseif strcmpi(methodology,'gev')
   mu = nonStatEvaParams(1).parameters.mu(timeIndex);
    scaleParam = nonStatEvaParams(1).parameters.sigma(timeIndex);
    shapeParam = nonStatEvaParams(1).parameters.epsilon;
    monteCarloRsmpls = gevinv(resampleProb, shapeParam, scaleParam, mu);
end

end

