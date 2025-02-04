
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

args.timeIndex = 'middle'; %to activate reading of "timeindex" parameter if set by the user
args.nResample=1000;
args.nonStationarity = 'marginsandcoupling'; %margins ; coupling ; marginsandcoupling  
args = tsEasyParseNamedArgs(varargin, args);
nonStationarity=args.nonStationarity;
timeIndex = args.timeIndex;
nResample=args.nResample;
methodology=copulaAnalysis.methodology;
%obtain some details from the copulaAnalysis input variable
copulaParam=copulaAnalysis.copulaParam;
marginalAnalysis=copulaAnalysis.marginalAnalysis;
copulaFamily = copulaParam.family;
timeVaryingCopula=copulaAnalysis.timeVaryingCopula;

rng default
%differentiate between the way time-varying and time-invariant copula need
%to be dealt with
%resampling from the copula function using the copularnd function
switch timeVaryingCopula
    case false
        resampleProb=cell(1,length(copulaFamily));
        for iFamily=1:length(copulaFamily)
            if strcmpi(copulaFamily{iFamily}, 'Gaussian')

                resampleProb{iFamily} = copularnd('gaussian', copulaParam.rho{iFamily}, nResample);
            elseif strcmpi(copulaFamily{iFamily}, 't')
                resampleProb{iFamily} = copularnd('t', copulaParam.rho{iFamily}, copulaParam.nu{iFamily}, nResample);
            elseif strcmpi(copulaFamily{iFamily}, 'Gumbel') || strcmpi(copulaFamily{iFamily}, 'Clayton') || strcmpi(copulaFamily{iFamily}, 'Frank')
                
                    resampleProb{iFamily} = copularnd(copulaFamily{iFamily}, copulaParam.rho{iFamily}, nResample);
               
            else
                error(['copulaFamily not supported: ' copulaFamily]);
            end
        end
    case true
        %for the case of a time-varying copula
        
        
        rhoCell=copulaParam.rho;
        nuCell=copulaParam.nu;
        resampleProb=cell(size(rhoCell));
        for iFamily=1:length(copulaFamily)
            if strcmpi(copulaFamily{iFamily}, 'Gaussian') || strcmpi(copulaFamily{iFamily}, 'Gumbel') || strcmpi(copulaFamily{iFamily}, 'Clayton') || strcmpi(copulaFamily{iFamily}, 'Frank')
               
                for ij=1:size(rhoCell,2)
                    if strcmpi(nonStationarity,'marginsandcoupling')
                        resampleProb{iFamily,ij} = copularnd(copulaFamily{iFamily}, rhoCell{iFamily,ij}, nResample);
                    elseif strcmpi(nonStationarity,'margins')
                        rhoCell(iFamily,:)={mean([rhoCell{iFamily,:}])};
                        resampleProb{iFamily,ij} = copularnd(copulaFamily{iFamily}, rhoCell{iFamily,ij}, nResample);
                    end
                end
                if strcmpi(nonStationarity,'margins')
                    copulaAnalysis.copulaParam.rhoMean=rhoCell;
                end
            elseif strcmpi(copulaFamily{iFamily}, 't')
                
                for ij=1:size(rhoCell,2)
                    resampleProb{iFamily,ij} = copularnd(copulaFamily{iFamily}, rhoCell{iFamily,ij}, nuCell{iFamily,ij},nResample);
                  
                end

            else
                error(['copulaFamily not supported: ' copulaFamily{iFamily}]);
            end

        end
end


% on the basis of the timeIndex, find the non-stationary values of the
% thresold and scale parameter
nSeries = length(marginalAnalysis);
switch timeVaryingCopula
    case false
        resampleLevel=cell(1,length(copulaFamily));
        for iFamily=1:length(copulaFamily)
            for ivar = 1:nSeries
                %if no timeindex is set by the user use time-index to assess
                % non-stationarity parameters at half the length of the time series
                nonStatEvaParams = marginalAnalysis{ivar}{1};
                % statTransData = marginalAnalysis{ivar}{2};

                if strcmpi(timeIndex,'first') & ivar==1
                    timeIndex=1;
                elseif strcmpi(timeIndex,'last') & ivar==1
                    timeIndex=(length(nonStatEvaParams(2).parameters.threshold));
                elseif strcmpi(timeIndex,'middle') & ivar==1
                    timeIndex=ceil(length(nonStatEvaParams(2).parameters.threshold)/2);
                elseif isnumeric(timeIndex) & ivar==1
                    if timeIndex<1 || timeIndex>length(nonStatEvaParams(2).parameters.threshold)
                        error('timeIndex parameter must be chosen from {"first","last","middle"} or a valid index')
                    end
                end

                % transfrom probabilities to data scale using inverse sampling law
                % no scaling is needed since thrshld parameter already transforms data with
                % lowest probability corresponding with thrshld value
                resampleLevel{iFamily}(:,ivar) = computeResampledLevels(resampleProb{iFamily}(:,ivar), nonStatEvaParams, timeIndex,methodology);

            end
        end
    case true
        %in case of a time-varying copula
        resampleLevelCell=cell(size(resampleProb));
        inputtimestampsWindowCell=copulaParam.inputtimestampsWindowCell;
        

           timeStamps = marginalAnalysis{1}{2}.timeStamps;
           timeStampsCell=repmat({timeStamps},1,size(rhoCell,2));
           iixCell=cellfun(@(x,y) find(x>=min(y)&x<=max(y)),timeStampsCell,inputtimestampsWindowCell,'UniformOutput',0);
            timeIndexArray=cellfun(@(x) x(round(length(x)/2)),iixCell);
           if ~any(strcmpi(varargin,'timeindex'))
               disp('no timeindex set - middle timeindex (for each time-window) selected automatically')
                          
           elseif isnumeric(timeIndex)
            
               sprintf(['numeric timeindex not accepted in case of \n',...
                   'a time-varying copula use first last or middle instead\n' ...
                   'middle timeindex for each time-window selected automatically'])
           elseif any(strcmpi(varargin,'first'))
                disp('first timeindex (for each time-window) selected')
              timeIndexArray=cellfun(@(x) x(1),iixCell);
           elseif any(strcmpi(varargin,'last'))
                disp('last timeindex (for each time-window) selected')
                timeIndexArray=cellfun(@(x) x(end),iixCell);
           elseif any(strcmpi(varargin,'middle'))
                disp('middle timeindex (for each time-window) selected')
           end
           for ik=1:size(rhoCell,1)
               for ij=1:size(rhoCell,2)
                   resampleLevel=[];

                   resampleProbTemp=resampleProb{ik,ij};

                   for ivar = 1:nSeries
                       nonStatEvaParams = marginalAnalysis{ivar}{1};
                       resampleLevelCell{ik,ij}(:,ivar) = computeResampledLevels(resampleProbTemp(:,ivar), nonStatEvaParams, timeIndexArray(ij),methodology);

                   end
                  

               end
           end
        resampleLevel=resampleLevelCell;

end
%append resampleLevel and resampleProb to the copulaAnalysis file of type
%structure
copulaAnalysis.resampleLevel=resampleLevel;

copulaAnalysis.resampleProb=resampleProb;
if exist("timeIndexArray")

copulaAnalysis.timeIndexArray=timeIndexArray;
end
end


function resampleLevels = computeResampledLevels(resampleProb, nonStatEvaParams, timeIndex,methodology)

if strcmpi(methodology,'gpd')
    thrshld = nonStatEvaParams(2).parameters.threshold(timeIndex);
    scaleParam = nonStatEvaParams(2).parameters.sigma(timeIndex);
    shapeParam = nonStatEvaParams(2).parameters.epsilon;
    resampleLevels = gpinv(resampleProb, shapeParam, scaleParam, thrshld);
elseif strcmpi(methodology,'gev')
   mu = nonStatEvaParams(1).parameters.mu(timeIndex);
    scaleParam = nonStatEvaParams(1).parameters.sigma(timeIndex);
    shapeParam = nonStatEvaParams(1).parameters.epsilon;
    resampleLevels = gevinv(resampleProb, shapeParam, scaleParam, mu);
end

end

