function [samplingAnalysis] = tsCopulaSampleJointPeaksMultiVariatePruning(inputtimestamps,inputtimeseries,varargin)

%tsCopulaSampleJointPeaksMultiVariatePruning       multivariate peak-over-threshold sampling of compound events

% [samplingAnalysis] = tsCopulaSampleJointPeaksMultiVariatePruning(inputtimestamps,inputtimeseries,varargin)
%                     returns a variable of type structure containing various parameters 
%                     related to the sampled compound extremes

% This function samples the compound events (multivariate peaks) from
% inputtimeseries, based on a multivariate peak over threshold
% approach.inputtimeseries must be stationarized beforehand. The function
% keeps record of both compound peak events (defined as events during which
% all monovariates exceed their respective thresholds) and compound
% non-peak events (defined as events during which at least one (but not
% all) monovariate exceed its threshold value).

% POT sampling is conducted on each stationarized series x(t)
% (inputtimeseries need to be stationarized prior to calling this
% function), selecting multivariate peaks within a defined maximum time
% interval 〖Δt〗_multivariate. A challenge with this approach is the
% potential for multiple combinations of univariate peaks within the
% interval 〖Δt〗_multivariate. In tsEVA 2.0, this issue is addressed by
% prioritizing joint peaks with the largest mean values (average of
% univariate peak values), iteratively removing all other peak
% combinations.

% The compound (joint) extremes are sampled using the stationarized series.
% Transformation of the non-stationary series to stationarized series and
% the calculation of marginal distributions are performed using method of
% Mentaschi, et al., 2016 [1], applied on each margin separately.

% input:
%  inputtimestamps                           - 1d array with length nt, time stamps for the input
%                                              time series. must be the same for all the time series
%  inputtimeseries                           - 2d array with size [nt x n], where n is the number of variables.

% other (optional) inputs:

%  samplingThresholdPrct                     - 1d array with length n, percentile threshold to be used for each time
%                                              series. 
%  minPeakDistanceInDaysMonovarSampling      - 1d array with length n, minimum time distance (in days) among peaks of the same
%                                              variable used for sampling
%  maxPeakDistanceInDaysMultivarSampling     - maximum time distance (in days) among peaks of different variables
%                                              for the peaks to be considered joint. 1d array with  
%                                              length(maxPeakDistanceInDaysMultivarSampling) =size(nchoosek([1:n],2),1),
%                                              where nchoosek([1:n],2) shows the format in which maxPeakDistanceInDaysMultivarSampling
%                                              will be interpreted. Alternatively, it can be 1d array with length 1                                             
% marginalAnalysis                           - a 1d cell array with length n, where each cell contains result of applying tsEvaNonStationary
%                                              on each variable. This will be used in calculation of non-stationary thresholds.



% output:
%  samplingAnalysis:                           - A variable of type structure containing:
%                                               jointextremes                        -- 3d array with size [njointevents, n, 2] , where 
%                                                                                      njointevents is the number of found events, 
%                                                                                      and n is the number of variables. 
%                                                                                      jointextremes[ijointevent, ivariable, 1] contains 
%                                                                                      the time stamp of the peak related to the ijointevent-th
%                                                                                      event and ivariable-th variable, while 
%                                                                                      jointextremes[ijointevent, ivariable,2] contains the value
%                                                                                      of the peak.
%                                               jointextremes2                       -- 3d array with size [njointnonpeakevents, n, 2] , where 
%                                                                                      njointnonpeakevents is the number of found non-peak joint events, 
%                                                                                      and n is the number of variables. Stacking in the third dimension 
%                                                                                      is similar with that of jointextremes
%                                                                                      
%                                               thresholdsC                          -- 1d array of length n, containing thresholds in selecting peaks of
%                                                                                      monovariates
%                                               peaksjointidx                        -- indices of compound peak events
%                                               peaksjointidx2                       -- indices of compound non-peak events
                                             

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

args.samplingThresholdPrct=[99,99];                      
args.minPeakDistanceInDaysMonovarSampling=[3,3];         
args.maxPeakDistanceInDaysMultivarSampling=3;           
args.marginalAnalysis=cell(1,size(inputtimeseries,2));   % used for calculation of non-stationary parameters
args.samplingOrder=0;
args.peakType='allExceedThreshold';

% parsing of input parameters, overrides if different with the default
args = tsEasyParseNamedArgs(varargin, args);


samplingThresholdPrct=args.samplingThresholdPrct;
minPeakDistanceInDaysMonovarSampling=args.minPeakDistanceInDaysMonovarSampling;
maxPeakDistanceInDaysMultivarSampling=args.maxPeakDistanceInDaysMultivarSampling;
marginalAnalysis=args.marginalAnalysis;
peakType=args.peakType;
samplingOrder=args.samplingOrder;

% adjust maxPeakDistanceInDaysMultivarSampling to a compatible format 
if size(maxPeakDistanceInDaysMultivarSampling,2)~=size(nchoosek([1:size(inputtimeseries,2)],2),1)
    maxPeakDistanceInDaysMultivarSampling=...
        repmat(maxPeakDistanceInDaysMultivarSampling,1,size(nchoosek([1:size(inputtimeseries,2)],2),1));
end

numVar=size(inputtimeseries,2);         
dt = tsEvaGetTimeStep(inputtimestamps); 

inputtimeseriesCell=mat2cell(inputtimeseries,size(inputtimeseries,1),ones(1,size(inputtimeseries,2)));
thresholdsArray=cellfun(@(x,y) prctile(x,y),...
    inputtimeseriesCell,...
    transpose(num2cell(samplingThresholdPrct(1:numVar))) );

minPeakDistanceMonovarSampling = minPeakDistanceInDaysMonovarSampling/dt;

%find peaks of monovariate 
[pksCell,indxCell]=...
    cellfun(@(x,y) findpeaks(x,'minpeakdistance',y),inputtimeseriesCell,num2cell(minPeakDistanceMonovarSampling),'UniformOutput',0);
pksTimeCell=cellfun(@(x,y) x(y),repmat({inputtimestamps},1,numVar),indxCell,'UniformOutput',0);

%combine all monovariate peaks
combinedPeaksTime=cell2mat(pksTimeCell'); %
combinedPeaksIndex=cell2mat(indxCell');
combinedPeaks=cell2mat(pksCell');

% sort combined peaks in a time-ascending order
[combinedPeaksTime,sortIndexCombinedPeaks]=sort(combinedPeaksTime,'ascend'); 
combinedPeaks=combinedPeaks(sortIndexCombinedPeaks);
combinedPeaksIndex=combinedPeaksIndex(sortIndexCombinedPeaks);

% create an array indicating each peak's belonging to certain monovariate
varIdCell=num2cell(1:numVar);
idPeakCell=cellfun(@(x,y) y.*ones(size(x,1),1),pksTimeCell,varIdCell,'UniformOutput',0);
combinedPeaksId=cell2mat(idPeakCell');
combinedPeaksId=combinedPeaksId(sortIndexCombinedPeaks);

%pre-assign arrays to store information regarding peaks
indexJointPeaks=[]; 
indexJointNonPeaks=[]; 

% determines the order with which to proceed through existing combinations (e.g., first 1-2; next 2-3, then 1-3)
numVarId=unique(combinedPeaksId);
combinationVariates=nchoosek([1:size(numVarId,1)],2);
[~,sortIndexCombinations]=sort(diff(combinationVariates,[],2));
combinationVariates=combinationVariates(sortIndexCombinations,:);

for combinedPeaksCount=1:length(combinedPeaksTime)-1

    combinedPeaksIdPair=numVarId(numVarId~=combinedPeaksId(combinedPeaksCount));
    indicesPerEventCell={};
    for combinedPeaksIdPairCount=1:length(combinedPeaksIdPair)
        eventPeaksId(1,1,1:length([combinedPeaksId(combinedPeaksCount),combinedPeaksIdPair(combinedPeaksIdPairCount)]))=...
            [combinedPeaksId(combinedPeaksCount),combinedPeaksIdPair(combinedPeaksIdPairCount)];
        indexToCombinationVariates=find(all(~all(bsxfun(@minus,combinationVariates,eventPeaksId),2),3));
        passIndices=find(combinedPeaksTime(combinedPeaksCount+1:end)-combinedPeaksTime(combinedPeaksCount)<=...
            maxPeakDistanceInDaysMultivarSampling(indexToCombinationVariates)&combinedPeaksId(combinedPeaksCount+1:end)==...
            combinedPeaksIdPair(combinedPeaksIdPairCount));
        indicesPerEventCell=[indicesPerEventCell,{combinedPeaksCount+passIndices}];
    end
    if any(cellfun(@(x) isempty(x),indicesPerEventCell))
        continue
    end
    %generate indices of all possible combinations that could be formed by
    %peaks that are sampled
    indicesCombination=nchoosek([cell2mat(indicesPerEventCell');combinedPeaksCount],numVar); 

     %generate id list of compound peaks
    compoundPeaksId=combinedPeaksId(indicesCombination);   
    if size(compoundPeaksId,2)==1 %a column vector needs to change to a row vector for proper functioning of next lines
        compoundPeaksId=compoundPeaksId';
    end

    % ideally, each variate must only contribute one value to the compound
    % sample;
    desiredPeaksId(1,1,1:numVar) = 1:numVar;    
    indexToDesiredPeaksId = all(any(bsxfun(@eq,compoundPeaksId,desiredPeaksId),2),3); 
    indexToDesiredPeaksId = find(indexToDesiredPeaksId);
    indicesCombination=indicesCombination(indexToDesiredPeaksId,:);
    compoundPeaks=combinedPeaks(indicesCombination);
    compoundPeaksId=compoundPeaksId(indexToDesiredPeaksId,:);

    if size(compoundPeaks,2)==1 %a column vector needs to change to a row vector for proper functioning of next lines
        compoundPeaks=compoundPeaks';
    end

    % among the compound peaks, first we check if at least one combination
    % matches definition of a joint peak event; if two or more events match
    % the definition of a joint peak event, all are kept; if no combination
    % matches with definition of a joint peak event, we then search for
    % combinations that match with definition of a joint non-peak event; if
    % two or more events match with the definition of a joint non-peak
    % event, all are kept;
    thresholdsArrayReGrouped=thresholdsArray(compoundPeaksId);
    if any(all(compoundPeaks>=thresholdsArrayReGrouped,2)) 
        indextojointpeaks=find((all(compoundPeaks>=thresholdsArrayReGrouped,2)));
        if size(indextojointpeaks,1)==1  
            indexJointPeaks=[indexJointPeaks,indicesCombination(indextojointpeaks,:)];
        elseif size(indextojointpeaks,1)>1 
            indjointpeaks=indicesCombination(indextojointpeaks,:);
            indicesJointPeaksT=indjointpeaks';
            indicesJointPeaksT=indicesJointPeaksT(:);
            indexJointPeaks=[indexJointPeaks,indicesJointPeaksT'];
        end
    elseif any(any(compoundPeaks>=thresholdsArrayReGrouped,2))  
        indicesToJointNonPeaks=find((any(compoundPeaks>=thresholdsArrayReGrouped,2)));
        if size(indicesToJointNonPeaks,1)==1  
            indexJointNonPeaks=[indexJointNonPeaks,indicesCombination(indicesToJointNonPeaks,:)];
        elseif size(indicesToJointNonPeaks,1)>1 
            indicesJointNonPeaks=indicesCombination(indicesToJointNonPeaks,:);
            indicesJointNonPeaksT=indicesJointNonPeaks';
            indicesJointNonPeaksT=indicesJointNonPeaksT(:);
            indexJointNonPeaks=[indexJointNonPeaks,indicesJointNonPeaksT'];
        end
    end


end

%sampling of all joint peak and non-peak events from the total series based
%on the indices found from the above loop

jointPeaksId=combinedPeaksId(indexJointPeaks);
jointPeaksTime=combinedPeaksTime(indexJointPeaks);  
jointPeaks=combinedPeaks(indexJointPeaks);
jointNonPeaksId=combinedPeaksId(indexJointNonPeaks);
jointNonPeaksTime=combinedPeaksTime(indexJointNonPeaks);  
jointNonPeaks=combinedPeaks(indexJointNonPeaks);

jointPeaksTimeColumnWise=[];
jointPeaksColumnWise=[];
indexJointPeaksColumnWise=[];
jointNonPeaksTimeColumnWise=[];
jointNonPeaksColumnWise=[];
indexJointNonPeaksColumnWise=[];


for jj=1:numVar
    jointPeaksIdJJ=find(jointPeaksId==jj);
    jointPeaksTimeJJ=jointPeaksTime(jointPeaksIdJJ);
    jointPeaksJJ=jointPeaks(jointPeaksIdJJ);
    indexJointPeaksJJ=indexJointPeaks(jointPeaksIdJJ);
    jointPeaksTimeColumnWise=[jointPeaksTimeColumnWise,jointPeaksTimeJJ];%eventstime
    jointPeaksColumnWise=[jointPeaksColumnWise,jointPeaksJJ];%eventpeaks
    indexJointPeaksColumnWise=[indexJointPeaksColumnWise,indexJointPeaksJJ'];%eventid

    jointNonPeaksIdJJ=find(jointNonPeaksId==jj);
    jointNonPeaksTimeJJ=jointNonPeaksTime(jointNonPeaksIdJJ);
    jointNonPeaksJJ=jointNonPeaks(jointNonPeaksIdJJ);
    indexJointNonPeaksJJ=indexJointNonPeaks(jointNonPeaksIdJJ);
    jointNonPeaksTimeColumnWise=[jointNonPeaksTimeColumnWise,jointNonPeaksTimeJJ];%eventstime2
    jointNonPeaksColumnWise=[jointNonPeaksColumnWise,jointNonPeaksJJ];%eventpeaks2
    indexJointNonPeaksColumnWise=[indexJointNonPeaksColumnWise,indexJointNonPeaksJJ'];%eventid2
end


%pruning of the overlapping events

jointPeaksTimeTotal=[jointPeaksTimeColumnWise;jointNonPeaksTimeColumnWise];
jointPeaksTotal=[jointPeaksColumnWise;jointNonPeaksColumnWise];
jointIdTotal=[indexJointPeaksColumnWise;indexJointNonPeaksColumnWise];

%abide-by ordering of variables
if find(samplingOrder)
    nonrealisticIndices=find(jointPeaksTimeTotal(:,samplingOrder(2))-jointPeaksTimeTotal(:,samplingOrder(1))<0);
    jointPeaksTimeTotal(nonrealisticIndices,:)=[];
    jointPeaksTotal(nonrealisticIndices,:)=[];
    jointIdTotal(nonrealisticIndices,:)=[];
else
    nonrealisticIndices=[];
end
% peaks need to be sorted in a descending manner
[~,idSort]=sort(mean(jointPeaksTotal,2),'descend');
jointPeaksTimeTotal=jointPeaksTimeTotal(idSort,:);
jointPeaksTotal=jointPeaksTotal(idSort,:);
jointIdTotal=jointIdTotal(idSort,:);

idPeaksArtificial=[ones(size(jointPeaksTimeColumnWise,1),1);2*ones(size(jointNonPeaksTimeColumnWise,1),1)];
idPeaksArtificial(nonrealisticIndices,:)=[];
idPeaksArtificial=idPeaksArtificial(idSort);

minArrayPeaksTime=min(jointPeaksTimeTotal,[],2);
maxArrayPeaksTime=max(jointPeaksTimeTotal,[],2);
indicesToRemove=[];

for jx=1:size(jointPeaksTimeTotal,1) %loop through all peaks starting from
                                     % the largest ones


    if(any(indicesToRemove==jx))     %because that particular joint event was
                                     %ruled out by a prior check
        continue
    end
    eventTime=jointPeaksTimeTotal(jx,:);

    indicesNonOverlap=(min(eventTime)>maxArrayPeaksTime|max(eventTime)<minArrayPeaksTime);
    indicesNonOverlap(1:jx)=1;       % this line ensures that all the peaks that
                                     % were previously checked, are no longer accounted by
    if find(jointPeaksTimeTotal(~indicesNonOverlap))
        indicesToRemove=[indicesToRemove,find(~indicesNonOverlap)'];
    end

end

if isempty(indicesToRemove)
    disp('No case for pruning found')
else
    disp([num2str(round(100*(size(unique(indicesToRemove),2)/size(jointPeaksTimeTotal,1))*10)/10),' %',' of sampled events pruned due to overlapping of events'])
end

jointPeaksTimeTotal(indicesToRemove,:)=[];
jointPeaksTotal(indicesToRemove,:)=[];
jointIdTotal(indicesToRemove,:)=[];
idPeaksArtificial(indicesToRemove)=[];

jointPeaksTimeColumnWise=jointPeaksTimeTotal(idPeaksArtificial==1,:);
jointNonPeaksTimeColumnWise=jointPeaksTimeTotal(idPeaksArtificial==2,:);
jointPeaksColumnWise=jointPeaksTotal(idPeaksArtificial==1,:);
jointNonPeaksColumnWise=jointPeaksTotal(idPeaksArtificial==2,:);
indexJointPeaksColumnWise=jointIdTotal(idPeaksArtificial==1,:);
indexJointNonPeaksColumnWise=jointIdTotal(idPeaksArtificial==2,:);

jointExtremeIndices=combinedPeaksIndex(indexJointPeaksColumnWise);
if size(combinedPeaksIndex(indexJointNonPeaksColumnWise),2)~=1
    peakIndicesAll=[combinedPeaksIndex(indexJointPeaksColumnWise);combinedPeaksIndex(indexJointNonPeaksColumnWise)];
else
    peakIndicesAll=[combinedPeaksIndex(indexJointPeaksColumnWise);combinedPeaksIndex(indexJointNonPeaksColumnWise)'];
end
disp([num2str(((size((jointPeaksTimeColumnWise),1)))),' ','Compound peak events found'])
% disp([num2str(((size((jointNonPeaksTimeColumnWise),1)))),' ','Compound non-peak events found'])
disp(['average time among peaks of ',...
    num2str(mean(abs(diff(jointPeaksTimeColumnWise,[],2)))),' days',...
    ' with minimum of ',num2str(min(abs(diff(jointPeaksTimeColumnWise,[],2)))),' days',...
    ' and maximum of ', num2str(max(abs(diff(jointPeaksTimeColumnWise,[],2)))),' days'])
jointextremes=cat(3,jointPeaksTimeColumnWise,jointPeaksColumnWise);
jointextremes2=cat(3,jointNonPeaksTimeColumnWise,jointNonPeaksColumnWise);
if ~isempty(marginalAnalysis)

    nonStatSeriesC=cellfun(@(x) x{2}.nonStatSeries,marginalAnalysis,'UniformOutput',0);
    nonStatSeries=[nonStatSeriesC{:}];
    jointExtremesNS=[];

    for ii=1:size(jointextremes,2)
       % jointExtremesNS=[jointExtremesNS,nonStatSeries(jointExtremeIndices(:,ii),ii)];
       if strcmpi(peakType,'allExceedThreshold')
           jointExtremesNS=[jointExtremesNS,nonStatSeries(jointExtremeIndices(:,ii),ii)];
       elseif strcmpi(peakType,'anyExceedThreshold')
           jointExtremesNS=[jointExtremesNS,nonStatSeries(peakIndicesAll(:,ii),ii)];
       end

    end
else
    jointExtremesNS=[];
end
if ~isempty(marginalAnalysis)
    trendSeries=cellfun(@(x) x{2}.trendSeries,marginalAnalysis,'UniformOutput',0);
    stdDevSeries=cellfun(@(x) x{2}.stdDevSeries,marginalAnalysis,'UniformOutput',0);

    thresholdsNonStation=cellfun(@(x,y,z) y*z+x,trendSeries,stdDevSeries,num2cell(thresholdsArray),'UniformOutput',0);%thresholdPotNS
    samplingAnalysis.thresholdsNonStation=thresholdsNonStation;
else
    samplingAnalysis.thresholdsNonStation=[];

end

samplingAnalysis.jointextremes=jointextremes;
samplingAnalysis.jointextremes2=jointextremes2;
samplingAnalysis.thresholdsC=thresholdsArray;
samplingAnalysis.jointExtremeIndices=jointExtremeIndices;
samplingAnalysis.peakIndicesAll=peakIndicesAll;
samplingAnalysis.jointExtremesNS=jointExtremesNS;


