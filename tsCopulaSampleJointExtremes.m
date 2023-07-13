function [jointExtremes, thresholds] = tsCopulaSampleJointExtremes(inputTimeStamps, ...
    inputTimeSeries, ...
    thresholdPercentiles, ...
    minPeakDistanceInDays, ...
    maxDistanceMultivariatePeaksInDays)

% this function samples the joint extremes from a multivariate
% inputTimeSeries, based on a multivariate POT.
% This function assumes the input time series are stationary, so you need
% to feed it with the transformed time series.

% input data:
% - inputTimeStamps: 1d array with length nt, time stamps for the input
%   time series. Must be the same for all the time series
% - inputTimeSeries: data of each time series. 2d array with size [nt x N], where N is the number of
%   variables we need to analyze.
% - thresholdPercentiles: percentile threshold to be used for each time series.
%   1d array with length N, where N is the number of variables to be
%   considered.
% - minPeakDistanceInDays: minimum time distance among peaks of the same
%   variable. 1d array with length N, where N is the number of variables to
%   be considered.
% - maxDistanceMultivariatePeaksInDays: maximum time distance among peaks
%   of different variables for the peaks to be considered joint.

% output data:
% - jointExtremes: array containing the joint extremes. 3d array with shape
%   [nJointEvents, N, 2], where nJointEvents is the number of found events,
%   N is the number of variables. jointExtremes[iJointEvent, iVariable, 1]
%   contains the time stamp of the peak related to the iJointEvent-th event
%   and iVariable-th variable, while jointExtremes[iJointEvent, iVariable, 2]
%   contains the value of the peak.
% - thresholds: 1d array with size N, containing the value of the threshold
%   for each variable.

thrsh1=prctile(inputTimeSeries(:,1),thresholdPercentiles(1));
thrsh2=prctile(inputTimeSeries(:,2),thresholdPercentiles(2));
thresholds=[thrsh1,thrsh2];


thrsdt = prctile(inputTimeSeries(:,1),thresholdPercentiles(1));
dt = tsEvaGetTimeStep(inputTimeStamps);
minPeakDistance = minPeakDistanceInDays(1)/dt;
[pks1,indx1] = findpeaks(inputTimeSeries(:,1),'MinPeakDistance',minPeakDistance,'MinPeakHeight',thrsdt);
TimeStamps1=inputTimeStamps(indx1);


thrsdt = prctile(inputTimeSeries(:,2),thresholdPercentiles(2));
dt = tsEvaGetTimeStep(inputTimeStamps);
minPeakDistance = minPeakDistanceInDays(2)/dt;
[pks2,indx2] = findpeaks(inputTimeSeries(:,2),'MinPeakDistance',minPeakDistance,'MinPeakHeight',thrsdt);
TimeStamps2=inputTimeStamps(indx2);
idc1=ones(size(TimeStamps1,1),1);
idc2=ones(size(TimeStamps2,1),1)*2;
idtotal=[idc1;idc2];
TimeStampstotal=[TimeStamps1;TimeStamps2];
pkstotal=[pks1;pks2];
[TimeStampstotal,iTimeStampstotal]=sort(TimeStampstotal,'ascend');
pkstotal=pkstotal(iTimeStampstotal);
idtotal=idtotal(iTimeStampstotal);
indexk=[];
aa=0;
for ii=1:length(TimeStampstotal)-1
    if aa==1
        aa=0;
        continue
    end
    if idtotal(ii)==1
        iix= find(idtotal(ii+1:end)==2,1);
        iix2= find(idtotal(ii+1:end)==1,1);
        if isempty(iix) || isempty(iix2)
            continue
        end
        if (TimeStampstotal(ii+iix2)-TimeStampstotal(ii))>minPeakDistanceInDays(1)...
                && (TimeStampstotal(ii+iix)-TimeStampstotal(ii))<maxDistanceMultivariatePeaksInDays(1)
            if (pkstotal(ii))>thresholds(1)...
                    && (pkstotal(ii+iix))>thresholds(2)
                indexk=[indexk,[ii,ii+iix]];
                aa=1;
            end
        else
            
            continue
        end
    elseif idtotal(ii)==2
        iix= find(idtotal(ii+1:end)==1,1);
        iix2= find(idtotal(ii+1:end)==2,1);
        if isempty(iix2)
            continue
        end
        if (TimeStampstotal(ii+iix2)-TimeStampstotal(ii))>minPeakDistanceInDays(2)...
                && (TimeStampstotal(ii+iix)-TimeStampstotal(ii))<maxDistanceMultivariatePeaksInDays(1)
            if (pkstotal(ii))>thresholds(2)...
                    && (pkstotal(ii+iix))>thresholds(1)
                indexk=[indexk,[ii,ii+iix]];
                aa=1;
            end
        else
            
            continue
        end
    end
    
    
end

idg=idtotal(indexk);
Timeg=TimeStampstotal(indexk);
pksg=pkstotal(indexk);

id1=find(idg==1);
TimeStamps1f=Timeg(id1);
pks1f=pksg(id1);

id2=find(idg==2);
TimeStamps2f=Timeg(id2);
pks2f=pksg(id2);


jointExtremes=cat(3,[TimeStamps1f,TimeStamps2f],[pks1f,pks2f]);





