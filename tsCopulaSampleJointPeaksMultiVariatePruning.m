function [jointextremes,jointextremes2,thresholdsC,jointExtremeIndices,peakIndicesAll] = ...
    tsCopulaSampleJointPeaksMultiVariatePruning(inputtimestamps,...
    inputtimeseries,...
    thresholdpercentiles,...
    minpeakdistanceindays,...
    maxdistancemultivariatepeaksindays)

% This function samples the compound events (multivariate peaks) from
% inputtimeseries, based on a multivariate peak over threshold approach.
% input time series must be stationarized beforehand. A compound peak event
% is defined as an event whereby all monovariates exceed their respective
% thresholds while a compound non-peak event occurs whereby not all of the
% monovariates exceed their respecive threshold but at least one
% monovariate does
%The algorithm of the code is such that initially every possible
%combination of compound peak and non-peak events are sampled. In this
%sampling the there are two important parameters to be
%set:minpeakdistanceindays which sets the distance between monovariate
%peaks and maxdistancemultivariatepeaksindays which sets the distance (for
%every possible combination of multivariates which has to have a format and
%size matching size(nchoosek([1:numvar],2),1) where numvar is 2 in
%bivariate case, 3 in trivariate case, and so on). For pruning of
%overlapping events, events are sorted from the largest mean value to the
%lowest mean value and any overlapping event gets discarded as the code
%works out through in the decreasing order of events magnitude

% input data:
% - inputtimestamps: 1d array with length nt, time stamps for the input
%   time series. must be the same for all the time series
% - inputtimeseries: data of each time series. 2d array with size [nt x n],
% where n is the number of
%   variables.
% - thresholdpercentiles: percentile threshold to be used for each time
% series.
%   1d array with length n, where n is the number of variables to be
%   considered.
% - minpeakdistanceindays: minimum time distance among peaks of the same
%   variable. 1d array with length n, where n is the number of variables to
%   be considered.
% - maxdistancemultivariatepeaksindays: maximum time distance among peaks
%   of different variables for the peaks to be considered joint. 1d array
%   with length size(nchoosek([1:n],2),1), where n is the number of
%   variables to be considered. nchoosek([1:n],2) shows the format in which
%   maxdistancemultivariatepeaksindays will be interpreted.

% output data:
% - jointextremes: array containing the joint extremes. 3d array with shape
%   [njointevents, n, 2], where njointevents is the number of found events,
%   n is the number of variables. jointextremes[ijointevent, ivariable, 1]
%   contains the time stamp of the peak related to the ijointevent-th event
%   and ivariable-th variable, while jointextremes[ijointevent, ivariable,
%   2] contains the value of the peak.
% - jointextremes2: Same as jointextremes but for the non-peak joint
%   events, with size[njointnonpeakevents, n, 2], where njointnonpeakevents
%   is the number of found non-peak joint events
% - thresholdsC: 1d array of thresholds used for selecting peaks of
%   monovariates, with size 1*n
% - peaksjointidx: indices of compound peak events
% - peaksjointidx2: indices of compound non-peak events

if size(maxdistancemultivariatepeaksindays,2)~=size(nchoosek([1:size(inputtimeseries,2)],2),1)
    maxdistancemultivariatepeaksindays=repmat(maxdistancemultivariatepeaksindays,1,size(nchoosek([1:size(inputtimeseries,2)],2),1));
end
inputtimeseriesCell=mat2cell(inputtimeseries,size(inputtimeseries,1),ones(1,size(inputtimeseries,2)));
thresholdsC=cellfun(@(x,y) prctile(x,y),inputtimeseriesCell,num2cell(thresholdpercentiles(1:size(inputtimeseriesCell,2))));
dt = tsEvaGetTimeStep(inputtimestamps);
minpeakdistancearray = minpeakdistanceindays(1:size(inputtimeseriesCell,2))/dt;
[pkscell,indxcell]=cellfun(@(x,y) findpeaks(x,'minpeakdistance',y),inputtimeseriesCell,num2cell(minpeakdistancearray),'UniformOutput',0);
timestampscell=cellfun(@(x,y) x(y),repmat({inputtimestamps},1,size(inputtimeseriesCell,2)),indxcell,'UniformOutput',0);

timestampstotal=cell2mat(timestampscell');
indxcelltotal=cell2mat(indxcell');
pkstotal=cell2mat(pkscell');
[timestampstotal,itimestampstotal]=sort(timestampstotal,'ascend'); %this ensures a time-sorted combination of series
pkstotal=pkstotal(itimestampstotal);
indxcelltotal=indxcelltotal(itimestampstotal);
varIdcell=num2cell(1:size(inputtimeseriesCell,2));
idccell=cellfun(@(x,y) y.*ones(size(x,1),1),timestampscell,varIdcell,'UniformOutput',0);

idtotal=cell2mat(idccell');
idtotal=idtotal(itimestampstotal);

indexk=[]; %to store information of joint peak events
indexk2=[]; %to store information of joint non-peak events
Numvar=unique(idtotal);
combvars=nchoosek([1:size(unique(idtotal),1)],2);
[~,iss]=sort(diff(combvars,[],2));
combvars=combvars(iss,:);

for ii=1:length(timestampstotal)-1


    Numvar2=Numvar(Numvar~=idtotal(ii));
    iixt={};
    for jj=1:length(Numvar2)
        A(1,1,1:length([idtotal(ii),Numvar2(jj)]))=[idtotal(ii),Numvar2(jj)];
        ixy=find(all(~all(bsxfun(@minus,combvars,A),2),3));
        iix3=find(timestampstotal(ii+1:end)-timestampstotal(ii)<=maxdistancemultivariatepeaksindays(ixy)&idtotal(ii+1:end)==Numvar2(jj));
        iixt=[iixt,{ii+iix3}];
    end
    if any(cellfun(@(x) isempty(x),iixt))
        continue
    end
    indcombination=nchoosek([cell2mat(iixt');ii],size(inputtimeseries,2)); %all possible combinations that could be formed by indices sampled
    indseries=idtotal(indcombination); %this refer to indices of the series, i.e.,1 1 for first series, 2 for second series, and so on. however, not all combinations necessarily have one sample taken from all series and it is possible that one series have contributed more than one value in the combination (i.e., 2 2 3 or 1 1 2)
    want(1,1,1:size(inputtimeseries,2)) = [1:size(inputtimeseries,2)]; %this is to make sure that sampled combinations have one value taken from each series and discarding combinations that have two or more samples from one series only
    if size(indseries,2)==1 %a column vector needs to change to a row vector for proper functioning of next lines
        indseries=indseries';
    end
    indextodesiredrows = all(any(bsxfun(@eq,indseries,want),2),3); %this ensure that only the combinations that have one value from each series are kept
    rownumbers = find(indextodesiredrows);
    indcombination=indcombination(rownumbers,:);
    pksgroup=pkstotal(indcombination);
    if size(pksgroup,2)==1 %a column vector needs to change to a row vector for proper functioning of next lines
        pksgroup=pksgroup';
    end
    indseries=indseries(rownumbers,:);

    thrgr=thresholdsC(indseries);
    if any(all(pksgroup>=thrgr,2)) %checks if at least one combination matches of a joint peak event
        indextojointpeaks=find((all(pksgroup>=thrgr,2)));
        if size(indextojointpeaks,1)==1  %meaning only one such combinations exists
            indexk=[indexk,indcombination(indextojointpeaks,:)];
        elseif size(indextojointpeaks,1)>1 %mening more than one such combinations exist
            indjointpeaks=indcombination(indextojointpeaks,:);
            [~,maxindex]=max(mean(pksgroup(indextojointpeaks,:),2)); %the combination with largest mean values is kept (among the joint peak events)
            indexk=[indexk,indjointpeaks(maxindex,:)];
        end
    elseif any(any(pksgroup>=thrgr,2))  %checks if at least one combination matches a joint non-peak event
        indextojointnpeaks=find((any(pksgroup>=thrgr,2)));
        if size(indextojointnpeaks,1)==1  %meaning only one such combinations exists
            indexk2=[indexk2,indcombination(indextojointnpeaks,:)];
        elseif size(indextojointnpeaks,1)>1 %mening more than one such combinations exist
            indjointnpeaks=indcombination(indextojointnpeaks,:);
            [~,maxindex]=max(mean(pksgroup(indextojointnpeaks,:),2)); %the combination with largest mean values is kept (among the joint non peak events)
            indexk2=[indexk2,indjointnpeaks(maxindex,:)];
        end
    end


end

idg=idtotal(indexk);
timeg=timestampstotal(indexk);  %this samples all joint peak events from the total series
pksg=pkstotal(indexk);
idg2=idtotal(indexk2);
timeg2=timestampstotal(indexk2);   %this samples all joint non-peak events from the total series
pksg2=pkstotal(indexk2);
eventstime=[];
eventpeaks=[];
eventstime2=[];
eventpeaks2=[];
eventid=[];eventid2=[];

indexks=reshape(indexk,size(inputtimeseries,2),length(indexk)/size(inputtimeseries,2));
[~,idgs]=sort(reshape(idg,size(inputtimeseries,2),length(indexk)/size(inputtimeseries,2)),1,'ascend');
indexkss=[];
for ij=1:size(idgs,2)
    indexkss=[indexkss,(indexks(idgs(:,ij),ij))'];
end
peaksjointid=reshape(indxcelltotal(indexkss),size(inputtimeseries,2),length(indxcelltotal(indexkss))/size(inputtimeseries,2));
peaksjointid=peaksjointid';

for jj=1:size(inputtimeseries,2)
    id1=find(idg==jj);
    timestamps1f=timeg(id1);
    pks1f=pksg(id1);
    idgf=indexk(id1);
    eventstime=[eventstime,timestamps1f];
    eventpeaks=[eventpeaks,pks1f];
    eventid=[eventid,idgf'];

    id1=find(idg2==jj);
    timestamps1f2=timeg2(id1);
    pks1f2=pksg2(id1);
    idg2f=indexk2(id1);
    eventstime2=[eventstime2,timestamps1f2];
    eventpeaks2=[eventpeaks2,pks1f2];
    eventid2=[eventid2,idg2f'];
end


%the next lines are concerned wih pruning of overlapping events

eventstimetotal=[eventstime;eventstime2];
eventpeakstotal=[eventpeaks;eventpeaks2];
idjn=[ones(size(eventstime,1),1);2*ones(size(eventstime2,1),1)];
eventidtotal=[eventid;eventid2];
[~,ids]=sort(mean(eventpeakstotal,2),'descend');
eventstimetotal=eventstimetotal(ids,:);
eventpeakstotal=eventpeakstotal(ids,:);
idjn=idjn(ids);
eventidtotal=eventidtotal(ids,:);
minarraypeakx=min(eventstimetotal,[],2);
maxarraypeakx=max(eventstimetotal,[],2);
indextoremove3=[];

for jx=1:size(eventstimetotal,1)


    if(any(indextoremove3==jx))
        continue
    end
    event0=eventstimetotal(jx,:);

    ind0=(min(event0)>maxarraypeakx|max(event0)<minarraypeakx);
    ind0(1:jx)=1;
    if find(eventstimetotal(~ind0))
        indextoremove3=[indextoremove3,find(~ind0)'];
    end

end
if isempty(indextoremove3)
    disp('no case for pruning found')
else
    disp([num2str(round(100*(size(unique(indextoremove3),2)/size(eventstimetotal,1))*10)/10),' %',' of sampled events pruned due to overlapping of events'])
end
eventstimetotal(indextoremove3,:)=[];
eventpeakstotal(indextoremove3,:)=[];
eventidtotal(indextoremove3,:)=[];
idjn(indextoremove3)=[];
eventstime=eventstimetotal(idjn==1,:);
eventstime2=eventstimetotal(idjn==2,:);
eventpeaks=eventpeakstotal(idjn==1,:);
eventpeaks2=eventpeakstotal(idjn==2,:);
eventid1=eventidtotal(idjn==1,:);
jointExtremeIndices=indxcelltotal(eventid1);
eventid22=eventidtotal(idjn==2,:);
peakIndicesAll=[indxcelltotal(eventid1);indxcelltotal(eventid22)];
disp([num2str(((size((eventstime),1)))),' ','compound peak events found'])
disp([num2str(((size((eventstime2),1)))),' ','compound non-peak events found'])
jointextremes=cat(3,eventstime,eventpeaks);
jointextremes2=cat(3,eventstime2,eventpeaks2);











