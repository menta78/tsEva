function [jointextremes,jointextremes2, thresholds,timestampstotal,pkstotal] = tsCopulaSampleJointPeaks(inputtimestamps, ...
    inputtimeseries, ...
    thresholdpercentiles, ...
    minpeakdistanceindays, ...
    maxdistancemultivariatepeaksindays)

% this function samples the joint extremes (or peaks) from a multivariate
% inputtimeseries, based on a multivariate pot.
% this function assumes the input time series are stationary, so you need
% to feed it with the transformed time series.
% Note that a joint event is defined as an event that occurs almost concurrently
% among the considered series. The code specifically handles overlaps. 
% If an overlap of joint peak events (meaning an event where all the series exhibit 
% exceedance above their respective thresholds) occur with joint non-peak events 
% (meaning an event that not all of the series exhibit exceedance above their respective thresholds), the overlapping non-peak events get discarded. If an 
%overlap of two or more joint non-peak events occur, the one having the largest mean value is kept and the others are discarded. 
%Also, if an overlap of two or more joint peak events occur, the one having the largest mean value is kept and others get discarded. 
% The functioning of the code is in a way that the search radius to define overlaps is about two times larger than average value set for maxdistancemultivariatepeaksindays.
%maxdistancemultivariatepeaksindays is organized as first, second, and
%third value refer to distance between combination of first and second,
%third and second, and first and third series in case of trivariate input.
%in case of bivariate input, both entries have to be identical.

% input data:
% - inputtimestamps: 1d array with length nt, time stamps for the input
%   time series. must be the same for all the time series
% - inputtimeseries: data of each time series. 2d array with size [nt x n], where n is the number of
%   variables we need to analyze. Can handle bivariate or trivariate cases.
% - thresholdpercentiles: percentile threshold to be used for each time series.
%   1d array with length n, where n is the number of variables to be
%   considered.
% - minpeakdistanceindays: minimum time distance among peaks of the same
%   variable. 1d array with length n, where n is the number of variables to
%   be considered.
% - maxdistancemultivariatepeaksindays: maximum time distance among peaks
%   of different variables for the peaks to be considered joint. 1d array with length
%   n, where n is the number of variable to be considered. In case of bivariate inputtimeseries,
%   the two entries for this parameter have to be identical. In case of trivariate inputtimeseries, 
%   the first entry represent time distance between first and second series, the second entry represents the
%   time distance between second and third series, and the third entry represents the time distance between
%   first and third series. 

% output data:
% - jointextremes: array containing the joint extremes. 3d array with shape
%   [njointevents, n, 2], where njointevents is the number of found events,
%   n is the number of variables. jointextremes[ijointevent, ivariable, 1]
%   contains the time stamp of the peak related to the ijointevent-th event
%   and ivariable-th variable, while jointextremes[ijointevent, ivariable, 2]
%   contains the value of the peak.
% - jointextremes2: array containing the joint extremes for non-peak events. 3d array with shape
%   [njointevents, n, 2], where njointevents is the number of found events,
%   n is the number of variables. jointextremes[ijointevent, ivariable, 1]
%   contains the time stamp of the peak related to the ijointevent-th event
%   and ivariable-th variable, while jointextremes[ijointevent, ivariable, 2]
%   contains the value of the peak.
% - thresholds: 1d array with size n, containing the value of the threshold
%   for each variable.
% - timestampstotal: 1d array with size n, containing the time stamps of all 
%   peaks of the monovariates.
% - pkstotal: 1d array with size n, containing the values of all peaks of the
%   monovariates. 




if size(inputtimeseries,2)==3  %in the case of trivariate input

    thrsh1=prctile(inputtimeseries(:,1),thresholdpercentiles(1));
    thrsh2=prctile(inputtimeseries(:,2),thresholdpercentiles(2));
    thrsh3=prctile(inputtimeseries(:,3),thresholdpercentiles(3));
    thresholds=[thrsh1,thrsh2,thrsh3];


    dt = tsEvaGetTimeStep(inputtimestamps); %this is only calculated once, assuming inputtimestamps are similar in between the input series
    minpeakdistance1 = minpeakdistanceindays(1)/dt;
    [pks1,indx1] = findpeaks(inputtimeseries(:,1),'minpeakdistance',minpeakdistance1); %minpeakheight was not set to keep information about peaks below threshold for its later use
    timestamps1=inputtimestamps(indx1);

    % plot(datetime(datevec(timestamps1)),pks1,'k.')
    minpeakdistance2= minpeakdistanceindays(2)/dt;
    [pks2,indx2] = findpeaks(inputtimeseries(:,2),'minpeakdistance',minpeakdistance2);
    timestamps2=inputtimestamps(indx2);

    % plot(datetime(datevec(timestamps2)),pks2,'g.')
    minpeakdistance3 = minpeakdistanceindays(3)/dt;
    [pks3,indx3] = findpeaks(inputtimeseries(:,3),'minpeakdistance',minpeakdistance3);
    timestamps3=inputtimestamps(indx3);
    % plot(datetime(datevec(timestamps3)),pks3,'r.')

    timestampstotal=[timestamps1;timestamps2;timestamps3];
    pkstotal=[pks1;pks2;pks3];
    [timestampstotal,itimestampstotal]=sort(timestampstotal,'ascend'); %this ensures a time-sorted combination of series
    pkstotal=pkstotal(itimestampstotal);

    idc1=ones(size(timestamps1,1),1);   %this is a tag to assign to series, 1 for first series, 2 for second, and 3 for third
    idc2=ones(size(timestamps2,1),1)*2;
    idc3=ones(size(timestamps3,1),1)*3;
    idtotal=[idc1;idc2;idc3];
    idtotal=idtotal(itimestampstotal);

    indexk=[]; %to store information of joint peak events
    indexk2=[]; %to store information of joint non-peak events where not all series exhibit exceedance above their respective thresholds but at least one does

    for ii=1:length(timestampstotal)-1

        if idtotal(ii)==1
            iix3=find(timestampstotal(ii+1:end)-timestampstotal(ii)<=maxdistancemultivariatepeaksindays(3)&idtotal(ii+1:end)==3);%looks for the next peaks belonging to the third series that are also within time distance not further than maximum allowed time distance
            iix2=find(timestampstotal(ii+1:end)-timestampstotal(ii)<=maxdistancemultivariatepeaksindays(1)&idtotal(ii+1:end)==2);% same as above but for the second series
           
        elseif idtotal(ii)==2
            iix3=find(timestampstotal(ii+1:end)-timestampstotal(ii)<=maxdistancemultivariatepeaksindays(1)&idtotal(ii+1:end)==1);%looks for the next peaks belonging to the first series that are also within time distance not further than maximum allowed time distance
            iix2=find(timestampstotal(ii+1:end)-timestampstotal(ii)<=maxdistancemultivariatepeaksindays(2)&idtotal(ii+1:end)==3);% same as above but for the third series
            
        elseif idtotal(ii)==3
            iix3=find(timestampstotal(ii+1:end)-timestampstotal(ii)<=maxdistancemultivariatepeaksindays(3)&idtotal(ii+1:end)==1);%looks for the next peaks belonging to the first series that are also within time distance not further than maximum allowed time distance
            iix2=find(timestampstotal(ii+1:end)-timestampstotal(ii)<=maxdistancemultivariatepeaksindays(2)&idtotal(ii+1:end)==2);% same as above but for the second series
            
        end

        if isempty(iix2) || isempty(iix3) %if either iix2 or iix3 is empty, the next iteration needs to be followed
            continue
        end
        indcombination=nchoosek([ii+iix3;ii+iix2;ii],3); %all possible combinations that could be formed by indices sampled
        indseries=idtotal(indcombination); %this refer to indices of the series, i.e.,1 1 for first series, 2 for second series, and so on. however, not all combinations necessarily have one sample taken from all series and it is possible that one series have contributed more than one value in the combination (i.e., 2 2 3 or 1 1 2)
        want(1,1,1:3) = [1 2 3]; %this is to make sure that sampled combinations have one value taken from each series and discarding combinations that have two or more samples from one series only
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
        thrgr=thresholds(indseries);
        if any(all(pksgroup>=thrgr,2)) %checks if at least one combination exhibit features of a joint peak event
            indextojointpeaks=find((all(pksgroup>=thrgr,2)));
            if size(indextojointpeaks,1)==1  %meaning only one such combinations exists
                indexk=[indexk,indcombination(indextojointpeaks,:)];
            elseif size(indextojointpeaks,1)>1 %mening more than one such combinations exist
                indjointpeaks=indcombination(indextojointpeaks,:);
                [~,maxindex]=max(mean(pksgroup(indextojointpeaks,:),2)); %the combination with largest mean values is kept (among the joint peak events)
                indexk=[indexk,indjointpeaks(maxindex,:)];
            end
        elseif any(any(pksgroup>=thrgr,2))  %checks if at least one combination exhibit features of a joint non-peak event
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

    id1=find(idg==1);
    timestamps1f=timeg(id1);
    pks1f=pksg(id1);

    id2=find(idg==2);
    timestamps2f=timeg(id2);
    pks2f=pksg(id2);

    id3=find(idg==3);
    timestamps3f=timeg(id3);
    pks3f=pksg(id3);

    idg2=idtotal(indexk2);
    timeg2=timestampstotal(indexk2);   %this samples all joint non-peak events from the total series
    pksg2=pkstotal(indexk2);

    id1=find(idg2==1);
    timestamps1f2=timeg2(id1);
    pks1f2=pksg2(id1);

    id2=find(idg2==2);
    timestamps2f2=timeg2(id2);
    pks2f2=pksg2(id2);

    id3=find(idg2==3);
    timestamps3f2=timeg2(id3);
    pks3f2=pksg2(id3);


    eventstime=[timestamps1f,timestamps2f,timestamps3f]; %creates an n*3 matrix of all combination of events (their timestamps)  that were sampled, however, it is possible that some combinations
                                                         % have overlaps (or duplicates). Then, these duplicates are detected and discarded through comparing each
                                                         % combination with another overlapping combination and selecting the one that has higher mean value
                                                         %this appears to require a while loop as overlaps does not necessarily occur only between two neighboring combinations so the loop continues until
                                                         %there is an overlap
    eventpeaks=[pks1f,pks2f,pks3f];
    SZremoved=[];
    Szbefore=size(eventpeaks,1);
    while   any((min(eventstime(2:end,:),[],2)-max(eventstime(1:end-1,:),[],2))<=0)
        indextoremove=[];

        for jx=2:length(eventstime)
            event0=eventstime(jx-1,:);
            event1=eventstime(jx,:);


            if min(event1,[],2)<=max(event0,[],2)
                % disp('joint peak event')
                ww0=numel(find(mean(eventpeaks(jx-1,:))>=mean(eventpeaks(jx,:))));
                ww1=numel(find(mean(eventpeaks(jx,:))>=mean(eventpeaks(jx-1,:))));

                if ww0>ww1
                    indextoremove=[indextoremove,jx];
                else
                    indextoremove=[indextoremove,jx-1];

                end

            end
        end
        SZremoved=[SZremoved,size(unique(indextoremove),2)];
        eventstime(indextoremove,:)=[];
        eventpeaks(indextoremove,:)=[];
    end
disp([num2str(round(((sum(SZremoved)/Szbefore)*100)*10)/10),'%',' discarded from joint peak events'])
    jointextremes=cat(3,eventstime,eventpeaks); %this holds all information about joint peak events


    eventstime2=[timestamps1f2,timestamps2f2,timestamps3f2];
    eventpeaks2=[pks1f2,pks2f2,pks3f2];
    SZremoved=[];
    Szbefore=size(eventpeaks2,1);
    while   any((min(eventstime2(2:end,:),[],2)-max(eventstime2(1:end-1,:),[],2))<=0) %this part deletes overlaps for joint non-peak events the same way
        %it was deleted for joint peak events
        indextoremove2=[];
        for jx=2:length(eventstime2)
            event0=eventstime2(jx-1,:);
            event1=eventstime2(jx,:);


            if min(event1,[],2)<=max(event0,[],2)
                % disp('joint-nonpeak event')

                ww0=numel(find(mean(eventpeaks2(jx-1,:))>=mean(eventpeaks2(jx,:))));
                ww1=numel(find(mean(eventpeaks2(jx,:))>=mean(eventpeaks2(jx-1,:))));
                if ww0>ww1
                    indextoremove2=[indextoremove2,jx];
                else
                    indextoremove2=[indextoremove2,jx-1];

                end

            end
        end
              SZremoved=[SZremoved,size(unique(indextoremove2),2)];

        eventstime2(indextoremove2,:)=[];
        eventpeaks2(indextoremove2,:)=[];
    end

    disp([num2str(round(((sum(SZremoved)/Szbefore)*100)*10)/10),'%',' discarded from joint non peak events'])

 mintnonpeak=min(eventstime2,[],2); %an array with minimum time of each non-joint peak event
          maxtnonpeak=max(eventstime2,[],2);
indextoremove3=[];
SZremoved=[];
    for indpeak=1:size(eventpeaks,1) %this loop for every joint peak events look for overlapping non-peak joint events and discard
        %it was deleted for joint peak events
       
          event0=eventstime(indpeak,:);
          
        ind0=(min(event0)>maxtnonpeak|max(event0)<mintnonpeak);
        if find(eventstime2(~ind0))
%disp('overlap non-peak joint event with peak joint event')
indextoremove3=[indextoremove3,find(~ind0)'];
        end
      
    end
                  SZremoved=[SZremoved,size(unique(indextoremove3),2)];
    disp([num2str(round(((sum(SZremoved)/Szbefore)*100)*10)/10),'%',' further discarded from joint non peak events due to overlap with peak events'])

  eventstime2(indextoremove3,:)=[];
        eventpeaks2(indextoremove3,:)=[];


    jointextremes2=cat(3,eventstime2,eventpeaks2);

elseif size(inputtimeseries,2)==1
    disp('only one series- needs to be at least two')
elseif size(inputtimeseries,2)==2
    thrsh1=prctile(inputtimeseries(:,1),thresholdpercentiles(1));
    thrsh2=prctile(inputtimeseries(:,2),thresholdpercentiles(2));

    thresholds=[thrsh1,thrsh2];


    dt = tsEvaGetTimeStep(inputtimestamps);
    minpeakdistance1 = minpeakdistanceindays(1)/dt;
    [pks1,indx1] = findpeaks(inputtimeseries(:,1),'minpeakdistance',minpeakdistance1);
    timestamps1=inputtimestamps(indx1);
    % plot(datetime(datevec(timestamps1)),pks1,'k.')
    minpeakdistance2= minpeakdistanceindays(2)/dt;
    [pks2,indx2] = findpeaks(inputtimeseries(:,2),'minpeakdistance',minpeakdistance2);
    timestamps2=inputtimestamps(indx2);
    % plot(datetime(datevec(timestamps2)),pks2,'g.')
    timestampstotal=[timestamps1;timestamps2];
    pkstotal=[pks1;pks2];
    [timestampstotal,itimestampstotal]=sort(timestampstotal,'ascend');
    pkstotal=pkstotal(itimestampstotal);

    idc1=ones(size(timestamps1,1),1);
    idc2=ones(size(timestamps2,1),1)*2;
    idtotal=[idc1;idc2];
    idtotal=idtotal(itimestampstotal);

    indexk=[];
    indexk2=[];

    for ii=1:length(timestampstotal)-1

        if idtotal(ii)==1
            iix2=find(timestampstotal(ii+1:end)-timestampstotal(ii)<=maxdistancemultivariatepeaksindays(1)&idtotal(ii+1:end)==2);% same as above but for the second series
            if isempty(iix2)
                continue
            end
        elseif idtotal(ii)==2
            iix2=find(timestampstotal(ii+1:end)-timestampstotal(ii)<=maxdistancemultivariatepeaksindays(1)&idtotal(ii+1:end)==1);%looks for the next peaks belonging to the thirs series that are also within time distance not further than maximum allowed time distance
            if isempty(iix2)
                continue
            end



        end
        indcombination=nchoosek([ii+iix2;ii],2); %all possible combinations that could be formed by indices sampled
        indseries=idtotal(indcombination); %this refer to indices of the series, i.e.,1 1 for first series, 2 for second series, and so on. however, not all combinations necessarily have one sample taken from all series and it is possible that one series have contributed more than one value in the combination (i.e., 2 2 3 or 1 1 2) but this has to be like 1 2 3
        want(1,1,1:2) = [1 2]; %this is to make sure that combinations have one sample from each series and discarding combinations that have two or more samples from one series only
        if size(indseries,2)==1 %a column vector needs to change to a row vector for proper functioning of next lines
            indseries=indseries';
        end
        indextodesiredrows = all(any(bsxfun(@eq,indseries,want),2),3);
        rownumbers = find(indextodesiredrows);
        indcombination=indcombination(rownumbers,:);
        pksgroup=pkstotal(indcombination);
        if size(pksgroup,2)==1 %a column vector needs to change to a row vector for proper functioning of next lines
            pksgroup=pksgroup';
        end
        indseries=indseries(rownumbers,:);
        thrgr=thresholds(indseries);
        if any(all(pksgroup>=thrgr,2)) %checks if all peaks or at least one combination are above their respective thresholds (i.e., joint peak event)
            indextojointpeaks=find((all(pksgroup>=thrgr,2)));
            if size(indextojointpeaks,1)==1  %meaning only one combination exists
                indexk=[indexk,indcombination(indextojointpeaks,:)];
            elseif size(indextojointpeaks,1)>1 %mening more than one combination exist
                indjointpeaks=indcombination(indextojointpeaks,:);
                [~,maxindex]=max(mean(pksgroup(indextojointpeaks,:),2)); %the combination with largest mean values is kept (among the joint peaks)
                indexk=[indexk,indjointpeaks(maxindex,:)];
            end
        elseif any(any(pksgroup>=thrgr,2))  %non-peak joint events
            indextojointnpeaks=find((any(pksgroup>=thrgr,2)));
            if size(indextojointnpeaks,1)==1  %meaning only one combination exists
                indexk2=[indexk2,indcombination(indextojointnpeaks,:)];
            elseif size(indextojointnpeaks,1)>1 %mening more than one combination exist
                indjointnpeaks=indcombination(indextojointnpeaks,:);
                [~,maxindex]=max(mean(pksgroup(indextojointnpeaks,:),2)); %the combination with largest mean values is kept (among the joint peaks)
                indexk2=[indexk2,indjointnpeaks(maxindex,:)];
            end
        end





    end

    idg=idtotal(indexk);
    timeg=timestampstotal(indexk);
    pksg=pkstotal(indexk);

    id1=find(idg==1);
    timestamps1f=timeg(id1);
    pks1f=pksg(id1);

    id2=find(idg==2);
    timestamps2f=timeg(id2);
    pks2f=pksg(id2);


    idg2=idtotal(indexk2);
    timeg2=timestampstotal(indexk2);
    pksg2=pkstotal(indexk2);

    id1=find(idg2==1);
    timestamps1f2=timeg2(id1);
    pks1f2=pksg2(id1);



    id2=find(idg2==2);
    timestamps2f2=timeg2(id2);
    pks2f2=pksg2(id2);


    eventstime=[timestamps1f,timestamps2f];
    eventpeaks=[pks1f,pks2f];
    SZremoved=[];
    Szbefore=size(eventpeaks,1);
    while   any((min(eventstime(2:end,:),[],2)-max(eventstime(1:end-1,:),[],2))<=0)
        indextoremove=[];
        for jx=2:length(eventstime)
            event0=eventstime(jx-1,:);
            event1=eventstime(jx,:);


            if min(event1,[],2)<=max(event0,[],2)
                %  disp('joint peak event')
                ww0=numel(find(mean(eventpeaks(jx-1,:))>=mean(eventpeaks(jx,:))));
                ww1=numel(find(mean(eventpeaks(jx,:))>=mean(eventpeaks(jx-1,:))));
                if ww0>ww1
                    indextoremove=[indextoremove,jx];
                else
                    indextoremove=[indextoremove,jx-1];

                end

            end
        end
        SZremoved=[SZremoved,size(unique(indextoremove),2)];
        eventstime(indextoremove,:)=[];
        eventpeaks(indextoremove,:)=[];
    end
    disp([num2str(round(((sum(SZremoved)/Szbefore)*100)*10)/10),'%',' discarded from joint peak events'])

    jointextremes=cat(3,eventstime,eventpeaks);

    eventstime2=[timestamps1f2,timestamps2f2];
    eventpeaks2=[pks1f2,pks2f2];
    SZremoved=[];
    Szbefore=size(eventpeaks2,1);
    while   any((min(eventstime2(2:end,:),[],2)-max(eventstime2(1:end-1,:),[],2))<=0) %this part deletes overlaps for joint non-peak events the same way
        %it was deleted for joint peak events
        indextoremove2=[];
        for jx=2:length(eventstime2)
            event0=eventstime2(jx-1,:);
            event1=eventstime2(jx,:);


            if min(event1,[],2)<=max(event0,[],2)
                %  disp('joint-nonpeak event')
                ww0=numel(find(mean(eventpeaks2(jx-1,:))>=mean(eventpeaks2(jx,:))));
                ww1=numel(find(mean(eventpeaks2(jx,:))>=mean(eventpeaks2(jx-1,:))));
                if ww0>ww1
                    indextoremove2=[indextoremove2,jx];
                else
                    indextoremove2=[indextoremove2,jx-1];

                end

            end
        end
        SZremoved=[SZremoved,size(unique(indextoremove2),2)];
        eventstime2(indextoremove2,:)=[];
        eventpeaks2(indextoremove2,:)=[];
    end
    disp([num2str(round(((sum(SZremoved)/Szbefore)*100)*10)/10),'%',' discarded from joint non peak events'])

    mintnonpeak=min(eventstime2,[],2); %an array with minimum time of each non-joint peak event
    maxtnonpeak=max(eventstime2,[],2);
    indextoremove3=[];
    SZremoved=[];
    for indpeak=1:size(eventpeaks,1) %this loop for every joint peak events look for overlapping non-peak joint events and discard
        %it was deleted for joint peak events

        event0=eventstime(indpeak,:);

        ind0=(min(event0)>maxtnonpeak|max(event0)<mintnonpeak);
        if find(eventstime2(~ind0))
            %disp('overlap non-peak joint event with peak joint event')
            indextoremove3=[indextoremove3,find(~ind0)'];
        end

    end
    SZremoved=[SZremoved,size(unique(indextoremove3),2)];
    disp([num2str(round(((sum(SZremoved)/Szbefore)*100)*10)/10),'%',' further discarded from joint non peak events due to overlap with peak events'])

    eventstime2(indextoremove3,:)=[];
    eventpeaks2(indextoremove3,:)=[];




    jointextremes2=cat(3,eventstime2,eventpeaks2);


end





