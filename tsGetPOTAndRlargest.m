function [POTdata,Rlargest_data]=tsGetPOTAndRlargest(ms,pcts,desiredEventsPerYear, varargin)
% function POTdata=tsGetPOTAndRlargest(ms,pcts,desiredEventsPerYear)
% Gets POT using an automatic threshold such that the mean number of events per year is equal to desiredEventsPerYear
% 
% INPUTS:
% ms  data in two columns sdate and values. IT IS ASSUMED
% pcts    vestor of percentiles tested
% desiredEventsPerYear    mean number of events per year
% 
% Michalis Vousdoukas, Evangelos Voukouvalas, Lorenzo Mentaschi 2015


% get number of years
% nyears=(nanmax(ms(:,1))-nanmin(ms(:,1)))/365;

args.minPeakDistanceInDays = -1;
args = tsEasyParseNamedArgs(varargin, args);
minPeakDistanceInDays = args.minPeakDistanceInDays;
if minPeakDistanceInDays == -1
    error('label parameter ''minPeakDistanceInDays'' must be set')
end

dt = tsEvaGetTimeStep(ms(:,1));
minPeakDistance = minPeakDistanceInDays/dt;

nyears=nanmin(diff(ms(:,1)))*length(ms(~isnan(ms(:,2)),1))/365;

if length(pcts) == 1
    % testing at least 2 percentages, to be able to compute the error on
    % the percentage.
    pcts = [pcts(1) - 3, pcts(1)];
    % if there is only one percentile that means that the user does not
    % want to look for n peaks every year. Setting therefore
    % desiredEventsPerYear to -1;
    desiredEventsPerYear = -1;
end

numperyear=nan(length(pcts),1);
minnumperyear=nan(length(pcts),1);
thrsdts=nan(length(pcts),1);

for ipp=1:length(pcts)

    disp(['Finding optimal threshold ' num2str(100*ipp/length(pcts)) '%...']);

    thrsdt=prctile(ms(:,2),pcts(ipp));
    thrsdts(ipp) = thrsdt;

    [pks,locs] = findpeaks(ms(:,2),'MinPeakDistance',minPeakDistance,'MinPeakHeight',thrsdt);

    %     for POT
    numperyear(ipp)=length(pks)/nyears;

%     for R-largest
    nperYear=tsGetNumberPerYear(ms,locs);
    minnumperyear(ipp)=nanmin(nperYear);

    if (ipp > 1) && (length(pks)/nyears<desiredEventsPerYear) && (nanmin(nperYear)<desiredEventsPerYear)

        break

    end
end

% evaluating the error on the threshold
diffNPerYear = nanmean(diff(flipud(numperyear)));
if diffNPerYear == 0
  diffNPerYear = 1;
end
thresholdError = nanmean(diff(thrsdts)/diffNPerYear)/2;

%% Use optimal for POT
indexp=nanmax(find(numperyear>desiredEventsPerYear==1));

try
   
if ~isempty(indexp)
    thrsd = prctile(ms(:,2),pcts(indexp)); % threshold is the 99.9th percentile of the timeseries;
else
    thrsd = 0;
end

% figure(1)
% plot(pcts,numperyear,pcts,minnumperyear,'r')
% 
% figure;plot(ms(:,1),ms(:,2));
% horizontal_lines(thrsd)


[pks,locs] = findpeaks(ms(:,2),'MinPeakDistance',minPeakDistance,'MinPeakHeight',thrsd);

POTdata.threshold=thrsd;
POTdata.thresholdError = thresholdError;
POTdata.percentile=pcts(indexp);
POTdata.peaks=pks;
POTdata.ipeaks=locs;
POTdata.sdpeaks=ms(locs,1);

%% Use optimal for R-largest

indexp=nanmax(find(minnumperyear>desiredEventsPerYear));

if ~isempty(indexp)
    thrsd=prctile(ms(:,2),pcts(indexp)); % threshold is the 99.9th percentile of the timeseries;

    % figure(1)
    % plot(pcts,numperyear)
    % 
    % figure;plot(ms(:,1),ms(:,2));
    % horizontal_lines(thrsd)


    [pks,locs] = findpeaks(ms(:,2),'MinPeakDistance',minPeakDistance,'MinPeakHeight',thrsd);
    
else
    
    [pks,locs] = findpeaks(ms(:,2),'MinPeakDistance',minPeakDistance);
    
end


%%  make table with n number events

% get all years
sdfull=nanmin(ms(:,1)):nanmax(ms(:,1));

% reset
sdfull2=sdfull-nanmin(sdfull);

% round to years
[sdround,dvec,sdunique,dvunique]=tsRoundSDate(sdfull2,1);


[yearss,~,~] = unique(dvec(:,1),'rows');

yearss=yearss+1;

    
sdp=ms(locs,1);
% make years vector
[sdround,dvecp,sdunique,dvunique]=tsRoundSDate(sdp-min(sdfull),1);
[yearssp,~,icm] = unique(dvecp(:,1),'rows');

yearssp=yearssp+1;

if desiredEventsPerYear > -1
  sdTable=nan(length(yearss),desiredEventsPerYear);
  valTable=nan(length(yearss),desiredEventsPerYear);
  indTable=nan(length(yearss),desiredEventsPerYear);

  for i=1:length(yearss)

      ii=dvecp(:,1)==yearss(i)-1;

      vals=pks(ii);
      sds=sdp(ii);

      [v2,isort]=sort(vals,'descend');

      s2=sds(isort);

      inds=locs(isort);

      if length(inds)>=desiredEventsPerYear

          sdTable(i,:)=s2(1:desiredEventsPerYear);
          valTable(i,:)=v2(1:desiredEventsPerYear);
          indTable(i,:)=inds(1:desiredEventsPerYear);

      end
  end

  %% Export

  Rlargest_data.threshold=thrsd;
  Rlargest_data.percentile=pcts(indexp);
  Rlargest_data.peaks=valTable;
  Rlargest_data.ipeaks=sdTable;
  Rlargest_data.sdpeaks=sdTable;
else
  Rlargest_data = [];
end
  
catch exc
    dbstop at 148;
    disp(getReport(exc));
end
