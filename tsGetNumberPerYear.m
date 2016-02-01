function nperYear=tsGetNumberPerYear(ms,locs)
% function nperYear=tsGetNumberPerYear(ms,locs)
% 
% Gives number of events per year from a time series ms (sd and values vertical vectors) and sporadic indices of events (locs)
% 
% Michalis Vousdoukas 2015


%% make the full time vector
% get all years
sdfull=nanmin(ms(:,1)):nanmax(ms(:,1));

% reset
sdfull2=sdfull-nanmin(sdfull);

% round to years
[sdroundfull,dvec,sdunique,dvunique]=tsRoundSDate(sdfull2,1);


[yearss,~,~] = unique(dvec(:,1),'rows');

yearss=yearss+1;

%% prepare the series time vector

% reset
sdser=ms(:,1)-nanmin(ms(:,1));

% round to years
[sdroundser,dvecser,sduniqueser,dvuniqueser]=tsRoundSDate(sdser,1);

[yearsser,~,icmser] = unique(dvecser(:,1),'rows');

yearsser=yearsser+1;


%% 


    
sdp=ms(locs,1);
% make years vector
[sdround,dvec,sdunique,dvunique]=tsRoundSDate(sdp-min(sdfull),1);
[yearssp,~,icm] = unique(dvec(:,1),'rows');

yearssp=yearssp+1;

% get number of events per year

ecoount=accumarray(icm,sdp,[],@(x) length(x)); 

% get number of measurements per year

ecoountSeries=accumarray(icmser,sdroundser,[],@(x) length(x)); 

% exclude years with few measurements (less than 50%)
iexcl=ecoountSeries/max(ecoountSeries)<0.5;


nperYear=nan(length(yearss),1);
nperYear(yearsser)=0;
nperYear(yearssp)=ecoount;
nperYear(yearsser(iexcl))=nan;

