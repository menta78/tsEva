function [sdround,dvecm,sdunique,dvunique]=tsRoundSDate(sd,sd_precision)
% function [sdround,dvec,sdunique,dvunique]=tsRoundSDate(sd,sd_precision)
% ROunds sd values according to desired precision, years, months, days, etc : 1-4
% 
% Michalis Vousdoukas 2015

dvecm=datevec(sd);

dvecm(:,sd_precision+1)=1;
dvecm(:,sd_precision+2:end)=0;

sdround=datenum(dvecm);

sdunique=unique(sdround);

[dvunique,~,~] = unique(dvecm,'rows');

