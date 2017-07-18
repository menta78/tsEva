function [nx,ny,addedy]=tsLinearExtrapolation(x,y,extrax,npoints)
% 
% Linearly extrapolates time series x,y to points extrax, according to the slope of the first/last number of points of the series (npoints)
% 
% Michalis Vousdoukas 2016

[x,ii]=sort(x);
y=y(ii);

x=x(:);
y=y(:);

extrax=extrax(:);

nx=x;
ny=y;

clear x y
if sum(extrax<nanmin(nx))>0


    dvv=nanmean(diff(ny(1:npoints+1)));
    dTr=nanmean(diff(nx(1:npoints+1)));

    addedx=nx(1)-extrax(extrax<nanmin(nx));

    addedy0=ny(1)-dvv.*addedx./dTr;


    ny=[addedy0 ; ny];
    nx=[extrax(extrax<nanmin(nx)) ; nx];

    addedy{1}=addedy0;
    
end

if sum(extrax>nanmax(nx))>0

    dvv=nanmean(diff(ny(end-npoints:end)));
    dTr=nanmean(diff(nx(end-npoints:end)));

    addedx=nx(end)-extrax(extrax>nanmax(nx));

    addedy0=ny(end)-dvv.*addedx./dTr;


    ny=[ny ; addedy0];
    nx=[nx ; extrax(extrax>nanmax(nx))];
    
    addedy{2}=addedy0;
end