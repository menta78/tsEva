function [nx,ny,addedy]=tsLoglogExtrapolation(x,y,extrax,npoints)
% 
% Extrapolates exponentiallytime series x,y to points extrax, according to the slope of the first/last number of points of the series (npoints)
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

addedy=[];
if sum(extrax<nanmin(nx))>0

    if sum(ny<=0)>0
        
        offs=-min(ny)+1;
        
    else
        
        offs=0;
        
    end
    
    [nx0,ny0,addedy0]=tsLinearExtrapolation(log(nx),log(ny+offs),log(extrax(extrax<nanmin(nx))),npoints);
    
%     plot(nx,ny,nx0,ny0,'r--')
    nx=exp(nx0);
    ny=exp(ny0)-offs;
    addedy{1}=exp(addedy0{1})-offs;

%     ny=[addedy ; ny];
%     nx=[extrax ; nx];

end


if sum(extrax>nanmax(nx))>0

    if sum(ny<=0)>0
        
        offs=-min(ny)+1;
        
    else
        
        offs=0;
        
    end
    [nx0,ny0,addedy0]=tsLinearExtrapolation(log(nx),log(ny+offs),log(extrax(extrax>nanmax(nx))),npoints);
    
    nx=exp(nx0);
    ny=exp(ny0)-offs;
    addedy{2}=exp(addedy0{2})-offs;

%     ny=[y ; addedy];
%     nx=[x ; extrax];
    
end

if iscell(addedy)
    
    if length(addedy)==1
        
        addedy=addedy{1};
        
    end
    
end