function nvals = tsInterp1Extrap(X,V,Xq,logExtrap)
% function nvals=interp1Extrap(Tr,values,RP,logExtrap)
% 
% This function interpolates a series and extrapolates linearly or exponeniantly if necessary.
% Can also handle matrices (only for dependent variable values)
% Will also sort the Tr series if not monotonically increasing
% Tr-     dependent variable
% values- independent variable if not a vector it is recommended that entries are in rows
% RP- desired value
% logExtrap-  switch for log log extrapolation
%         
% Michalis Vousdoukas 2017



% Making sure that vector is horizontal
if sum(size(V)==1)>0
    
    V=V(:)';
    
else
    
    if find(size(V)==length(X))==1
        
        V=V';
        
    elseif find(size(V)==length(X))==2
        
        disp('No need to transpose')
        
    elseif length(find(size(V)==length(X))==1)>1
        
        warning('Make sure dependent variable has the right orientation')
        
    end
    
end

minX = min(min(Xq(:)), min(X(:)));
if minX < 1
  xoffset = -minX + 1;
else
  xoffset = 0;
end

X = X + xoffset;
Xq = Xq + xoffset;

% sort x independent variable
[X,ii]=sort(X);

V=V(:,ii);

clear nvals

npoints=1;

for i=1:size((V),1)
 
    
    if sum(Xq<min(X))>0 | sum(Xq>max(X))>0
        
        if length(Xq)==1

            extrax=unique([Xq(:)/2 ; Xq(:) ; Xq(:)*2]);
            
        else
            
            extrax=Xq;
            
        end
        
        if logExtrap==1
            
            [nx,ny,addedy]=tsLoglogExtrapolation(X,V(i,:),extrax,npoints);
            
        else
            
            [nx,ny,addedy]=tsLinearExtrapolation(X,V(i,:),extrax,npoints);
            
        end
        
        nvals(i,:)=interp1(nx(:)',ny(:)',Xq(:)');
        
    elseif sum(Xq>=min(X))==length(Xq) & sum(Xq<=max(X))==length(Xq)
        
        nvals(i,:)=interp1(X(:)',V(i,:),Xq(:)');
        
    end
    
end
