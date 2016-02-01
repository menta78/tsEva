function [inds,rinds]=tsSameValuesSegmentation(iii,varargin)
% 
% function [inds,rinds]=tsSameValuesSegmentation(iii,val)
% separates segments of same value val
% Inputs: 
% iii-the series
% optional variable-the value
% Outputs:
% inds-cell with the values of each continuous set
% rinds-cell with indexes of each continuous set
% 
% Michalis Vousdoukas 2009


if nargin==2
    
    val=varargin{1};
    
elseif nargin==1
    
    val=1;
    
end

inds=[];
rinds=[];

ll=length(iii);

difs=iii(2:ll)-iii(1:ll-1);

indsep=find(difs~=0);

l1=[1 ; indsep(:)+1];
l2=[indsep(:); ll];

cc=0;
for i=1:length(l1)
    
    if iii(l1(i))==val
       
        cc=cc+1;
        
        inds{cc,1}=iii(l1(i):l2(i));
    
        rinds{cc,1}=[l1(i):l2(i)];
        
    end
end