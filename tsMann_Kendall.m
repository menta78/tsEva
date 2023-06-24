%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MANN-KENDALL TEST ORIGINAL %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[H,p_value]=Mann_Kendall(V,alpha)
%%%%%%%%%%%%%%%%%
%%% Performs original Mann-Kendall test of the null hypothesis of trend
%%% absence in the vector V,  against the alternative of trend. 
%%% The result of the test is returned in H = 1 indicates
%%% a rejection of the null hypothesis at the alpha significance level. H = 0 indicates
%%% a failure to reject the null hypothesis at the alpha significance level.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUTS
%V = time series [vector]
%alpha =  significance level of the test [scalar]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% From Matlab Help %%%%%%%%%%%%%%%
%The significance level of a test is a threshold of probability a agreed
%to before the test is conducted. A typical value of alpha is 0.05. If the p-value of a test is less than alpha,
%the test rejects the null hypothesis. If the p-value is greater than alpha, there is insufficient evidence 
%to reject the null hypothesis. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OUTPUTS
%H = test result [1] Reject of Null Hypthesis [0] Insufficient evidence to reject the null hypothesis
%p_value = p-value of the test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% From Matlab Help %%%%%%%%%%%%%%%
%The p-value of a test is the probability, under the null hypothesis, of obtaining a value
%of the test statistic as extreme or more extreme than the value computed from
%the sample.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% References 
%Mann, H. B. (1945), Nonparametric tests against trend, Econometrica, 13, 
%245– 259.
%Kendall, M. G. (1975), Rank Correlation Methods, Griffin, London.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simone Fatichi -- simonef@dicea.unifi.it
%   Copyright 2009
%   $Date: 2009/10/03 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V=reshape(V,length(V),1); 
alpha = alpha/2; %
n=length(V); 
i=0; j=0; S=0; 
for i=1:n-1
   for j= i+1:n 
      S= S + sign(V(j)-V(i)); 
   end
end
VarS=(n*(n-1)*(2*n+5))/18;
StdS=sqrt(VarS); 
%%%% Note: ties are not considered 
if S >= 0
   Z=((S-1)/StdS)*(S~=0);
else
   Z=(S+1)/StdS;
end
p_value=2*(1-normcdf(abs(Z),0,1)); %% Two-tailed test 
pz=norminv(1-alpha,0,1); 
H=abs(Z)>pz; %% 
return 