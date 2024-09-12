% test case 01 of Bahmanpour et al., 2024
% this script is used to demonstrate the suitability of the code for analysing
% non-stationary joint distribution of swh and riverine discharge
% wave data is based on Mentaschi, et al., 2023; this dataset covers a period of 1950 - 2023
% with a temporal resolution of 3 hours; Mentaschi et al., 2023 data includes SSH,SWH,T01,MWD,MWL

%river data belongs to a point in the northern parts of France in the strait
% of Dover (in the English Channel) off the port of Port of Boulogne-sur-Mer
% and in front of Liane river mouth; 1.5368 E/50.7543 N:
% 6-hourly river discharge data covering period of 1950 - 2020 Liane river
%

% M. H. Bahmanpour, 2024

close all;
clc
clearvars;
addpath('../');
addpath('../m_map');
% copyfile('/Users/hadi/Downloads/remote-work-italy/tsEva_dvlp-copula_git9/Alois_collaboration/lorenzodata.mat',pwd)
% copyfile('/Users/hadi/Downloads/remote-work-italy/tsEva_dvlp-copula_git9/Alois_collaboration/alloisdata.mat',pwd)

load 4PtsPortugal    %1 - 2 - 4 were good 
load 3PtsPortugal2
load latlonwaveportugal

timeAndSeries1=[tt,timeAndSeries4Pts(:,1)];
timeAndSeries2=[tt,timeAndSeries4Pts(:,2)];
timeAndSeries3=[tt,timeAndSeries4Pts(:,3)];
timeAndSeries4=[tt,timeAndSeries4Pts(:,4)];
timeAndSeries5=[tt,timeAndSeries3Pts(:,1)];
timeAndSeries6=[tt,timeAndSeries3Pts(:,2)];
timeAndSeries7=[tt,timeAndSeries3Pts(:,3)];
pbb=nchoosek([1:7],3);

for ix=30%1:size(pbb,1)  %30 was selected
selectedSeries=pbb(ix,:);
% set some parameters

ciPercentile = [99,99,99];      %to be used in tsEvaNonStationary
potPercentiles=[{99},{99},{99}];   %o be used in tsEvaNonStationary; better be set to only one value but a range is also possible
timeWindow = 365*80; %timeWindow parameter needs to be set in days; if smaller than length of series, a stationary copula would be used
minPeakDistanceInDaysMonovarSampling=[0.5,0.5,0.5];    %minimum distance between monovariate peaks of each series
maxPeakDistanceInDaysMultivarSampling=[0.5]; %7%can either take one value or has to have a format and size matching size(nchoosek([1:numvar],2),1) where numvar is 2 in bivariate case, 3 in trivariate case, and so on
copulaFamily={'gaussian'};  %clayton, t and gaussian are possible choices

% perform copula analysis
ser01=eval(['timeAndSeries',num2str(selectedSeries(1)),'(:,2)']);
ser02=eval(['timeAndSeries',num2str(selectedSeries(2)),'(:,2)']);
ser03=eval(['timeAndSeries',num2str(selectedSeries(3)),'(:,2)']);

[copulaAnalysis] = tsCopulaExtremes(timeAndSeries1(:,1), ...
    [ser01,ser02,ser03], ...
    'minPeakDistanceInDaysMonovarSampling',minPeakDistanceInDaysMonovarSampling, ...
    'maxPeakDistanceInDaysMultivarSampling',maxPeakDistanceInDaysMultivarSampling, ...
    'copulaFamily',copulaFamily,...
    'transfType','trendlinear','timeWindow',timeWindow,...
    'ciPercentile',ciPercentile,'potPercentiles',potPercentiles,...
    'peakType','allExceedThreshold');


% perform Monte-Carlo analysis
[copulaAnalysis] = tsCopulaCompoundGPDMontecarlo(copulaAnalysis,...
    'nResample',1000,'timeIndex','middle');

%assess goodness-of-fit
[gofStatistics] = tsCopulaGOF(copulaAnalysis,'pValSn',0);

%plot joint return periods
% figHnd = tsCopulaPlot(copulaAnalysis,gofStatistics, ...
%     'ylbl', {'SWH_{1} (m)','SWH_{2} (m)','SWH_{3} (m)'});
char01=arrayfun(@(x) num2str(x),selectedSeries);
char02=[];
for ii=1:length(selectedSeries)

char02=[char02,append("Loc",char01(ii))];
end

append(char01,char02);
figHnd = tsCopulaPlotTrivariate(copulaAnalysis,gofStatistics, ...
    'ylbl', {'SWH (m)','SWH (m)','SWH (m)'},'locString',char02,'latlon',latlon(:,selectedSeries));

 char01x=arrayfun(@(x) num2str(x),selectedSeries);
saveas(gcf, [char01x,'.png'], 'png');

% 
end



 
  

   
