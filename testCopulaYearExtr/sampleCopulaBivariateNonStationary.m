dt = load('testESL_closePts_bivariateNonStationary.mat');
timeAndSeries1 = dt.timeAndSeries1;
timeAndSeries2 = dt.timeAndSeries2;
clear dt;

timeWindow = 365.25*30; % 30 years
minPeakDistanceInDays = 3;
ciPercentile = 98.5;
returnPeriodsInYears = [1.5 2 3 4 5 7 10 15 20 30 50 70 100 150 250 350 500 700 1000 1500 2000];
years = (1970:2100)';

outTmStmp = datenum([years, ones(size(years)), ones(size(years))]);


disp('Performing the non-stationary extreme value analysis ...');
[nonStatEvaParams1, statTransfData1] = tsEvaNonStationary(timeAndSeries1, timeWindow, 'transfType', 'trendCiPercentile',... 
  'ciPercentile', ciPercentile, 'minPeakDistanceInDays', minPeakDistanceInDays);
[nonStatEvaParams1, statTransfData1] = tsEvaReduceOutputObjSize(nonStatEvaParams1, statTransfData1, outTmStmp, 'maxTimeStepDist', 365.25);
[rlev1, rlevErr1] = tsEvaComputeReturnLevelsGPDFromAnalysisObj(nonStatEvaParams1, returnPeriodsInYears);
yMax1 = tsEvaComputeAnnualMaxima(timeAndSeries1);


[nonStatEvaParams2, statTransfData2] = tsEvaNonStationary(timeAndSeries1, timeWindow, 'transfType', 'trendCiPercentile',... 
  'ciPercentile', ciPercentile, 'minPeakDistanceInDays', minPeakDistanceInDays);
[nonStatEvaParams2, statTransfData2] = tsEvaReduceOutputObjSize(nonStatEvaParams2, statTransfData2, outTmStmp, 'maxTimeStepDist', 365.25);
[rlev2, rlevErr2] = tsEvaComputeReturnLevelsGPDFromAnalysisObj(nonStatEvaParams2, returnPeriodsInYears);
yMax2 = tsEvaComputeAnnualMaxima(timeAndSeries2);


disp('Estimating the non-stationary joint distribution from the return levels ...');
retLev = cat(3, rlev1, rlev2);
yMax = [yMax1, yMax2];

[retLev, jpdf, cplParam] = tsCopulaYearExtrFit(returnPeriodsInYears, retLev, yMax, 'copulafamily', 'gaussian');

