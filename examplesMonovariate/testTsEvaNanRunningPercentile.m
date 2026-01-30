function testTsEvaNanRunningPercentile
% this function tests the function for computing of the running percentile

addpath('../');

load('timeAndSeriesHebrides.mat');
timeAndSeries = timeAndSeriesHebrides;
timeWindow = 365.25*6; % 6 years
percent = 80;

timeStamps = timeAndSeries(:,1);
series = timeAndSeries(:,2);

[ trendSeries, filledTimeStamps, filledSeries, nRunMn ] = tsEvaRunningMeanTrend( timeStamps, series, timeWindow);
plot(filledTimeStamps, filledSeries);
hold on;

percent = 80;
tic;
[ rnprcnt, err ] = tsEvaNanRunningPercentile( filledSeries, nRunMn, percent );
toc;
disp(['error = ' num2str(err/nanmean(rnprcnt)*100) ' %']);
plot(filledTimeStamps, rnprcnt);

percent = 90;
tic;
[ rnprcnt, err ] = tsEvaNanRunningPercentile( filledSeries, nRunMn, percent );
toc;
disp(['error = ' num2str(err/nanmean(rnprcnt)*100) ' %']);
plot(filledTimeStamps, rnprcnt);

percent = 95;
tic;
[ rnprcnt, err ] = tsEvaNanRunningPercentile( filledSeries, nRunMn, percent );
toc;
disp(['error = ' num2str(err/nanmean(rnprcnt)*100) ' %']);
plot(filledTimeStamps, rnprcnt);

percent = 98;
tic;
[ rnprcnt, err ] = tsEvaNanRunningPercentile( filledSeries, nRunMn, percent );
toc;
disp(['error = ' num2str(err/nanmean(rnprcnt)*100) ' %']);
plot(filledTimeStamps, rnprcnt);

percent = 99;
tic;
[ rnprcnt, err ] = tsEvaNanRunningPercentile( filledSeries, nRunMn, percent );
toc;
disp(['error = ' num2str(err/nanmean(rnprcnt)*100) ' %']);
plot(filledTimeStamps, rnprcnt);

datetick('x');
xlim([min(filledTimeStamps), max(filledTimeStamps)]);
