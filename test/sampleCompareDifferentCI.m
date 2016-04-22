function sampleCompareDifferentCI( varargin )

  args.timeWindow = 30*365.25; % 30 years
  args.timeAndSeries = [];
  args.minPeakDistanceInDays = 3;
  args.ylabel = 'Hs (m)';
  args = tsEasyParseNamedArgs(varargin, args);
  timeAndSeries = args.timeAndSeries;
  timeWindow = args.timeWindow;
  minPeakDistanceInDays = args.minPeakDistanceInDays;
  ylabel = args.ylabel;
  
  if isempty(timeAndSeries)
    load('timeAndSeries_waves_020_953E_-035_371N.mat');
    %load('timeAndSeries_waves_015_220E_055_509N.mat');
    %load('timeAndSeries_waves_023_688E_059_519N.mat');
  end
  
  disp('trying moving standard deviation ...');
  [nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, timeWindow, 'minPeakDistanceInDays', minPeakDistanceInDays);
  minext = (max(statTransfData.nonStatSeries) + 3*min(statTransfData.nonStatSeries))/4;
  maxext = max(statTransfData.nonStatSeries)*1.2;
  xext = minext:.01:maxext;
  disp('  plotting time varying GEV');
  tsEvaPlotGEVImageScFromAnalysisObj(xext, nonStatEvaParams, statTransfData, 'ylabel', ylabel);
  title('standard deviation');
  disp('  transformation diagnostic plot');
  tsEvaPlotTransfToStatFromAnalysisObj(nonStatEvaParams, statTransfData)
  title('standard deviation');
  
  disp('trying moving 98th percentile ...');
  [nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, timeWindow,...
    'transfType', 'trendCiPercentile','ciPercentile', 98, 'minPeakDistanceInDays', minPeakDistanceInDays);
  disp('  plotting time varying GEV');
  tsEvaPlotGEVImageScFromAnalysisObj(xext, nonStatEvaParams, statTransfData, 'ylabel', ylabel);
  title('98th percentile');
  disp('  transformation diagnostic plot');
  tsEvaPlotTransfToStatFromAnalysisObj(nonStatEvaParams, statTransfData)
  title('98th percentile');
  
  disp('trying moving 98.5th percentile ...');
  [nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, timeWindow,...
    'transfType', 'trendCiPercentile','ciPercentile', 98.5, 'minPeakDistanceInDays', minPeakDistanceInDays);
  disp('  plotting time varying GEV');
  tsEvaPlotGEVImageScFromAnalysisObj(xext, nonStatEvaParams, statTransfData, 'ylabel', ylabel);
  title('98.5th percentile');
  disp('  transformation diagnostic plot');
  tsEvaPlotTransfToStatFromAnalysisObj(nonStatEvaParams, statTransfData)
  title('98.5th percentile');
  
  disp('trying moving 99th percentile ...');
  [nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, timeWindow,...
    'transfType', 'trendCiPercentile','ciPercentile', 99, 'minPeakDistanceInDays', minPeakDistanceInDays);
  disp('  plotting time varying GEV');
  tsEvaPlotGEVImageScFromAnalysisObj(xext, nonStatEvaParams, statTransfData, 'ylabel', ylabel);
  title('99th percentile');
  disp('  transformation diagnostic plot');
  tsEvaPlotTransfToStatFromAnalysisObj(nonStatEvaParams, statTransfData)
  title('99th percentile');
  