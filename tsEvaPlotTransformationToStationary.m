function phandles = tsEvaPlotTransformationToStationary( timeStamps, statSeries, srsmean, stdDev, stdThirdMom, stdFourthMom, varargin )
  %stdThirdMom = third root of the third statistical momentum
  %stdFouthMom = fourth root of the fourth statistical momentum
  
  args.axisFontSize = 20;
  args.legendFontSize = 20;
  args = tsEasyParseNamedArgs(varargin, args);
  
  f = figure;
  f.Position = [0, 0, 1450, 700] + 10;

  psrs = plot(timeStamps, statSeries);
  pmean = plot(timeStamps, srsmean, 'k', 'linewidth', 3);
  pStdDev = plot(timeStamps, stdDev, 'color', [.5, 0, 0], 'linewidth', 3);
  pThirdMom = plot(timeStamps, stdThirdMom, 'color', [0, 0, .5], 'linewidth', 3);
  pFourthMom = plot(timeStamps, stdFourthMom, 'color', [0, .5, 0], 'linewidth', 3);
  
  legend([psrs pmean pStdDev pStdDev pThirdMom pFourthMom], {'Series', 'Mean', 'Std. dev.', 'Skewness (3rd root)', 'Kurtosis (4th root)'}, 'fontsize', args.legendFontSize);
  dateticks('x');
  set(gca, 'fontsize', args.axisFontSize);

  phandles = {f, psrs};
end

