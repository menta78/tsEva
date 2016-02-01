function phandles = tsEvaPlotTransfToStat( timeStamps, statSeries, srsmean, stdDev, thirdMom, fourthMom, varargin )
  %stdThirdMom = third root of the third statistical momentum
  %stdFouthMom = fourth root of the fourth statistical momentum
  
  args.axisFontSize = 20;
  args.legendFontSize = 20;
  args.xtick = [];
  args.figPosition = [0, 0, 1450, 700] + 10;
  args.minYear = -7000;
  args.maxYear = 11000;
  args.dateformat = 'yyyy';
  args.legendLocation = 'northeast';
  args = tsEasyParseNamedArgs(varargin, args);

  minTS = datenum([args.minYear, 1, 1]);
  maxTS = datenum([args.maxYear, 1, 1]);
  statSeries = statSeries( (timeStamps >= minTS) & (timeStamps <= maxTS) );
  srsmean = srsmean( (timeStamps >= minTS) & (timeStamps <= maxTS) );
  stdDev = stdDev( (timeStamps >= minTS) & (timeStamps <= maxTS) );
  thirdMom = thirdMom( (timeStamps >= minTS) & (timeStamps <= maxTS) );
  fourthMom = fourthMom( (timeStamps >= minTS) & (timeStamps <= maxTS) );
  timeStamps = timeStamps( (timeStamps >= minTS) & (timeStamps <= maxTS) );
  minTS = min(timeStamps);
  maxTS = max(timeStamps);
  
  f = figure;
  f.Position = args.figPosition;

  psrs = plot(timeStamps, statSeries);
  hold on;
  pmean = plot(timeStamps, srsmean, '--', 'color', 'k', 'linewidth', 3);
  pStdDev = plot(timeStamps, stdDev, '--', 'color', [.5, 0, 0], 'linewidth', 3);
  %axe1 = gca;
  %axe1pos = axe1.Position;
  %axe2 = axes('Position', axe1pos, 'xaxislocation', 'top', 'yaxislocation', 'right', 'color', 'none');
  %pThirdMom = line(timeStamps, stdThirdMom, 'parent', axe2, 'color', [0, 0, .5], 'linewidth', 3);
  %pFourthMom = line(timeStamps, stdFourthMom, 'color', [0, .4, 0], 'linewidth', 3, 'parent', axe2);
  %axes(axe1);
  pThirdMom = plot(timeStamps, thirdMom, 'color', [0, 0, .5], 'linewidth', 3);
  pFourthMom = plot(timeStamps, fourthMom, 'color', [0, .4, 0], 'linewidth', 3);
  %thirdMomThrdRoot = nthroot(thirdMom, 3.);
  %fourthMomFourthRoot = nthroot(fourthMom, 4.);
  %pThirdMomThrdRoot = plot(timeStamps, thirdMomThrdRoot, '--', 'color', [0, 0, .5], 'linewidth', 3);
  %pFourthMomFourthRoot = plot(timeStamps, fourthMomFourthRoot, '--', 'color', [0, .4, 0], 'linewidth', 3);
  
  datetick('x', args.dateformat);
  xlim([minTS maxTS]);
%   pleg = legend([psrs pmean pStdDev pThirdMom pThirdMomThrdRoot pFourthMom pFourthMomFourthRoot], {'Normalized series', 'Mean', 'Std. dev.', 'Skewness', 'Skewness (3rd root)', 'Kurtosis', 'Kurtosis (4th root)'},...
%       'fontsize', args.legendFontSize, 'location', args.legendLocation);
  pleg = legend([psrs pmean pStdDev pThirdMom pFourthMom], {'Normalized series', 'Mean', 'Std. dev.', 'Skewness', 'Kurtosis'},...
      'fontsize', args.legendFontSize, 'location', args.legendLocation);
  set(gca, 'fontsize', args.axisFontSize);
  if ~isempty(args.xtick)
      set(gca, 'xtick', args.xtick);
      set(gca, 'xticklabel', datestr(args.xtick, args.dateformat));
  end
  
  grid on;
  set(gcf, 'paperpositionmode', 'auto');
  hold off;

  phandles = {f, psrs, pmean, pStdDev, pStdDev, pThirdMom, pFourthMom, pleg};
end

