function phandles = tsEvaPlotGPDImageSc( Y, timeStamps, epsilon, sigma, threshold, varargin )
  avgYearLength = 365.2425;
  nyears = (max(timeStamps) - min(timeStamps))/avgYearLength;
  nelmPerYear = length(timeStamps)/nyears;
  
  args.nPlottedTimesByYear = min(360, round(nelmPerYear));
  args.ylabel = 'levels (m)';
  args.zlabel = 'pdf';
  args.minYear = -7000;
  args.maxYear = 11000;
  args.dateFormat = 'yyyy';
  args.axisFontSize = 22;
  args.labelFontSize = 28;
  args.figPosition = [0, 0, 1450, 700] + 10;
  args.xtick = [];
  args = tsEasyParseNamedArgs(varargin, args);

  minTS = datenum([args.minYear, 1, 1]);
  maxTS = datenum([args.maxYear, 1, 1]);
  sigma = sigma( (timeStamps >= minTS) & (timeStamps <= maxTS) );
  threshold = threshold( (timeStamps >= minTS) & (timeStamps <= maxTS) );
  timeStamps = timeStamps( (timeStamps >= minTS) & (timeStamps <= maxTS) );
 
  f = figure;
  phandles{1} = f;
  f.Position = args.figPosition;
  
  L = length(timeStamps);
  minTS = timeStamps(1);
  maxTS = timeStamps(end);
  %npdf = (year(maxTS) - year(minTS) + 1)*args.nPlottedTimesByYear;
  npdf = ceil( ((maxTS - minTS)/avgYearLength)*args.nPlottedTimesByYear );
  navg = ceil(L/npdf);

  plotSLength = npdf*navg;
  timeStamps_plot = linspace(minTS, maxTS, plotSLength);
  
  if length(epsilon) == 1
      epsilon0 = ones(npdf, 1)*epsilon;
  else
      epsilon_ = nan*ones(npdf*navg, 1);
      epsilon_(1:L) = epsilon(:);
      epsilonMtx = reshape(epsilon_, navg, []);
      epsilon0 = nanmean(epsilonMtx)';
  end
  
  sigma_ = interp1(timeStamps, sigma, timeStamps_plot);
  sigmaMtx = reshape(sigma_, navg, []);
  sigma0 = nanmean(sigmaMtx)';
  
  threshold_ = interp1(timeStamps, threshold, timeStamps_plot);
  thresholdMtx = reshape(threshold_, navg, []);
  threshold0 = nanmean(thresholdMtx)';
  
  [~, epsilonMtx] = meshgrid(Y, epsilon0);
  [~, sigmaMtx] = meshgrid(Y, sigma0);
  [XMtx, thresholdMtx] = meshgrid(Y, threshold0);
  
  gevvar = gppdf(XMtx, epsilonMtx, sigmaMtx, thresholdMtx);
  
  phandles{2} = imagesc(timeStamps_plot, Y, gevvar');
  colormap(flipud(hot));
  set(gca,'YDir','normal');
  datetick('x', args.dateFormat);
  if ~isempty(args.xtick)
      set(gca, 'xtick', args.xtick);
      set(gca, 'xticklabel', datestr(args.xtick, args.dateFormat));
  end
  xlim([min(timeStamps_plot) max(timeStamps_plot)]);
  grid on;
  clb = colorbar;
  ylabel(clb, args.zlabel, 'fontsize', args.labelFontSize);

  ylabel(args.ylabel, 'fontsize', args.labelFontSize);
  
  set(gca, 'fontsize', args.axisFontSize);
  
  set(f, 'paperpositionmode', 'auto');
end

