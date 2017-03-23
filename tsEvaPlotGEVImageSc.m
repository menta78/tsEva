function phandles = tsEvaPlotGEVImageSc( Y, timeStamps, epsilon, sigma, mu, varargin )
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
  args.colormap = flipud(hot(64));
  args.plotColorbar = true;
  args.figPosition = [0, 0, 1450, 700] + 10;
  args.xtick = [];
  args.ax = [];
  args = tsEasyParseNamedArgs(varargin, args);

  minTS = datenum([args.minYear, 1, 1]);
  maxTS = datenum([args.maxYear, 1, 1]);
  sigma = sigma( (timeStamps >= minTS) & (timeStamps <= maxTS) );
  mu = mu( (timeStamps >= minTS) & (timeStamps <= maxTS) );
  timeStamps = timeStamps( (timeStamps >= minTS) & (timeStamps <= maxTS) );
  
  if isempty(args.ax)
    f = figure;
    phandles{1} = f;
    f.Position = args.figPosition;
  else
    axes(args.ax);
    phandles{1} = args.ax;
    f = [];
  end
  
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
      %epsilon0 = epsilon;
  else
      epsilon_ = nan*ones(npdf*navg, 1);
      epsilon_(1:L) = epsilon(:);
      epsilonMtx = reshape(epsilon_, navg, []);
      epsilon0 = nanmean(epsilonMtx)';
  end
  
  sigma_ = interp1(timeStamps, sigma, timeStamps_plot);
  sigmaMtx = reshape(sigma_, navg, []);
  if size(sigmaMtx, 1) > 1
    sigma0 = nanmean(sigmaMtx)';
  else
    sigma0 = sigmaMtx';
  end
  
  mu_ = interp1(timeStamps, mu, timeStamps_plot);
  muMtx = reshape(mu_, navg, []);
  if size(muMtx, 1) > 1
    mu0 = nanmean(muMtx)';
  else
    mu0 = muMtx';
  end
  
  [~, epsilonMtx] = meshgrid(Y, epsilon0);
  [~, sigmaMtx] = meshgrid(Y, sigma0);
  [XMtx, muMtx] = meshgrid(Y, mu0);
  
  gevvar = gevpdf(XMtx, epsilonMtx, sigmaMtx, muMtx);
  
  phandles{2} = imagesc(timeStamps_plot, Y, gevvar');
  colormap(args.colormap);
  set(gca,'YDir','normal');
  datetick('x', args.dateFormat);
  if ~isempty(args.xtick)
      set(gca, 'xtick', args.xtick);
      set(gca, 'xticklabel', datestr(args.xtick, args.dateFormat));
  end
  xlim([min(timeStamps_plot) max(timeStamps_plot)]);
  grid on;
  if args.plotColorbar
    clb = colorbar;
    ylabel(clb, args.zlabel, 'fontsize', args.labelFontSize);
  end

  ylabel(args.ylabel, 'fontsize', args.labelFontSize);
  
  set(gca, 'fontsize', args.axisFontSize);
  
  if ~isempty(f)
    set(f, 'paperpositionmode', 'auto');
  end
end

