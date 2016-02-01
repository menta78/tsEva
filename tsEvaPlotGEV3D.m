function phandles = plotGEV3D( X, timeStamps, epsilon, sigma, mu, varargin )
  avgYearLength = 365.2425;
  args.nPlottedTimesByYear = 180;
  args.xlabel = 'levels (m)';
  args.zlabel = 'pdf';
  args.minYear = -7000;
  args.maxYear = 11000;
  args.dateFormat = 'yyyy';
  args.axisFontSize = 22;
  args.labelFontSize = 28;
  args = tsEasyParseNamedArgs(varargin, args);

  minTS = datenum([args.minYear, 1, 1]);
  maxTS = datenum([args.maxYear, 1, 1]);
  sigma = sigma( (timeStamps >= minTS) & (timeStamps <= maxTS) );
  mu = mu( (timeStamps >= minTS) & (timeStamps <= maxTS) );
  timeStamps = timeStamps( (timeStamps >= minTS) & (timeStamps <= maxTS) );
  
  f = figure;
  phandles{1} = f;
  f.Position = [0, 0, 1300, 700] + 10;
  
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
  
  mu_ = interp1(timeStamps, mu, timeStamps_plot);
  muMtx = reshape(mu_, navg, []);
  mu0 = nanmean(muMtx)';
  
  timeStamps_plot = linspace(min(timeStamps), max(timeStamps), length(mu0));
  
  [~, epsilonMtx] = meshgrid(X, epsilon0);
  [~, sigmaMtx] = meshgrid(X, sigma0);
  [XMtx, muMtx] = meshgrid(X, mu0);
  
  gevvar = gevpdf(XMtx, epsilonMtx, sigmaMtx, muMtx);
  
  phandles{2} = surf(X, timeStamps_plot, gevvar, 'linestyle', 'none');
  datetick('y', args.dateFormat);
  view(24.3, 48);

  xlabel(args.xlabel, 'fontsize', args.labelFontSize);
  zlabel(args.zlabel, 'fontsize', args.labelFontSize);
  
  set(gca, 'fontsize', args.axisFontSize);
  
  set(f, 'paperpositionmode', 'auto');
end

