function phandles = plotGPD3D( X, timeStamps, epsilon, sigma, threshold, varargin )
  args.xlabel = 'levels (m)';
  args.zlabel = 'pdf';
  args = tsEasyParseNamedArgs(varargin, args);

  f = figure;
  phandles{1} = f;
  f.Position = [0, 0, 1300, 700] + 10;
  
  L = length(timeStamps);
  npdf = min(100, L);
  navg = ceil(L/npdf);
  
  deltai = round(L/npdf);
  if length(epsilon) == 1
      epsilon0 = ones(npdf, 1)*epsilon;
  else
      epsilon_ = nan*ones(npdf*navg, 1);
      epsilon_(1:L) = epsilon(:);
      epsilonMtx = reshape(epsilon_, navg, []);
      epsilon0 = nanmean(epsilonMtx)';
  end
  
  sigma_ = nan*ones(npdf*navg, 1);
  sigma_(1:L) = sigma(:);
  sigmaMtx = reshape(sigma_, navg, []);
  sigma0 = nanmean(sigmaMtx)';
  
  threshold_ = nan*ones(npdf*navg, 1);
  threshold_(1:L) = threshold(:);
  thresholdMtx = reshape(threshold_, navg, []);
  threshold0 = nanmean(thresholdMtx)';
  
  tmstmps = timeStamps(1:deltai:end);
  tmstmps = tmstmps(1:length(threshold0));
  
  [~, epsilonMtx] = meshgrid(X, epsilon0);
  [~, sigmaMtx] = meshgrid(X, sigma0);
  [XMtx, thresholdMtx] = meshgrid(X, threshold0);
  
  gevvar = gppdf(XMtx, epsilonMtx, sigmaMtx, thresholdMtx);
  
  phandles{2} = surf(X, tmstmps, gevvar, 'linestyle', 'none');
  datetick('y');
  view(24.3, 48);

  xlabel(args.xlabel, 'fontsize', 22);
  zlabel(args.zlabel, 'fontsize', 22);
  
  set(gca, 'fontsize', 16);
  
  set(f, 'paperpositionmode', 'auto');
end

