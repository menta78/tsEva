function phandles = tsEvaPlotReturnLevelsGEV( epsilon, sigma, mu, epsilonStdErr, sigmaStdErr, muStdErr, varargin  )
%

args.minReturnPeriodYears = 5;
args.maxReturnPeriodYears = 1000;
args.confidenceAreaColor = [189 252	201]/255; % light green
args.confidenceBarColor = [34 139 34]/255; % dark green
args.returnLevelColor = 'k';
args.xlabel = 'return period (years)';
args.ylabel = 'return levels (m)';
args.ylim = [];
args.dtSampleYears = 1; % one year
args.ax = [];
args = tsEasyParseNamedArgs(varargin, args);
minReturnPeriodYears = args.minReturnPeriodYears;
maxReturnPeriodYears = args.maxReturnPeriodYears;

returnPeriodsInYears = logspace(log10(minReturnPeriodYears), log10(maxReturnPeriodYears));
retrunPeriodsInDts = returnPeriodsInYears/args.dtSampleYears;

[returnLevels, retrunLevelsErrs] = tsEvaComputeReturnLevelsGEV(epsilon, sigma, mu, epsilonStdErr, sigmaStdErr, muStdErr, retrunPeriodsInDts);

supRLCI = returnLevels + 2*retrunLevelsErrs;
infRLCI = returnLevels - 2*retrunLevelsErrs;
if ~isempty(args.ylim)
  maxRL = max(args.ylim);
  minRL = min(args.ylim);
else
  maxRL = max(supRLCI);
  minRL = min(infRLCI);
end

if isempty(args.ax)
  f = figure;
  phandles{1} = f;
  f.Position = [0, 0, 1300, 700] + 10;
else
  axes(args.ax);
  phandles{1} = args.ax;
  f = [];
end

h = area(returnPeriodsInYears, cat(1, infRLCI, supRLCI - infRLCI)', 'linestyle', 'none');
h(1).FaceColor = [1,1,1];
h(2).FaceColor = args.confidenceAreaColor;
hold on;
phandles{2} = h;

phandles{3} = plot(returnPeriodsInYears, returnLevels, 'color', args.returnLevelColor, 'linewidth', 3);
hold on;
phandles{4} = plot(returnPeriodsInYears, supRLCI, 'color', args.confidenceBarColor, 'linewidth', 2);
phandles{5} = plot(returnPeriodsInYears, infRLCI, 'color', args.confidenceBarColor, 'linewidth', 2);
set(gca, 'Xscale', 'log');

ylim([minRL maxRL]);
xlim([minReturnPeriodYears, maxReturnPeriodYears]);

grid on;
set(gca,'layer','top');

set(gca, 'fontsize', 20);
xlabel(args.xlabel, 'fontsize', 24);
ylabel(args.ylabel, 'fontsize', 24);

if ~isempty(f)
  set(f, 'paperpositionmode', 'auto');
end

end

