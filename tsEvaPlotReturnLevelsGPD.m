function [phandles, returnPeriods, returnLevels, retrunLevelsErrs] = tsEvaPlotReturnLevelsGPD( epsilon, sigma, threshold, epsilonStdErr, sigmaStdErr, thresholdStdErr, nPeaks, timeHorizonInYears, varargin  )
%

args.minReturnPeriodYears = 5;
args.maxReturnPeriodYears = 1000;
args.confidenceAreaColor = [189 252	201]/255; % light green
args.confidenceBarColor = [34 139 34]/255; % dark green
args.returnLevelColor = 'k';
args.xlabel = 'return period (years)';
args.ylabel = 'return levels (m)';
args.ylim = [];
args = tsEasyParseNamedArgs(varargin, args);
minReturnPeriodYears = args.minReturnPeriodYears;
maxReturnPeriodYears = args.maxReturnPeriodYears;

returnPeriods = logspace(log10(minReturnPeriodYears), log10(maxReturnPeriodYears));

[returnLevels, retrunLevelsErrs] = tsEvaComputeReturnLevelsGPD(epsilon, sigma, threshold, epsilonStdErr, sigmaStdErr, thresholdStdErr, nPeaks, timeHorizonInYears, returnPeriods);

supRLCI = returnLevels + 2*retrunLevelsErrs;
infRLCI = returnLevels - 2*retrunLevelsErrs;
if ~isempty(args.ylim)
  maxRL = max(args.ylim);
  minRL = min(args.ylim);
else
  maxRL = max(supRLCI);
  minRL = min(infRLCI);
end

f = figure;
phandles{1} = f;
f.Position = [0, 0, 1300, 700] + 10;
h = area(returnPeriods, cat(1, infRLCI, supRLCI - infRLCI)', 'linestyle', 'none');
h(1).FaceColor = [1,1,1];
h(2).FaceColor = args.confidenceAreaColor;
hold on;
phandles{2} = h;

phandles{3} = plot(returnPeriods, returnLevels, 'color', args.returnLevelColor, 'linewidth', 3);
hold on;
phandles{4} = plot(returnPeriods, supRLCI, 'color', args.confidenceBarColor, 'linewidth', 2);
phandles{5} = plot(returnPeriods, infRLCI, 'color', args.confidenceBarColor, 'linewidth', 2);
set(gca, 'Xscale', 'log');

ylim([minRL maxRL]);
xlim([minReturnPeriodYears, maxReturnPeriodYears]);

grid on;
set(gca,'layer','top');

set(gca, 'fontsize', 20);
xlabel(args.xlabel, 'fontsize', 24);
ylabel(args.ylabel, 'fontsize', 24);

set(f, 'paperpositionmode', 'auto');

end

