addpath('../');

load('testSpeiTmaxCopula.mat');

disp('testing stationary copula for yearly extremes of SPEI (droughts) and TMax');
disp(['at lon, lat = ' num2str(xtst) ', ' num2str(ytst)]);

copulaFamily = 'gumbel';

kendallTauInput = corr(yMax, 'type', 'kendall');
spearmanCorrInput = corr(yMax, 'type', 'spearman');

[retLev, cplParam] = tsCopulaYearExtrFit(retPeriod, retLev, yMax, 'copulafamily', copulaFamily);

nResample = 100000;
[resampleLevel, resampleProb] = tsCopulaYearExtrRnd(retPeriod, retLev, cplParam, nResample);
%skipping the computation of the kendall tau for the output, as it is very slow
% kendallTauOutput = corr(resampleLevel, 'type', 'kendall');
spearmanCorrOutput = corr(resampleLevel, 'type', 'spearman');

figHnd = tsCopulaYearExtrPlotSctrBivar(resampleLevel, yMax, 'xlbl', '-SPEI', 'ylbl', 'TMax');
title(['Milan, joint distribution ' newline ' of extreme SPEI and TMax']);

saveas(figHnd(1), 'milan_SpeiTmaxJointDist.png');

[jpdf, jcdf] = tsCopulaYearExtrDistribution(retPeriod, cplParam, 'computeCdf', true);



clear all;

copulaFamily = 'gumbel';

load('testESL_closePts_bivariate_Copula.mat');

disp('testing stat. copula for yearly extr. sea level, for 2 close locations in Portugal');
disp(['at lon, lat = ' num2str(xtst) ', ' num2str(ytst)]);

kendallTauInput = corr(yMax, 'type', 'kendall');
spearmanCorrInput = corr(yMax, 'type', 'spearman');

[retLev, cplParam] = tsCopulaYearExtrFit(retPeriod, retLev, yMax, 'copulafamily', copulaFamily);

nResample = 10000;
[resampleLevel, resampleProb] = tsCopulaYearExtrRnd(retPeriod, retLev, cplParam, nResample);
%skipping the computation of the kendall tau for the output, as it is very slow
% kendallTauOutput = corr(resampleLevel, 'type', 'kendall');
spearmanCorrOutput = corr(resampleLevel, 'type', 'spearman');

figHnd = tsCopulaYearExtrPlotSctrBivar(resampleLevel, yMax, 'xlbl', 'Location 1', 'ylbl', 'Location 2');
title(['joint distribution ' newline ' of extr. sea lvl at 2 close locations in Portugal']);

saveas(figHnd(1), 'testESL_closePts_bivariate_Copula.png');
