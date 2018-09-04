load('testSpeiTmaxCopula.mat');

disp('testing stationary copula for yearly extremes of SPEI (droughts) and TMax');
disp(['at lon, lat = ' num2str(xtst) ', ' num2str(ytst)]);

[retLev, jdist, cplParam] = tsCopulaYearExtrFit(retPeriod, retLev, retLevError, yMax, 'copulafamily', 'gumbel');

nResample = 10000;
[resampleLevel, resampleProb] = tsCopulaYearExtrRnd(retPeriod, retLev, cplParam, nResample);

figHnd = tsCopulaYearExtrPlotSctrBivar(resampleLevel, yMax, 'xlbl', '-SPEI', 'ylbl', 'TMax');
title(['Milan, joint distribution ' newline ' of extreme SPEI and TMax']);

saveas(figHnd(1), 'milan_SpeiTmaxJointDist.png');
