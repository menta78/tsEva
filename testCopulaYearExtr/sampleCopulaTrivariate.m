load('testESL_closePts_trivariate_Copula.mat');

disp('testing stat. copula for yearly extr. sea level, for 3 close locations in Portugal');
disp(['at lon, lat = ' num2str(xtst) ', ' num2str(ytst)]);

[retLev, jdist, cplParam] = tsCopulaYearExtrFit(retPeriod, retLev, retLevError, yMax, 'copulafamily', 'gaussian');

nResample = 10000;
[resampleLevel, resampleProb] = tsCopulaYearExtrRnd(retPeriod, retLev, cplParam, nResample);

ttl = ['Joint distribution ' newline ' of extr. sea lvl at 3 close locations in Portugal'];
figHnd = tsCopulaYearExtrPlotSctrTrivar(resampleLevel, yMax, 'xlbl', 'Location 1', 'ylbl', 'Location 2', 'zlbl', 'Location 3', 'title', ttl);

saveas(figHnd(1), 'testESL_closePts_trivariate_Copula.png');
