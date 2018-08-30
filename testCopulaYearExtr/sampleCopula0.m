load('testSpeiTmaxCopula.mat');

disp('testing stationary copula for yearly extremes of SPEI (droughts) and TMax');
disp(['at lon, lat = ' num2str(xtst) ', ' num2str(ytst)]);

[retLev, jdist] = tsCopulaYearExtr(retPeriod, retLev, retLevError, yMax, 'copulafamily', 'clayton');


