load('testSpeiTmaxCopula.mat');

disp('testing stationary copula for yearly extremes of SPEI (droughts) and TMax');
disp(['at lon, lat = ' num2str(xtst) ', ' num2str(ytst)]);

[retLev, jdist, cplParam] = tsCopulaYearExtrFit(retPeriod, retLev, retLevError, yMax, 'copulafamily', 'gumbel');

nResample = 1000;
[resampleLevel, resampleProb] = tsCopulaYearExtrRnd(retPeriod, retLev, cplParam, nResample);

% figure('position', [100, 100, 800, 800]);
% rsmplSctr = scatterhist(resampleLevel(:,1), resampleLevel(:,2), 'direction', 'out');
% rsmplSctrObj = findall(rsmplSctr(1), 'type', 'line');
% rsmplDist1Ax = rsmplSctr(2);
% rsmplDist2Ax = rsmplSctr(3);
% hold on;
% ymaxSctr = scatter(yMax(:,1), yMax(:,2), 'markerfacecolor', 'r');
% xlabel('-SPEI');
% ylabel('TMax');
% grid on;
% ax = gca;
% ax.FontSize = 15;
% lgnd = legend([ymaxSctr, rsmplSctrObj], {'Yearly Maxima', ['Joint Distribution' newline 'Montecarlo']}, 'fontsize', 15);
% lgndPos = lgnd.Position;
% newHigh = lgndPos(4)*2;
% lgnd.Position = [rsmplDist2Ax.Position(1), rsmplDist1Ax.Position(2) + rsmplDist1Ax.Position(4) - newHigh, lgndPos(3), newHigh];
% lgnd.PlotBoxAspectRatio

figHnd = tsCopulaYearExtrPlotSctrBivar(resampleLevel, yMax, 'xlbl', '-SPEI', 'ylbl', 'TMax');
title(['Milan, joint distribution ' newline ' of extreme SPEI and TMax']);

saveas(figHnd(1), 'milan_SpeiTmaxJointDist.png');
