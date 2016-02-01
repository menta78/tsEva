function sampleGenerateSeriesEVAGraphs( timeAndSeries, extremeWaveRange, seriesId, seriesDescr, destFolder )
    timeWindow = 365*20; % about 20 years
    minTS = min(timeAndSeries(:,1));
    maxTS = max(timeAndSeries(:,1));
    axisFontSize = 20;
    axisFontSize3d = 16;
    labelFontSize = 28;
    titleFontSize = 30;
    
    destFolder = fullfile(destFolder, seriesId);
    system(['mkdir -p ' destFolder]);
    wr = linspace(min(extremeWaveRange), max(extremeWaveRange), 1501);
    slice1 = {'1985-1990' 1985 1990};
    slice2 = {'2085-2090' 2085 2090};
    
    disp('trend only statistics');
    [nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, timeWindow, 'transfType', 'trend');
    disp('  plotting the series');
    hndl = tsEvaPlotSeriesTrendStdDev(statTransfData.timeStamps, statTransfData.nonStatSeries, statTransfData.trendSeries, statTransfData.stdDevSeries,...
        'ylabel', 'Hs (m)', 'title', seriesDescr, 'titleFontSize', titleFontSize);
    disp('  saving the series plot');
    saveas(hndl{1}, fullfile(destFolder, 'seriesTrendOnly.png'));    
    disp('  plotting and saving the 3D GEV graph');
    hndl = tsEvaPlotGEV3DFromAnalysisObj(wr, nonStatEvaParams, statTransfData, 'xlabel', 'Hs (m)', 'axisfontsize', axisFontSize3d);
    title('GEV 3D', 'fontsize', titleFontSize);
    saveas(hndl{1}, fullfile(destFolder, 'GEV3DTrendOnly.png'), 'png');
    disp('  plotting and saving the 2D GEV graph');
    hndl = tsEvaPlotGEVImageScFromAnalysisObj(wr, nonStatEvaParams, statTransfData, 'ylabel', 'Hs (m)');
    title('GEV', 'fontsize', titleFontSize);
    saveas(hndl{1}, fullfile(destFolder, 'GEV2DTrendOnly.png'), 'png');
    disp('  plotting and saving the 2D GPD graph');
    hndl = tsEvaPlotGPDImageScFromAnalysisObj(wr, nonStatEvaParams, statTransfData, 'ylabel', 'Hs (m)');
    title('GPD', 'fontsize', titleFontSize);
    saveas(hndl{1}, fullfile(destFolder, 'GPD2DTrendOnly.png'), 'png');
    
    disp('plotting and saving stationary series');
    hndl = plot(statTransfData.timeStamps, statTransfData.stationarySeries);
    datetick('x', 'yyyy');
    xlim([minTS maxTS]);
    title('Stationary series (with seasonality)', 'fontsize', titleFontSize);
    set(gca, 'fontsize', axisFontSize);
    grid on;
    saveas(hndl, fullfile(destFolder, 'statSeriesTrendOnly.png'), 'png');
    
    disp('seasonal statistics');
    [nonStatEvaParams, statTransfData] = tsEvaNonStationary(timeAndSeries, timeWindow, 'transfType', 'seasonal');
    
    disp('plotting and saving stationary series');
    hndl = plot(statTransfData.timeStamps, statTransfData.stationarySeries);
    datetick('x', 'yyyy');
    xlim([minTS maxTS]);
    title('Stationary series', 'fontsize', titleFontSize);
    set(gca, 'fontsize', axisFontSize);
    grid on;
    saveas(hndl, fullfile(destFolder, 'statSeriesSeasonal.png'), 'png');
    
    disp('  plotting slice 1');
    slice = slice1;
    slid = slice{1};
    disp('    plotting the series');
    hndl = tsEvaPlotSeriesTrendStdDev(statTransfData.timeStamps, statTransfData.nonStatSeries, statTransfData.trendSeries, statTransfData.stdDevSeries,...
        'ylabel', 'Hs (m)', 'title', slid, 'minyear', slice{2}, 'maxyear', slice{3});
    disp('    saving the series plot');
    saveas(hndl{1}, fullfile(destFolder, ['seriesSeasonal_' slid '.png']));    
    disp('    plotting and saving the 3D GEV graph');
    hndl = tsEvaPlotGEV3DFromAnalysisObj(wr, nonStatEvaParams, statTransfData, 'xlabel', 'Hs (m)', 'minyear', slice{2}, 'maxyear', slice{3}, 'axisfontsize', axisFontSize3d);
    title(['GEV 3D, ' slid], 'fontsize', titleFontSize);
    saveas(hndl{1}, fullfile(destFolder, ['GEV3DSeasonal_' slid '.png']), 'png');
    disp('    plotting and saving the 2D GEV graph');
    hndl = tsEvaPlotGEVImageScFromAnalysisObj(wr, nonStatEvaParams, statTransfData, 'ylabel', 'Hs (m)', 'minyear', slice{2}, 'maxyear', slice{3});
    title(['GEV ' slid], 'fontsize', titleFontSize);
    saveas(hndl{1}, fullfile(destFolder, ['GEV2DSeasonal_' slid '.png']), 'png');
    disp('    plotting and saving the 2D GPD graph');
    hndl = tsEvaPlotGPDImageScFromAnalysisObj(wr, nonStatEvaParams, statTransfData, 'ylabel', 'Hs (m)', 'minyear', slice{2}, 'maxyear', slice{3});
    title(['GPD ' slid], 'fontsize', titleFontSize);
    saveas(hndl{1}, fullfile(destFolder, ['GPD2DSeasonal_' slid '.png']), 'png');
    
    disp('  plotting slice 2');
    slice = slice2;
    slid = slice{1};
    disp('    plotting the series');
    hndl = tsEvaPlotSeriesTrendStdDev(statTransfData.timeStamps, statTransfData.nonStatSeries, statTransfData.trendSeries, statTransfData.stdDevSeries,...
        'ylabel', 'Hs (m)', 'title', slid, 'minyear', slice{2}, 'maxyear', slice{3});
    disp('    saving the series plot');
    saveas(hndl{1}, fullfile(destFolder, ['seriesSeasonal_' slid '.png']));    
    disp('    plotting and saving the 3D GEV graph');
    hndl = tsEvaPlotGEV3DFromAnalysisObj(wr, nonStatEvaParams, statTransfData, 'xlabel', 'Hs (m)', 'minyear', slice{2}, 'maxyear', slice{3}, 'axisfontsize', axisFontSize3d);
    title(['GEV 3D, ' slid], 'fontsize', titleFontSize);
    saveas(hndl{1}, fullfile(destFolder, ['GEV3DSeasonal_' slid '.png']), 'png');
    disp('    plotting and saving the 2D GEV graph');
    hndl = tsEvaPlotGEVImageScFromAnalysisObj(wr, nonStatEvaParams, statTransfData, 'ylabel', 'Hs (m)', 'minyear', slice{2}, 'maxyear', slice{3});
    title(['GEV ' slid], 'fontsize', titleFontSize);
    saveas(hndl{1}, fullfile(destFolder, ['GEV2DSeasonal_' slid '.png']), 'png');
    disp('    plotting and saving the 2D GPD graph');
    hndl = tsEvaPlotGPDImageScFromAnalysisObj(wr, nonStatEvaParams, statTransfData, 'ylabel', 'Hs (m)', 'minyear', slice{2}, 'maxyear', slice{3});
    title(['GPD ' slid], 'fontsize', titleFontSize);
    saveas(hndl{1}, fullfile(destFolder, ['GPD2DSeasonal_' slid '.png']), 'png');
end

