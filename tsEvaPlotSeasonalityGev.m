function phandles = tsEvaPlotSeasonalityGev(extremesRange, referenceYear, timeStamps, epsilon, sigma, mu, monthlyMaxIndexes, series, trend, stddev, varargin)
%tsEvaPlotSeasonalityGev
%   plots a single year of data adding the series of monthly maxima
%   renormalized to the considered year and superimposed to the series.
%   I furtherly plots the time varying location and scale parameters.
    
    %args.colormap = flipud(gray);
    args.colormap = flipud(pink(64));
    args = tsEasyParseNamedArgs(varargin, args);

    minTS = datenum([referenceYear - 1, 12, 31]);
    maxTS = datenum([referenceYear + 1, 1, 2]);
    sigma = sigma( (timeStamps >= minTS) & (timeStamps <= maxTS) );
    mu = mu( (timeStamps >= minTS) & (timeStamps <= maxTS) );
    
    timeStampsRefYear = timeStamps( (timeStamps >= minTS) & (timeStamps <= maxTS) );

    plotVarArgIn = {'nPlottedTimesByYear', 360, 'dateFormat', 'mmm', 'colormap', args.colormap};
    plotVarArgIn = [plotVarArgIn varargin{:}];
    phandles = tsEvaPlotGEVImageSc(extremesRange, timeStampsRefYear, epsilon, sigma, mu,...
        plotVarArgIn{:});
    L = length(phandles);
    hold on;
    
    statSeries = (series - trend)./stddev;
    refYearTrend = mean(trend( (timeStamps >= minTS) & (timeStamps <= maxTS) ));
    refYearStdDev = mean(stddev( (timeStamps >= minTS) & (timeStamps <= maxTS) ));
    rescaledSeries = statSeries*refYearStdDev + refYearTrend;
    monthlyMax = rescaledSeries(monthlyMaxIndexes);
    monthlyMaxTimeStamp = timeStamps(monthlyMaxIndexes);
    
    mtsdatevec = datevec(monthlyMaxTimeStamp);
    mtsdatevec(:,1) = referenceYear;
    mts = datenum(mtsdatevec);
    monthlyMaxPlot = plot(mts, monthlyMax, 'o', 'color', 'b');
    phandles{L + 1} = monthlyMaxPlot;

    mupl = plot(timeStampsRefYear, mu, 'color', [.7, 0, 0], 'linewidth', 3);
    sigpl = plot(timeStampsRefYear, mu + sigma, 'color', [0, .5, 0], 'linewidth', 3);
    
    phandles{L + 2} = mupl;
    phandles{L + 3} = sigpl;
    
    legend([mupl, sigpl, monthlyMaxPlot], {'\mu', '\mu + \sigma', 'monthly max'}, 'fontsize', 22);
end

