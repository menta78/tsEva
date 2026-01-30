function copulaParam = tsCopulaFit(copulaFamily, uProb)

    if strcmpi(copulaFamily, 'gaussian')
        copulaParam = corr(uProb, 'type', 'Spearman'); % apparently works bettern than copulafit
    elseif strcmpi(copulaFamily, 'gumbel') || strcmpi(copulaFamily, 'clayton') || strcmpi(copulaFamily, 'frank')
        nSeries = size(uProb, 2);
        copulaParam = ones(nSeries);
        kendalT = corr(uProb, 'type', 'Kendall');
        if strcmpi(copulaFamily, 'gumbel')
            kendalT(triu(kendalT) < 0) = 0;

            % Replace negatives in the lower triangle
            kendalT(tril(kendalT) < 0) = 0;
        end
        for iSeries1 = 1:nSeries
            for iSeries2 = iSeries1+1:nSeries
                copulaParam(iSeries1, iSeries2) = copulaparam(copulaFamily,kendalT(iSeries1,iSeries2));
                copulaParam(iSeries2, iSeries1) = copulaParam(iSeries1, iSeries2);
            end
        end
        if strcmpi(copulaFamily, 'gumbel')
            copulaParam(triu(copulaParam) == Inf) = 1;

            % Replace negatives in the lower triangle
            copulaParam(tril(copulaParam) == Inf) = 1;
        end
    else
        error('copulaFamily not supported: ' + copulaFamily);
    end

end