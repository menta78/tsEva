function [tau, z, p, H] = Modified_MannKendall_test(t, X, alpha, alpha_ac)
    
    %% FUNCTION INPUTS AND OUTPUTS
    
    % INPUTS
    % t         - time vector corresponding to the timeseries.
    % X         - The timeseries for which kendall's tau value is to be evaluated.
    % alpha     - significance level above which kendall tau values are statistically significant. One-tailed test.
    % alpha_ac  - predetermined level for selecting only those autocorrelation lag values that are statistically significant. Two-tailed test.

    % OUTPUTS
    % tau   - The kendall rank correlation coefficient for the timeseries.
    % z     - z-score of the obtained tau value.
    % p     - p-value of the obtained tau value.
    % H     - Denotes whether to reject the null hypothesis or not. Null Hypothesis: There is no trend in the data. 0 -> retain the null hypothesis, 1 -> reject and increasing trend, -1 -> reject and decreasing trend, 2 -> if variance turns out to be negative.


    %% REFERENCES

    % 1) Kendall, M. G. (1938). A new measure of rank correlation. Biometrika, 30(1/2), 81-93.
    % 2) Kendall, M. G. (1948). Rank correlation methods.
    % 3) Hamed, K. H., & Rao, A. R. (1998). A modified Mann-Kendall trend test for autocorrelated data. Journal of hydrology, 204(1-4), 182-196.
    % 4) Hamed, K. H. (2008). Trend detection in hydrologic data: the Mann–Kendall trend test under the scaling hypothesis. Journal of hydrology, 349(3-4), 350-363.
    % 5) Yue, S., & Wang, C. (2004). The Mann-Kendall test modified by effective sample size to detect trend in serially correlated hydrological series. Water resources management, 18(3), 201-218.


    %% CONVERT DATA INTO ROW VECTOR
    
    % Convert input time and timeseries vectors to row vectors
    X = reshape(X, 1, length(X));
    t = reshape(t, 1, length(t));
    n = length(X);


    %% CALCULATE THE KENDALL TAU VALUE
    
    % Notation taken from: 3) Hamed, K. H., & Rao, A. R. (1998)
    S = 0;
    
    for i = 1: n-1
        for j = i+1: n
            S = S + sign(X(j) - X(i));
        end
    end
    
    % Calculate kendall's rank correlation coefficient - tau
    tau = S / (n * (n-1) / 2);


    %% CALCULATE VARIANCE OF KENDALL TAU UNCORRECTED FOR AUTOCORRELATION
    
    % Correction for tied ranks taken from:
    % 2) Kendall, M. G. (1948) edition 5, chapter 4, page 66
    % 4) Hamed, K. H. (2008)

    % Variance of kendall tau when the following conditions hold: 1) No autocorrelation, 2) No ties
    var_S_notie = n * (n - 1) * (2*n + 5) / 18;
    
    % Correcting variance in the case of tied ranks    
    tied_ranks = [];
    Y = sort(X);
    % Add dummy value to Y. Otherwise loop doesn't count tied ranks that exist at the end of Y.
    Y = [Y, Y(end) + 1];
    prev_counter = 1;
    current_counter = 1;
    for i_Y = 2: n+1
        if Y(i_Y) == Y(i_Y - 1)
            current_counter = current_counter + 1;
        else
            prev_counter = current_counter;
            current_counter = 1;
        end

        if current_counter == 1 && prev_counter ~= 1
            tied_ranks = [tied_ranks, prev_counter];
        end

    end

    var_S_tie_correction = sum(tied_ranks .* (tied_ranks - 1) .* (2*tied_ranks + 5) / 18);

    % Calculate variance corrected for ties but without considering autocorrelation
    var_S_noAC = var_S_notie - var_S_tie_correction;

    
    %% REMOVE SEN TREND ESTIMATE FROM THE DATA
    
    % Procedure and rationale for trend removal taken from:
    % 5) Yue, S., & Wang, C. (2004)
    % Wikipedia: Theil–Sen estimator, https://en.wikipedia.org/wiki/Theil-Sen_estimator
    % Basically, the presence of a trend leads to wrong estimation of the actual autocorrelation present. Therefore, the trend must first be removed before estimating the autocorrelation.
    m_list = zeros(1, (n * (n-1) / 2));
    b_list = [];
    
    element_counter = 1;
    for i = 1: n-1
        for j = i+1: n
            m_list(element_counter) = (X(j) - X(i)) / (t(j) - t(i));
            element_counter = element_counter + 1;
        end
    end
    
    m_sen = median(m_list);

    b_list = X - m_sen*t;
    b_sen = median(b_list);
    
    % Remove sen trend estimate from the data
    X = X - m_sen*t - b_sen;


    %% CALCULATE AUTOCORRELATION VALUES FOR STATISTICALLY SIGNIFICANT LAGS
    
    X_rank_order = tiedrank(X);
    z_ac = abs(norminv(alpha_ac / 2));  % norminv() is the inverse of the normcdf() function
    [acf, ~, acf_bounds] = autocorr(X_rank_order, NumLags = n - 1, NumSTD = z_ac);
    
    % Retain only those lags for which the autocorrelation value is statistically significant
    rho_lags = [];
    rho = [];

    for i = 2: n
        if acf(i) > acf_bounds(1) || acf(i) < acf_bounds(2)
            rho = [rho, acf(i)];
            rho_lags = [rho_lags, i-1];
        end
    end
    

    %% CALCULATE AUTOCORRELATION CORRECTED VARIANCE OF KENDALL TAU

    % Calculate variance correction factor
    const_factor = 2 / (n * (n-1) * (n-2));
    rho_factor_sum = 0;
    for i = 1: length(rho)
        rho_factor_sum = rho_factor_sum + (n - rho_lags(i)) * (n - rho_lags(i) - 1) * (n - rho_lags(i) - 2) * rho(i);
    end

    var_AC_correction_factor = 1 + const_factor * rho_factor_sum;

    % Calculate variance corrected for autocorrelation
    var_S = var_S_noAC * (var_AC_correction_factor);


    %% CHECK FOR STATISTICAL SIGNIFICANCE OF THE KENDALL TAU VALUE
    
    % Since the correction factor for the true variance is an approximation, in rare cases it may turn out to be negative. In that scenario abort the function and return H = 2 as an exception value.
    if var_S < 0
        z = 0;
        p = 0.5;
        H = 2;
        return
    end

    % Calculate z-score.
    % z-score = (value - mean) / std_dev. That is, how far is the value from mean as a multiple of the standard deviation.
    % Below the +1 and -1 are for "continuity correction".
    % Continuity correction taken from: 2) Kendall, M. G. (1948), edition 5, chapter 4, page 65
    if S > 0
        z = (S - 1) / sqrt(var_S);
    elseif S == 0
        z = 0;
    elseif S < 0
        z = (S + 1) / sqrt(var_S);
    end

    % Calculate p-value.
    % p-value is the probability of obtaining values further from the obtained value, on either the negative or positive tail.
    % p-value found here is for one-tailed test. Right tail values denote positive trend probabilities and left tail values represent negative trend probabilities
    if z >= 0
        p = 1 - normcdf(z);   % normcdf() returns the cdf value from the standard normal.
    elseif z < 0
        p = normcdf(z);
    end
    
    % Set H0 (Null Hypothesis) rejection value.
    % H = 0 -> retained. H = 1 -> rejected and increasing. H = -1 -> rejected and decreasing.
    if p <= alpha
        % That is there is a trend
        if z > 0
            % Positive trend
            H = 1;
        elseif z < 0
            % Negative trend
            H = -1;
        end
    else
        % That is there is no trend
        H = 0;
    end


end