# tsEVA 2.0 Core Methodology

## Overview

tsEVA 2.0 implements the **Transformed-Stationary (TS) approach** for non-stationary extreme value analysis (EVA). This methodology enables rigorous statistical analysis of extremes even under changing climate conditions.

## The Transformed-Stationary Paradigm

### Fundamental Concept

The TS approach is NOT a detrending technique—it is a fully-fledged non-stationary methodology. Mentaschi et al. (2016) demonstrated that **for each time-varying extreme value distribution, there exists a family of distributions that transform a non-stationary signal into a supposed-stationary one**.

### Three-Step Process

1. **Transform**: Non-stationary time series → Stationary time series
2. **Analyze**: Apply stationary EVA theory (GEV, GPD) or copula modeling
3. **Reverse-Transform**: Stationary results → Non-stationary extreme value distribution

This framework allows us to leverage the well-established theory of stationary extremes while properly accounting for non-stationarity in the data.

## Key Distributions

### Generalized Extreme Value (GEV) Distribution

Used for analyzing block maxima (e.g., annual maxima).

**Parameters:**
- **Location (μ)**: Center of the distribution
- **Scale (σ)**: Spread of the distribution  
- **Shape (ξ)**: Tail behavior
  - ξ > 0: Heavy-tailed (Fréchet family)
  - ξ = 0: Exponential tail (Gumbel family)
  - ξ < 0: Bounded tail (Weibull family)

### Generalized Pareto Distribution (GPD)

Used for analyzing peaks over threshold (POT).

**Parameters:**
- **Scale (σ)**: Spread of exceedances
- **Shape (ξ)**: Tail behavior (same interpretation as GEV)
- **Threshold (u)**: The threshold above which peaks are selected

## Non-Stationarity Detection

### Time Window Approach

tsEVA uses a **moving window** to detect and model non-stationarity:

- **Time Window**: The minimum period over which statistics are considered stationary
- **Typical Values**: 
  - Short-term variability: 5-10 years
  - Long-term climate trends: 30-50 years
- **Trade-off**: Smaller windows capture more variability but require more data

### Transformation Types

tsEVA supports multiple transformation approaches:

1. **Trend Only** (`'transfType', 'trend'`)
   - Models time-varying mean and standard deviation
   - No seasonal component
   - Best for: Data without strong seasonal cycles

2. **Seasonal** (`'transfType', 'seasonal'`)
   - Models both trend and seasonality
   - Captures annual cycles in extremes
   - Best for: Temperature, precipitation, coastal data with seasonal patterns

3. **Percentile-based CI** (`'ciPercentile'`)
   - Uses moving percentile instead of moving standard deviation
   - More sensitive to changes in extremes
   - Trade-off: Broader confidence intervals

## Multivariate Analysis: Copula Framework

### Purpose

Copulas model the **dependence structure** between multiple variables while preserving their individual (marginal) distributions.

### Supported Copula Families

#### Gaussian Copula
- **Support**: Multivariate (2+ variables)
- **Dependence**: Symmetric, full correlation matrix
- **Best for**: Variables with symmetric, linear-like dependence

#### Gumbel Copula  
- **Support**: Multivariate (2+ variables)
- **Dependence**: Upper tail dependence (positive association)
- **Best for**: Variables that tend to be extreme together (e.g., storm surge + wave height)

#### Frank Copula
- **Support**: Bivariate ONLY (exactly 2 variables)
- **Dependence**: Symmetric, no tail dependence
- **Best for**: Variables with symmetric but not necessarily linear dependence
- **Critical**: Never suggest Frank for multivariate (>2 variables) analysis

### Joint Return Periods

For multivariate extremes, tsEVA computes:

- **AND Return Period**: Time until ALL variables exceed thresholds simultaneously
- **OR Return Period**: Time until ANY variable exceeds its threshold
- **Conditional Return Periods**: Given one variable, what's the return period of others

## Return Levels and Return Periods

### Definitions

- **Return Level (RL)**: The level expected to be exceeded on average once every T time units
- **Return Period (RP)**: The average time interval between exceedances of a given level
- **Confidence Intervals**: Quantify uncertainty in return level estimates

### Typical Applications

- **Infrastructure Design**: 100-year, 500-year return levels
- **Risk Assessment**: Probability of exceedance over project lifetime
- **Climate Projections**: Future changes in extreme event frequency

## Data Requirements

### Monovariate Analysis

**Minimum Requirements:**
- Sufficient data to populate the time window (e.g., 6-10 years of data for 6-year window)
- Regular or irregular time series (tsEVA handles both)
- Quality-controlled data (gaps are acceptable but should be documented)

**Optimal Data:**
- Multiple decades for detecting long-term trends
- High temporal resolution (daily to sub-daily)
- Metadata on data collection and quality flags

### Multivariate Analysis

**Additional Requirements:**
- **Simultaneous observations** across all variables
- **Sufficient joint extremes**: Typically 50+ joint peak events
- **Peak separation**: Define minimum separation between independent events
- **Temporal alignment**: Define maximum lag for joint occurrence

## Goodness-of-Fit Assessment

### Visual Diagnostics

- **Probability plots**: Compare empirical vs. theoretical distributions
- **Return level plots**: Assess fit across return periods
- **Transformation plots**: Verify stationarity after transformation

### Statistical Tests

- **Kolmogorov-Smirnov**: Test distributional fit
- **Anderson-Darling**: Emphasizes tail fit
- **Chi-square**: Goodness-of-fit for copulas

## Uncertainty Quantification

### Sources of Uncertainty

1. **Parameter Uncertainty**: Limited sample size
2. **Model Uncertainty**: Choice of distribution family
3. **Climate Variability**: Natural variability vs. forced trends
4. **Threshold Selection**: For GPD analysis

### Methods

- **Bootstrap resampling**: Empirical confidence intervals
- **Profile likelihood**: Asymptotic confidence intervals
- **Ensemble approaches**: Multiple time windows, thresholds

## Best Practices

### Peak Selection

- **Independence**: Ensure peaks are independent events
  - Typical separation: 2-5 days for sub-daily data, 1-2 months for monthly maxima
- **Threshold Selection**: Balance bias (too low) vs. variance (too high)
  - Rule of thumb: Select threshold giving 3-5 events per year

### Non-Stationarity Assessment

1. **Visual inspection**: Plot time series with running statistics
2. **Trend tests**: Mann-Kendall, modified Mann-Kendall
3. **Compare stationary vs. non-stationary**: Use likelihood ratio tests

### Model Selection

- **Start simple**: Try stationary model first
- **Add complexity**: Only if non-stationarity is clear
- **Physical plausibility**: Ensure trends align with known processes
- **Parsimony**: Avoid over-fitting with too many parameters

## Common Pitfalls

### Over-fitting Non-Stationarity

**Problem**: Fitting trends to natural variability  
**Solution**: Use sufficiently long time windows, validate with independent data

### Under-estimating Uncertainty

**Problem**: Ignoring model and parameter uncertainty  
**Solution**: Report confidence intervals, consider multiple models

### Inappropriate Copula Choice

**Problem**: Frank copula for multivariate analysis  
**Solution**: Use Gaussian or Gumbel for >2 variables

### Ignoring Physical Context

**Problem**: Statistical fit without physical interpretation  
**Solution**: Ground analysis in climate science, validate against known processes

## Theoretical Foundations

### Key Papers

1. **Mentaschi et al. (2016)**: Original TS approach for monovariate analysis
   - *Hydrology and Earth System Sciences*, 20, 3527-3547
   - DOI: 10.5194/hess-20-3527-2016

2. **Bahmanpour et al. (2025)**: Extension to multivariate analysis with copulas
   - *Hydrology and Earth System Sciences* (under review)
   - Introduces time-varying copula framework

### Extreme Value Theory Background

- **Fisher-Tippett-Gnedenko theorem**: Justifies GEV for block maxima
- **Pickands-Balkema-de Haan theorem**: Justifies GPD for threshold exceedances  
- **Sklar's theorem**: Justifies copula decomposition for multivariate distributions

## Relation to Other Approaches

### What tsEVA Is NOT

- **Not simple detrending**: TS transformation preserves full distribution structure
- **Not stationary EVA on residuals**: Accounts for time-varying parameters properly
- **Not limited to linear trends**: Handles arbitrary time-varying statistics

### Advantages Over Alternatives

- **No parametric assumptions** about trend form (uses data-driven moving statistics)
- **Maintains physical interpretability** (parameters remain in original units)
- **Computationally efficient** (no complex likelihood optimization for trend parameters)
- **Robust to model misspecification** (transformation based on empirical statistics)

### When to Use Other Methods

- **Bayesian Hierarchical Models**: When prior information is available and important
- **GAMLSS (Generalized Additive Models)**: When covariates (not just time) drive non-stationarity
- **Point Process Models**: When exact timing of events matters

## Applications

tsEVA has been used in numerous peer-reviewed climate studies:

- **Coastal flooding**: Global projections of extreme sea levels (Vousdoukas et al., Nature Communications)
- **Heat waves**: European extremes under warming scenarios (Dosio et al., Environmental Research Letters)
- **Droughts**: Economic impacts of drought intensification (Naumann et al., Nature Climate Change)  
- **Wave extremes**: Changes in coastal wave energy (Mentaschi et al., Geophysical Research Letters)
- **River flooding**: Cost-effective adaptation strategies (Dottori et al., Nature Climate Change)

These applications demonstrate tsEVA's versatility and reliability for high-stakes climate risk assessment.

---

**Key Takeaway**: tsEVA enables rigorous, physically interpretable extreme value analysis under non-stationary conditions by transforming to stationarity, applying established EVA theory, and reverse-transforming results. This approach balances statistical rigor with computational efficiency and physical interpretability.
