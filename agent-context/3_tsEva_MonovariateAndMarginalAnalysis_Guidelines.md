# tsEVA 2.0 — Monovariate and Marginal Extreme Value Analysis Guidelines

## Scope and purpose of these guidelines

This document defines **methodological guidelines** for monovariate extreme value analysis (EVA) and marginal modelling in **tsEVA 2.0**.

It is intentionally **not limited to a single parameter or function**. Instead, it provides a structured way to document, justify, and apply key modelling choices that influence:

* the definition of extremes,
* the number of events entering tail models,
* the stability of fitted distributions,
* and the interpretability of results used in **risk assessment, infrastructure design, and climate adaptation**.

Within this document, individual sections address **specific tsEVA arguments and design choices** (e.g. POT threshold selection, declustering, event rates). Each section should be read as one component of a broader, coherent modelling framework.

---

## 0. Sampling of the extremes

tsEVA executes the sampling of the extremes on its own, so you don't need to pre-extract the annual (or monthly) maxima or implement POT manually.

If a user explicitly asks for seasonal / month-specific analysis, then you can suggest using the seasonal transformations (`seasonal`, `seasonalCiPercentile`). Otherwise ***NEVER*** suggest seasonal analysis on your own.


## 1. Peak-Over-Threshold (POT) modelling in tsEVA

Peak-Over-Threshold (POT) modelling is a central component of monovariate and marginal extreme value analysis in tsEVA. POT is used to define and extract extreme events from continuous time series and to fit tail distributions (typically GPD) in both stationary and non-stationary contexts.

In tsEVA, POT is not treated as a purely statistical optimisation problem, but as an **operational definition of extremes** that must remain physically interpretable, reproducible, and consistent across applications.

The following subsections describe how POT-related arguments are designed, interpreted, and applied in tsEVA.

---

### 1.1 Conceptual role of POT in tsEVA

In tsEVA, POT is not a purely statistical device but an **operational representation of extreme events**.

A POT threshold defines:

* which observations are considered *extreme events*;
* how many independent events enter the tail model;
* the effective balance between bias (threshold too low) and variance (threshold too high).

For this reason, tsEVA treats POT selection as a **controlled design choice**, not a black-box optimisation.

---

### 1.2 The `potPercentiles` parameter

#### Definition

`potPercentiles` specifies one or more **percentiles of the data distribution** used to define candidate POT thresholds.

Examples:

```matlab
potPercentiles = 99;              % single fixed threshold
potPercentiles = 97:0.5:99.5;     % scan over multiple candidate thresholds
```

Percentiles are always applied to the **original (non-transformed) series**.

---

### 1.3 POT threshold selection modes

#### 1.3.1 Multiple-percentile mode (adaptive selection)

When `potPercentiles` is a **vector**, tsEVA:

1. Computes a POT threshold for each percentile.
2. Extracts **declustered peaks** above each threshold.
3. Estimates the **mean number of events per year** for each candidate.
4. Selects the *winning percentile* according to a target event rate.

This mode is designed for **robust, adaptive threshold selection**.

Typical usage (coastal and marine applications):

```matlab
potPercentiles = 97:0.5:99.5;
meanEventsPerYear = 5;
```

---

#### 1.3.2 Single-percentile mode (fixed threshold)

When `potPercentiles` is a **scalar**, tsEVA:

* uses that percentile **directly**;
* performs **no automatic threshold selection**;
* treats the threshold as an explicit modelling assumption.

This mode is recommended when:

* strict reproducibility is required (e.g. continental or global analyses);
* a policy or literature-driven threshold must be enforced;
* inter-comparison across many sites is more important than local optimisation.

Example:

```matlab
potPercentiles = 99;
```

---

### 1.4 How the winning percentile is selected

When multiple percentiles are provided, tsEVA selects the percentile that yields a **mean number of independent POT events per year closest to** the user-defined target:

```matlab
meanEventsPerYear
```

#### Selection rule

1. For each candidate percentile `p`:

   * extract declustered peaks;
   * compute mean events/year.
2. Compute the absolute difference from the target event rate.
3. Select the percentile that **minimises this difference**.
4. If multiple candidates are equally close:

   * select the **highest percentile** (i.e. highest threshold, fewer events).

This rule is deterministic, transparent, and insensitive to the spacing of the percentile grid.

Candidate thresholds that yield too few total peaks are excluded to prevent unstable tail fits.

---

### 1.5 Uncertainty associated with the POT threshold

In tsEVA, POT threshold uncertainty reflects the **resolution of the percentile scan**, not statistical estimation error.

* It represents how much the threshold would change if a nearby percentile were selected.
* It is estimated locally around the selected percentile.

This uncertainty is propagated consistently into return-level uncertainty estimates.

---

### 1.6 Consistency between `ciPercentile` and `potPercentiles`

Good practice in tsEVA is to ensure **conceptual alignment** between the percentile used to characterise non-stationary amplitude (`ciPercentile`) and the percentile(s) used to define extremes (`potPercentiles`).

These two parameters act at different stages of the methodology, but they describe *related parts of the distribution*:

* `ciPercentile` controls how variability of the **upper tail** is tracked through time during the transformed–stationary step.
* `potPercentiles` control which observations are treated as **extreme events** in the tail model.

For this reason, they should be chosen so that they are **compatible in magnitude and intent**.

#### Single-percentile POT analyses

When a **single fixed POT percentile** is used (e.g. for large-scale or reproducible analyses), best practice is:

```matlab
ciPercentile  = potPercentiles;
```

This ensures that:

* the non-stationary normalization targets the *same tail level* that defines extremes;
* no artificial mismatch is introduced between transformation and event selection;
* results remain interpretable and reproducible across many time series.

#### Multiple-percentile POT scans

When `potPercentiles` is a **vector**, `ciPercentile` should:

* lie within the same percentile range, or
* be representative of the upper-tail variability relevant to the scan.

Typical, a tsEVA-consistent choice for coastal hazard is:

```matlab
ciPercentile  = 99;
potPercentiles = 97:0.5:99.5;
```
but the choice of these values is hazard-dependent.

In this configuration, for coastal hazard:

* the transformation tracks changes in extreme variability;
* POT selection remains adaptive but coherent with the normalization step.

Using a `ciPercentile` far below the POT range is discouraged, as it may:

* under-represent changes in extreme variability;
* leave residual non-stationarity in the tail;
* or distort the interpretation of fitted extreme-value parameters.

---

### 1.6 Changing POT behaviour

Users can fully control POT behaviour via:

* `potPercentiles` — range or fixed value
* `meanEventsPerYear` — target event rate
* `minPeakDistanceInDays` — declustering control

All parameters are passed transparently through tsEVA functions and can be modified without altering the codebase.

---

### 1.7 Recommended practices

#### Coastal and marine water levels / waves

* Use **multiple-percentile mode** during method development.
* Typical settings:

  ```matlab
  potPercentiles = 97:0.5:99.5;
  meanEventsPerYear = 5;
  ```
* Once validated, a **single fixed percentile** may be adopted for large-scale applications.

#### Large-scale or operational analyses

* Use **single-percentile mode** for reproducibility.
* Document the chosen percentile explicitly.

If one needs to choose a percentile beforehand, a commonly used reference for coastal applications is:

> Wahl et al. (2015). *Increasing risk of compound flooding from storm surge and rainfall for major US cities*. Nature Climate Change.

which states that a good choice is 99.

## 2. Non-stationary normalization and the `ciPercentile` parameter

Non-stationary extreme value analysis in tsEVA relies on the **Transformed–Stationary (TS) approach**, in which non-stationarity is first isolated and removed through a transformation, and extreme value theory is then applied to the transformed (approximately stationary) series.

Within this framework, the parameter `ciPercentile` plays a central role in how tsEVA estimates the **time-varying amplitude** of the process.

---

### 2.1 What `ciPercentile` is

`ciPercentile` defines the **percentile of the local distribution** used to estimate the amplitude (or scale) of the process within a moving time window during the transformation step.

In practical terms, for each time window tsEVA:

* estimates a slowly varying mean or trend component;
* estimates a slowly varying amplitude component;
* normalizes the series by this amplitude before applying EVA.

When `ciPercentile` is used, the amplitude is derived from a **high percentile of the data**, rather than from variance or standard deviation.

This mechanism is used consistently for:

* GPD-based POT analyses;
* GEV-based block maxima analyses;
* both stationary and non-stationary marginal modelling.

---

### 2.2 Why tsEVA uses a percentile instead of standard deviation

Using the standard deviation to characterize variability implicitly assumes:

* symmetry of the distribution;
* finite second moments;
* weak sensitivity to extremes.

These assumptions are often violated in environmental and climate time series, especially when the focus is on extremes.

By contrast, using a high percentile:

* is robust to skewness and heavy tails;
* directly reflects changes in the **upper tail**, which is most relevant for extremes;
* remains meaningful even when variance is unstable or poorly defined.

For this reason, percentile-based amplitude estimation is generally preferred in tsEVA when analysing extreme-value behaviour.

---

### 2.3 Relationship to the Transformed–Stationary philosophy

The goal of the TS approach is to ensure that **non-stationarity is captured by the transformation**, not absorbed into the extreme value parameters.

A well-chosen `ciPercentile`:

* removes long-term changes in variability from the transformed series;
* leads to more stable GPD or GEV parameters through time;
* improves the interpretability of non-stationary return levels.

An inappropriate choice may leave residual non-stationarity in the tail or over-normalize genuine changes.

---

### 2.4 How to choose `ciPercentile`

There is no universal value for `ciPercentile`. It should be chosen based on the **process under study** and the **type of extremes of interest**.

General guidance is as follows:

#### Coastal water levels and waves

* Use high percentiles (typically `99`).
* Rationale: variability changes are driven by storminess and extreme events rather than by the bulk of the distribution.

Example:

```matlab
ciPercentile = 99;
```

#### Large-scale or reproducible analyses

* Align `ciPercentile` with the fixed POT threshold when a single `potPercentiles` value is used.

Example:

```matlab
ciPercentile   = 99;
potPercentiles = 99;
```

This ensures coherence between non-stationary normalization and extreme-event definition.

#### Other environmental variables

* For indices or variables with bounded or sparse extremes (e.g. drought indices), lower percentiles (e.g. `80–90`) may be more appropriate.
* The choice should reflect where meaningful variability for extremes resides in the distribution.

---

### 2.5 Recommended practice

* Treat `ciPercentile` as a **modelling assumption**, not a tuning parameter.
* Choose it deliberately and document the rationale.
* Prefer percentile-based normalization over standard deviation when analysing extremes.
* Ensure conceptual consistency between `ciPercentile` and other tail-related parameters (e.g. `potPercentiles`).

---





## 3 Analysis output inspection

The output of tsEvaNonStationary consists of 2 structures: nonStatEvaParams, statTransfData (or similar name). We'll explain here how to handle them.

### 3.1 statTransfData
This is the 2nd output of tsEvaNonStationary.

To inspect/analyze statTransfData (aka stationaryTransformData), start with a structural audit: list fields, check array sizes, data types, and quantify missing values (NaN counts). Treat timeStamps as the only time axis saved: it is reconstructed by the algorithm and is the reference grid for all stored series (even if some components originate from running statistics). The core output is stationarySeries (post-transform) and it is the primary validity diagnostic: assess whether it looks “supposed-stationary” by checking for residual trend, residual seasonality, changing spread, and distributional drift. trendSeries represents the time-varying central tendency removed by the transform. stdDevSeries is a time-varying scale/tail-size proxy: historically it came from running standard deviation (hence the name), but under ciPercentile-style transforms it instead reflects the distance between the mean/trend and high percentiles (i.e., a robust tail amplitude measure). If seasonal analysis was applied, *NonSeasonal* fields store the long-term (non-seasonal) change separately from the seasonal component. Use statSer3Mom and statSer4Mom (3rd/4th moments of stationarySeries) to diagnose changing skewness/kurtosis after transformation. runningStatMultiplicity documents the running-stat support (windowing/multiplicity). Finally, interpret \*Error fields cautiously: some uncertainties are time-dependent vectors, others may be scalar constants depending on the estimation step.

When the trendlinear transformation is used, statTransfData may also include summary-change diagnostics computed on the annual series. These can include p-values assessing change (e.g., pValueChange\* such as annual/other aggregations) and percent change metrics over the analyzed period, notably the percentage change in the annual ciPercentile-based scale/tail proxy (e.g., percentChangePercentile) and a corresponding trend metric (e.g., percentChangeTrend). These fields are not always present: they are only computed/stored when trendlinear is used, so their absence should not be treated as an error—just an indication that a different transformation was applied.
In particular:
pValueChange: p-value of the annual cipercentile
pValueChangeStat: p-value of the high percentiles of the transformed series (should be non significant)
pValueChangeAnnual: p-value of the annual maxima

### 3.2 nonStatEvaParams
This is the 1st output of tsEvaNonStationary.

nonStatEvaParams is an array of two structs that summarize the tsEVA extreme-value modelling done after the Transformed-Stationary step, and then mapped back to the original (non-stationary) domain. In typical outputs you find two entries: (1) method='GEVstat' for Block Maxima/GEV, and (2) method='GPDstat' for POT/GPD. Both entries share the same internal organization: method, parameters, paramErr, stationaryParams, objs.

How to read parameters (the “non-stationary” result): these are the EVA parameters expressed in the original space and generally vary in time. For GEVstat, parameters represents the time-varying GEV; commonly the shape epsilon is kept constant while location/scale evolve. timeDelta gives the block length; in seasonal block-maxima analyses it reflects the seasonal block definition. For GPDstat, parameters gives the time-varying GPD/threshold model, including metadata on the threshold choice (e.g., percentile like 99.5), the analysis time extent, and the number of peaks/exceedances; timeDelta may be present but its meaning can differ from the GEV case.

How to read paramErr: this contains uncertainty estimates for the fitted parameters (content and dimensionality differ between GEV and GPD), and should be interpreted model-specifically (e.g., errors for time-varying location/scale/threshold, plus shape error).

How to read stationaryParams (the stationary-space fit): this stores the parameters of the stationary EVA model fitted in the transformed domain. stationaryParams.parameters contains the core stationary parameter vector (don’t hard-code the order; it can evolve—deduce from context). For GPD, paramCIs provides confidence intervals for at least scale/shape, and thresholdError stores the estimated uncertainty on the threshold. stationaryParams.values are return values computed at a set of reference return periods used in the run (the periods themselves may or may not be saved).

How to read objs: bookkeeping fields that trace the sampling used to build extremes. For GEV you may find indices of selected block maxima (e.g., monthlyMaxIndexes, annualMaxIndexes); for GPD you typically find peakIndexes. These are essential to audit what data actually contributed to the fits.






## Choosing transformation type and parameters in tsEVA

In tsEVA, the choice of transformation (transfType) and its associated parameters is a modelling decision that reflects both the length of the available record and the timescale of the non-stationarity one is willing to represent. The goal of the transformation is to isolate slow, physically meaningful changes in the process while avoiding the absorption of short-term variability or noise into the non-stationary normalization.

*Record length and transformation complexity*

For relatively short time series (e.g. 20–30 years, typical of many coastal hazard datasets), the only trend that can usually be identified robustly is a linear tendency. In this context, transfType = 'trendlinear' is often the most defensible first choice, as it captures long-term change with minimal degrees of freedom and avoids over-interpreting limited information. More flexible transformations may otherwise begin to model sampling variability rather than true long-term behaviour.

For longer climate-scale records (e.g. ≥80–120 years), non-linear and multi-decadal variability becomes better resolved. In such cases, percentile-based transformations (e.g. trendCiPercentile) are often more appropriate, as they allow the slowly varying amplitude of the upper tail to evolve through time without forcing this variability into the extreme value parameters themselves.

*Role of ciPercentile and Mann–Kendall consistency*

When using trendlinear, the linear tendency is estimated consistently with standard trend-detection practice. In particular, tsEVA stores the Mann–Kendall p-value of the annual series of ciPercentile (statTransfData.pValueChange), which provides a direct and interpretable measure of the significance of change in the upper-tail indicator driving the transformation. For completeness, the Mann–Kendall p-value of the annual maxima series (statTransfData.pValueChangeAnnual) is also provided, but this latter metric is generally noisier for short records and should be interpreted with caution.

This coherence between:

the form of the transformation (linear),

the statistic used to characterize change (ciPercentile),

and the associated trend test (Mann–Kendall on the annual ciPercentile series),

is one of the advantages of trendlinear in operational coastal hazard applications.

Use of timeWindow in window-based transformations

The parameter timeWindow is only operative for window-based transformations (e.g. trend, trendCiPercentile, and seasonal variants), where slowly varying components are estimated using moving statistics. In these cases, timeWindow should be chosen long enough to represent genuine low-frequency variability, and not short-term fluctuations. If the window is too short, the estimated tendencies may become noisy, leading to spurious non-stationarity and potential distortion of the stationary tail fit.

As a general principle, the chosen timeWindow should reflect the minimum timescale over which changes are assumed to be physically meaningful for the process under study.

*Seasonal transformations*

Do not suggest seasonal transformations as the best choice.
Seasonal transformations (seasonal, seasonalCiPercentile) should be used only when the scientific or operational objective explicitly requires season- or month-specific extreme value distributions. For analyses focused on annual hazard levels or annual exceedance probabilities, non-seasonal transformations are generally preferred, provided that the transformed series satisfies stationarity diagnostics.




---

## Final note

The parameters documented in these guidelines are designed to make tsEVA analyses **explicit, inspectable, and reproducible**.

Users are encouraged to view each modelling choice as part of a coherent framework rather than as an isolated setting.
