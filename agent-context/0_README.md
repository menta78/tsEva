# tsEVA 2.0 Reference Pack — What to Read When

This folder contains the tsEVA 2.0 reference documents used to guide analysis and code-writing.
The goal is to keep workflows **scientifically rigorous** and **strictly aligned with documented tsEVA functions/examples**.


---

## Must-read (every chat)
These documents define the “allowed” workflows and the methodological guardrails.

### 1) `1_tsEVA_Monovariate_Examples_Reference.md`  ✅ MUST READ
**Role:** Canonical monovariate templates + function usage patterns.  
**Contains:** 8 monovariate examples (stationary/non-stationary; GEV/GPD; trend/percentile/linear; temperature example), plus “Common Workflow Patterns”.

Use this as:
- the function-availability registry for monovariate EVA
- the source of the canonical call patterns

---

### 2) `2_tsEVA_Copula_Examples_Reference.md` ✅ MUST READ
**Role:** Canonical multivariate/copula templates + function usage patterns.  
**Contains:** Common copula workflow and 3 case studies:
- CaseStudy01: bivariate compound flooding (GPD)
- CaseStudy02: trivariate spatial wave extremes (GPD)
- CaseStudy03: bivariate temperature & drought (GEV)

Use this as:
- the function-availability registry for multivariate EVA
- the source of canonical copula workflows (Gaussian/Gumbel multivariate; Frank bivariate only)

---

### 3) `3_tsEva_MonovariateAndMarginalAnalysis_Guidelines.md` ✅ MUST READ
**Role:** Methodological rules-of-thumb and parameter logic for monovariate/marginal EVA.  
**Core topics:**
- Sampling of extremes
- POT/GPD guidance and `potPercentiles` logic (including winning percentile, uncertainty)
- Non-stationary normalization and `ciPercentile` (why percentile-based normalization)
- Output inspection: `statTransfData`, `nonStatEvaParams`
- Choosing transformation type and parameters

Use this when:
- users ask “how do I choose threshold / ciPercentile / transform type?”
- interpreting outputs and uncertainties

---

### 4) `4_tsEva_MultivariateAnalysis_Guidelines.md` ✅ MUST READ
**Role:** Practical guidance for multivariate event sampling and compound-event interpretation.  
**Core topics:**
- Limits/choice of POT–GPD vs block-maxima–GEV in multivariate settings
- Copula family choice guidance
- Declustering & event pairing parameters:
  - `minPeakDistanceInDaysMonovarSampling`
  - `maxPeakDistanceInDaysMultivarSampling`
- Hazard-type guidance (storm-episode vs impact-based compound events)

Use this when:
- users ask about joint extremes, compound hazards, event pairing windows, declustering logic

---

### 5) `5_other_recommendations.txt` ✅ MUST READ
**Role:** Operational “do/don’t” checklist.  
**Contains:** Critical practical reminders (e.g., inspect NetCDF; ensure time units/calendar awareness).

Use this as:
- the default operational safety checklist at the start of any applied workflow

---

## Read only if needed (situational)

### 6) `6_tsEVA_Monovariate_LargescaleAnalysis_gudelines.md` 🔁 READ ONLY FOR LARGE-SCALE
**Trigger:** “over the whole domain”, gridded products, many stations, `parfor`, map outputs, memory/performance constraints.  
**Role:** Output tensor layout, evaluation-time axis strategy, preallocation rules, missing masks, large-scale loop patterns.

Use this when:
- the user requests domain-wide processing or any large-scale product design

---

## Code assets (open when you need implementation alignment)
### `tsEvaMonovariateLargeScale.zip`
Contains MATLAB scripts used for large-scale monovariate execution (e.g., loop/run helpers).
Only open when:
- aligning an implementation to the provided large-scale scripts, or debugging those scripts.

---

## Quick routing map (which doc drives what)
- **“Is function X available?”** → Examples first (`1_`, `2_`)
- **“How choose threshold / ciPercentile / transform?”** → `3_`
- **“How pair/decluster compound events?”** → `4_`
- **“NetCDF time axis / calendar / inspect structure?”** → `5_` (+ inspect file)
- **“Run over full domain / map product / parfor / big output arrays?”** → `6_`

---

## Notes on overlap (intentional)
Some fundamentals (sampling extremes, POT vs block maxima) appear in both `3_` and `4_`.
Resolution rule:
- **Examples (`1_`/`2_`) define valid calls and workflows**
- **Guidelines (`3_`/`4_`/`5_`) define recommended defaults, QC, and interpretation**
