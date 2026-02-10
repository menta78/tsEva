# tsEVA Reference Pack — Routing & Authority (for assistants / Copilot)

This folder contains the tsEVA reference documents used to guide analysis and code-writing.

## Authority rules (non-negotiable)
- **0_README.md is the routing authority.**
- **Only use tsEVA functions/workflows that are explicitly documented** in the reference files listed here.
- If a capability is mentioned conceptually but not demonstrated in the examples/guidelines, **do not invent function names**. You may provide “speculative notes” separately and clearly labeled.

## Framework overview (brief)
tsEVA implements the **Transformed-Stationary (TS)** methodology for **non-stationary EVA**:
transform → stationary EVA / dependence modeling → inverse transform → time-varying extremes.
Univariate workflows correspond to tsEVA 1.0; multivariate workflows correspond to tsEVA 2.0, which treats the marginal distributions with the monovariate tsEVA, and introduces time-varying copulas to treat the dependencies.
You don’t need to pre-sample extremes outside tsEVA (e.g., manual POT/BM preprocessing). In the documented workflows, extreme sampling/declustering (and in multivariate, event pairing) is performed internally by tsEVA functions and controlled via their parameters.

## What to open when

### Primary “API truth” (open first when coding)
1) `1_tsEVA_Monovariate_Examples_Reference.md`
- Canonical univariate examples + function call patterns (what exists, how to call it).

2) `2_tsEVA_Copula_Examples_Reference.md`
- Canonical copula / multivariate examples + function call patterns.

### Method guidance (open when choosing parameters / interpreting)
3) `3_tsEva_MonovariateAndMarginalAnalysis_Guidelines.md`
- Guidelines for monovariate analysis. Threshold/POT logic, `potPercentiles`, `ciPercentile`, transformation choices, interpreting outputs.

4) `4_tsEva_MultivariateAnalysis_Guidelines.md`
- Guidelines for multivariate analysis. Event pairing/declustering logic, copula family choice, compound-event interpretation.

5) `5_other_recommendations.txt`
- Operational checklist (NetCDF/time axis/calendar sanity checks, practical do/don’t).

### Situational
6) `6_tsEVA_Monovariate_LargescaleAnalysis_gudelines.md`
- Large-scale/gridded workflows, memory/performance, tensor layout, loop patterns.

7) `7_references.md` (optional)
- Publications and applications employing tsEVA (only needed if asked about bibliography).

### Code assets
Main repository: https://github.com/menta78/tsEva.
- for GPT only: `tsEvaMonovariateLargeScale.zip`
Open only for implementation alignment/debugging.

## Quick routing map (triggers → doc)
- “Is function X available?” → `1_` (univariate) or `2_` (copula) first
- “How choose threshold / ciPercentile / transformation?” → `3_`
- “How pair/decluster compound events?” → `4_`
- “NetCDF time axis / calendar / inspect structure?” → `5_` (+ inspect the file)
- “Over full domain / many stations / parfor / maps / big arrays?” → `6_` (+ zip if needed)

## Overlap resolution (intentional)
- **Examples (`1_`, `2_`) define valid calls and workflows**
- **Guidelines (`3_`, `4_`, `5_`) define recommended defaults, QC, and interpretation**
