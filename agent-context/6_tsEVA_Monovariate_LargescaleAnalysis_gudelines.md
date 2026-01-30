## 1. Output organization (large-scale monovariate tsEVA)

Large-scale monovariate tsEVA runs (gridded or many stations) should produce outputs that are:
- **map-friendly** (easy to reshape to lon×lat),
- **subset-friendly** (easy to extract a given return period and/or year),
- **self-consistent** (fixed dimensions, missing values handled uniformly),
- **self-describing** (axes and metadata stored with the product).

The guiding principle is that return levels are a function of:
- **space** (grid cell / station / point index),
- **return period** (years),
- **evaluation time** (only for non-stationary runs).

---

### 1.1 Axes and tensor layout

For each point `p`, the non-stationary return level can be written as:

`RL(p, T, t)`

where:
- `T` is the return period in years (`returnPeriodsInYears`)
- `t` is the evaluation time (`retLevTimeStamps`, often ≤ 1 timestamp per year)

Recommended storage layouts:

- **Unstructured point set** (most general):  
  `RL(npt, nRP, nEval)`

- **Gridded domain** (if you store lon/lat explicitly):  
  `RL(nLon, nLat, nRP, nEval)`  
  (or store flattened `npt = nLon*nLat` and reshape later)

with:
- `npt = number of points`
- `nRP = length(returnPeriodsInYears)`
- `nEval = length(retLevTimeStamps)` (or `length(retLevYears)`)

---

### 1.2 Evaluation-time axis: `retLevYears` → `retLevTimeStamps`

In the TS (Transformed-Stationary) framework, the fitted distribution evolves smoothly. For large-scale products, it is usually sufficient to store return levels at **≤ 1 timestamp per year**.

A common convention is to evaluate at **January 1st** of each year:

```matlab
retLevYears = (1970:2100)';  % example
retLevDtvec = [retLevYears, ones(size(retLevYears)), ones(size(retLevYears))];
retLevTimeStamps = datenum(retLevDtvec);  % Jan-01 of each year
```

Notes:
- `retLevYears` defines the **output time grid** (storage + temporal resolution).
- If stationary, you can keep compatibility by setting `nEval = 1` (single timestamp).

---

### 1.3 Return-period axis: `returnPeriodsInYears`

Choose return periods based on the intended application (risk mapping, design levels, adaptation planning). Store them as a 1D axis:

```matlab
returnPeriodsInYears = [1.5 2 3 4 5 7 10 15 20 30 50 70 100 150 250 350 500 700 ...
                        1000 1500 2000 5000 10000];
```

Guidelines:
- Keep them **sorted ascending**.
- Very large return periods (e.g., 5000–10000y) are possible but represent strong extrapolation; they should always be accompanied by uncertainty layers.

---

### 1.4 Stationary vs non-stationary products

Typically, large-scale applications target the **non-stationary** product. If support for the **stationary** analysis is needed, keep the output structure consistent.

**Non-stationary run**
- Return levels depend on evaluation time:
  - `RL(…, nRP, nEval)`
  - `RLerr(…, nRP, nEval)` (recommended)

**Stationary run**
- Return levels do not depend on time.
- For interoperability with non-stationary outputs, prefer a unified structure:
  - use `nEval = 1` (single evaluation timestamp) so arrays remain `(..., nRP, nEval)`.

---

### 1.5 Preallocation of output arrays (mandatory for large scale / `parfor`)

Before entering the main loop (and especially before any `parfor`), preallocate all output arrays at their final size. This:
- avoids dynamic array growth,
- enables `parfor` sliced assignments,
- ensures failed points remain `NaN` without special handling,
- guarantees a consistent product shape across the domain.

Use `NaN` initialization with the exact output axes:

```matlab
nRP   = length(returnPeriodsInYears);
nEval = length(retLevTimeStamps);

% Return levels and uncertainties: (point × returnPeriod × evalTime)
retLevGEV         = ones([npt, nRP, nEval]) * nan;
retLevErrGEV      = ones([npt, nRP, nEval]) * nan;
retLevErrFitGEV   = ones([npt, nRP, nEval]) * nan;
retLevErrTransGEV = ones([npt, nRP, nEval]) * nan;

retLevGPD         = ones([npt, nRP, nEval]) * nan;
retLevErrGPD      = ones([npt, nRP, nEval]) * nan;
retLevErrFitGPD   = ones([npt, nRP, nEval]) * nan;
retLevErrTransGPD = ones([npt, nRP, nEval]) * nan;

% Parameters:
% - time-invariant (in this product convention): (point × 1)
shapeGEV    = ones([npt, 1]) * nan;
shapeGEVErr = ones([npt, 1]) * nan;

shapeGPD    = ones([npt, 1]) * nan;
shapeGPDErr = ones([npt, 1]) * nan;

% - time-varying (evaluated at each timestamp): (point × evalTime)
scaleGEV    = ones([npt, nEval]) * nan;
scaleGEVErr = ones([npt, nEval]) * nan;
locGEV      = ones([npt, nEval]) * nan;
locGEVErr   = ones([npt, nEval]) * nan;

scaleGPD        = ones([npt, nEval]) * nan;
scaleGPDErr     = ones([npt, nEval]) * nan;
thresholdGPD    = ones([npt, nEval]) * nan;
thresholdGPDErr = ones([npt, nEval]) * nan;
```

Recommendations:
- Keep naming explicit: `Err` (total), `ErrFit`, `ErrTrans` when you separate uncertainty contributions.
- If a point is invalid (too many missing values, too few extremes, failed fit), leave outputs as `NaN` and (optionally) store a `fit_ok(npt,1)` flag.

---

### 1.6 Missing values and masks

- Use `NaN` for missing/invalid values in MATLAB arrays.
- In NetCDF, use a consistent `_FillValue` and propagate the mask.
- A point that fails should be represented by all-`NaN` outputs (plus `fit_ok = 0` if included).

---

## 2. Large-scale execution loop (parallel over points) and per-point products

This section describes the recommended *execution pattern* for large-scale monovariate tsEVA: a `parfor` loop over points, where each worker reads one time series, runs tsEVA, optionally reduces the analysis-object size for storage, computes return levels on the requested axes, and writes results into preallocated arrays.

> **Separation of roles**
> - **tsEVA core steps (documented workflow):** build `timeAndSeries`, run `tsEvaNonStationary` (or `tsEvaStationary`), compute return levels from the analysis object, extract parameters.
> - **Project utilities (I/O helpers):** functions like `getJoinedSeries`, `getSeriesFileName`, `parsave`, and the output-reduction helper `tsEvaReduceOutputObjSize` are part of the *large-scale pipeline utilities*, not the scientific core. They are used to make runs feasible at scale (I/O, storage, diagnostics).

---

### 2.1 Why `parfor` over points

Large-scale EVA is “embarrassingly parallel”: each grid cell/station can be processed independently. The recommended pattern is:

- **parallelize over points** (`parfor ipt = 1:npt`)
- keep the **time dimension internal** to the point’s tsEVA call
- write into **preallocated** arrays using sliced indexing (`retLevGEV(ipt,:,:) = ...`)

This is the most robust approach for performance and memory control.

**Debug vs production.** It is recommended to debug with a standard `for` loop (deterministic behavior, easier breakpoints, simpler logs). For production runs, switching to `parfor` is typically much faster. In this workflow the two are intentionally interchangeable: you can switch between `for` and `parfor` by commenting/uncommenting a single line, provided that outputs are **preallocated** and assignments use **sliced indexing** (e.g., `out(ipt,:,:) = ...`), which is required for `parfor`.

---

### 2.2 Per-point workflow inside the loop

Each iteration follows the same steps:

1) **Identify point location** (lon/lat or station id)  
2) **Read the time series** for that point (project-specific I/O helper)  
3) **Apply optional time horizon subsetting**  
4) **Apply quality filters** (remove invalid values; skip empty/constant series)  
5) **Run EVA**
   - `tsEvaNonStationary(timeAndSeries, timeWindow, ...)` or  
   - `tsEvaStationary(timeAndSeries, ...)`
6) **Skip invalid fits** (`if ~isValid, continue; end`)
7) **Optionally reduce analysis-object size** (for saving diagnostics at scale)
8) **Compute return levels** for both GEV and GPD from the analysis object
9) **Store results** into preallocated output arrays
10) **Optionally save per-point reduced analysis details** for later inspection/debug

---

### 2.3 Output reduction for scalable diagnostics (`tsEvaReduceOutputObjSize`)

Large-scale runs often need the option to **save per-point analysis objects** (`nonStatEvaParams`, `statTransfData`) for later diagnosis, without exploding disk usage.

A practical strategy is to call the helper:

```matlab
[nonStatEvaParams, statTransfData] = tsEvaReduceOutputObjSize( ...
    nonStatEvaParams, statTransfData, retLevTimeStamps);
```

**Intent of this step**
- retain only what is needed for:
  - recomputing return levels at `retLevTimeStamps`,
  - plotting/debugging at selected points,
  - lightweight checking of parameters and diagnostics
- drop large intermediate arrays that are not needed for later QA/QC

This reduction is typically performed **only for non-stationary runs** and **after** checking `isValid`.

---

### 2.4 Computing return levels (GEV and GPD) on the requested axes

After (optional) reduction, return levels are computed for each return period in `returnPeriodsInYears`, and for each evaluation timestamp in `retLevTimeStamps` (implicit in the analysis object).

Typical pattern:

```matlab
% GEV return levels
[RLgev, RLgevErr, RLgevErrFit, RLgevErrTrans] = ...
    tsEvaComputeReturnLevelsGEVFromAnalysisObj(nonStatEvaParams, returnPeriodsInYears);

retLevGEV(ipt,:,:)         = RLgev';
retLevErrGEV(ipt,:,:)      = RLgevErr';
retLevErrFitGEV(ipt,:,:)   = RLgevErrFit';
retLevErrTransGEV(ipt,:,:) = RLgevErrTrans';

% GPD return levels
[RLgpd, RLgpdErr, RLgpdErrFit, RLgpdErrTrans] = ...
    tsEvaComputeReturnLevelsGPDFromAnalysisObj(nonStatEvaParams, returnPeriodsInYears);

retLevGPD(ipt,:,:)         = RLgpd';
retLevErrGPD(ipt,:,:)      = RLgpdErr';
retLevErrFitGPD(ipt,:,:)   = RLgpdErrFit';
retLevErrTransGPD(ipt,:,:) = RLgpdErrTrans';
```

**Shape convention**
- `returnLevels` are typically returned as `(nEval × nRP)` and are transposed into `(nRP × nEval)` to match the preallocated `retLev*(ipt, nRP, nEval)` convention.

---

### 2.5 Storing parameters (shape vs time-varying parameters)

The large-scale product often stores:
- a **single** shape parameter per point (e.g., `shapeGEV(ipt)`), and
- time-varying parameters per point and evaluation time (e.g., `scaleGEV(ipt,:)`).

Example extraction (as used in the pipeline):

```matlab
% GEV parameters
shapeGEV(ipt)    = nonStatEvaParams(1).parameters.epsilon;
shapeGEVErr(ipt) = nonStatEvaParams(1).paramErr.epsilonErr;

scaleGEV(ipt,:)    = nonStatEvaParams(1).parameters.sigma;
scaleGEVErr(ipt,:) = nonStatEvaParams(1).paramErr.sigmaErr;

locGEV(ipt,:)    = nonStatEvaParams(1).parameters.mu;
locGEVErr(ipt,:) = nonStatEvaParams(1).paramErr.muErr;

% GPD parameters
shapeGPD(ipt)    = nonStatEvaParams(2).parameters.epsilon;
shapeGPDErr(ipt) = nonStatEvaParams(2).paramErr.epsilonErr;

scaleGPD(ipt,:)    = nonStatEvaParams(2).parameters.sigma;
scaleGPDErr(ipt,:) = nonStatEvaParams(2).paramErr.sigmaErr;

thresholdGPD(ipt,:)    = nonStatEvaParams(2).parameters.threshold;
thresholdGPDErr(ipt,:) = nonStatEvaParams(2).paramErr.thresholdErr;
```

---

### 2.6 Optional per-point saving of reduced analysis details (debug/QC)

To support future checking without rerunning the full analysis, store the (reduced) analysis objects per point:

```matlab
if saveAnalysisDetails
  fprintf(1, 'saving ts eva data to ...\n');
  fprintf(1, ['  ' evaFlPath '\n']);
  parsave(evaFlPath, nonStatEvaParams, statTransfData);
end
```

Notes:
- In `parfor`, direct `save(...)` patterns are often replaced by a helper like `parsave(...)` to avoid issues with workspace scoping and to standardize file naming.
- Prefer **one file per point** (or per chunk of points) to prevent write contention and to allow selective post-mortem inspection.

---

### 2.X Parallel pool lifecycle and safe cleanup (`parpool`, `try/catch`, `delete`)

Large-scale runs often create a parallel pool explicitly to control how many workers are used and to ensure the pool is properly shut down even if an error occurs. A robust pattern is:

```matlab
if createParObj && (nParWorker > 1)
  parObj = parpool(nParWorker);
else
  parObj = [];
end

try
  % ... large-scale processing (possibly using parfor) ...

catch exc
  % If the run fails, close the pool we created (do not leave resources dangling)
  if createParObj
    delete(parObj);
  end
  rethrow(exc);
end

% Normal completion: close the pool we created
if createParObj
  delete(parObj);
end
```

**Meaning of the flags and variables**
- `createParObj`: boolean that indicates whether this script should **create and own** the parallel pool.
  - `true`: the script creates the pool and is responsible for deleting it.
  - `false`: the script assumes the user/batch environment manages the pool (or no pool is desired).
- `nParWorker`: requested number of workers. The pool is created only if `nParWorker > 1`.
- `parObj`: handle to the pool object (empty if no pool was created).

**Why this pattern is used**
- **Control:** `parpool(nParWorker)` enforces the requested worker count instead of relying on whatever pool state happens to exist.
- **Resource hygiene:** on HPC/batch runs, leaving a pool open after an exception can waste allocated resources; the `catch` block prevents that.
- **Non-interference:** conditioning cleanup on `createParObj` avoids deleting a pool that was started outside this script.

**Important note**
- `rethrow(exc)` preserves the original error and stack trace, which is essential for debugging large-scale failures.

---

### 2.7 Practical runtime hygiene inside `parfor`

- Minimize heavy console I/O (`fprintf`) inside `parfor` (it can become a bottleneck).
- Keep point-level failures non-fatal: if a fit fails, `continue` and leave outputs as `NaN`.
- Periodically call `drawnow('update')` (or similar) only if you observed stdout buffering issues in your batch environment.

---

## 3. Saving the large-scale results (NetCDF product)

After the parallel loop completes, the analysis results exist as **preallocated MATLAB arrays** (mostly `NaN`-filled, with valid points populated). Section 3 defines a **portable, self-describing output format** and writes the product to disk. In this workflow, the chosen format is **NetCDF** because it is:
- language-agnostic (MATLAB/Python/R),
- efficient for large multidimensional arrays,
- standard in geosciences and easy to post-process and map.

---

### 3.1 Output file handling (overwrite safely)

The script announces the output and ensures a clean write by deleting any existing file:

```matlab
disp('parallel loop complete. Saving the return levels to the output netcdf');
disp(retLevNcOutFilePath);

if exist(retLevNcOutFilePath, 'file')
  delete(retLevNcOutFilePath);
end
```

Guideline:
- Always **overwrite** the target file explicitly to avoid mixing old and new dimensions/variables.

---

### 3.2 Define axis sizes once

The product axes are:
- `pt` (point index): `npt`
- `return_period`: `nretper`
- `rlyear` (evaluation year): `nRlTime`

```matlab
nRlTime = length(retLevYears);
nretper = length(returnPeriodsInYears);
```

These sizes must match the preallocated arrays from Section 1 (e.g., `retLevGEV(npt,nretper,nRlTime)`).

---

### 3.3 Write model outputs as NetCDF variables

The script creates NetCDF variables and writes the corresponding MATLAB arrays. The dimensions are explicitly encoded in the `nccreate(... 'dimensions', {...})` call so the file is self-describing.

#### 3.3.1 Return levels and uncertainty (GEV)

Variables (all dimensioned as `('pt','return_period','rlyear')`):
- `returnlevelGEV`
- `returnlevelErrorGEV` (total)
- `returnlevelFitErrorGEV` (fit component)
- `returnlevelTransErrorGEV` (transformation component)

```matlab
nccreate(retLevNcOutFilePath, 'returnlevelGEV', ...
  'dimensions', {'pt', npt, 'return_period', nretper, 'rlyear', nRlTime});
ncwrite(retLevNcOutFilePath, 'returnlevelGEV', retLevGEV);

nccreate(retLevNcOutFilePath, 'returnlevelErrorGEV', ...
  'dimensions', {'pt', npt, 'return_period', nretper, 'rlyear', nRlTime});
ncwrite(retLevNcOutFilePath, 'returnlevelErrorGEV', retLevErrGEV);

nccreate(retLevNcOutFilePath, 'returnlevelFitErrorGEV', ...
  'dimensions', {'pt', npt, 'return_period', nretper, 'rlyear', nRlTime});
ncwrite(retLevNcOutFilePath, 'returnlevelFitErrorGEV', retLevErrFitGEV);

nccreate(retLevNcOutFilePath, 'returnlevelTransErrorGEV', ...
  'dimensions', {'pt', npt, 'return_period', nretper, 'rlyear', nRlTime});
ncwrite(retLevNcOutFilePath, 'returnlevelTransErrorGEV', retLevErrTransGEV);
```

#### 3.3.2 Parameters (GEV)

Parameter conventions:
- `shapeParamGEV(pt)` is stored as time-invariant in this product.
- `scaleParamGEV(pt,rlyear)` and `locParamGEV(pt,rlyear)` vary with evaluation time.

```matlab
nccreate(retLevNcOutFilePath, 'shapeParamGEV', 'dimensions', {'pt', npt});
ncwrite(retLevNcOutFilePath, 'shapeParamGEV', shapeGEV);

nccreate(retLevNcOutFilePath, 'shapeParamGEVErr', 'dimensions', {'pt', npt});
ncwrite(retLevNcOutFilePath, 'shapeParamGEVErr', shapeGEVErr);

nccreate(retLevNcOutFilePath, 'scaleParamGEV', 'dimensions', {'pt', npt, 'rlyear', nRlTime});
ncwrite(retLevNcOutFilePath, 'scaleParamGEV', scaleGEV);

nccreate(retLevNcOutFilePath, 'scaleParamGEVErr', 'dimensions', {'pt', npt, 'rlyear', nRlTime});
ncwrite(retLevNcOutFilePath, 'scaleParamGEVErr', scaleGEVErr);

nccreate(retLevNcOutFilePath, 'locParamGEV', 'dimensions', {'pt', npt, 'rlyear', nRlTime});
ncwrite(retLevNcOutFilePath, 'locParamGEV', locGEV);

nccreate(retLevNcOutFilePath, 'locParamGEVErr', 'dimensions', {'pt', npt, 'rlyear', nRlTime});
ncwrite(retLevNcOutFilePath, 'locParamGEVErr', locGEVErr);
```

#### 3.3.3 Return levels and uncertainty (GPD)

Same structure as GEV, but stored in the `GPD` variables:

```matlab
nccreate(retLevNcOutFilePath, 'returnlevelGPD', ...
  'dimensions', {'pt', npt, 'return_period', nretper, 'rlyear', nRlTime});
ncwrite(retLevNcOutFilePath, 'returnlevelGPD', retLevGPD);

nccreate(retLevNcOutFilePath, 'returnlevelErrorGPD', ...
  'dimensions', {'pt', npt, 'return_period', nretper, 'rlyear', nRlTime});
ncwrite(retLevNcOutFilePath, 'returnlevelErrorGPD', retLevErrGPD);

nccreate(retLevNcOutFilePath, 'returnlevelFitErrorGPD', ...
  'dimensions', {'pt', npt, 'return_period', nretper, 'rlyear', nRlTime});
ncwrite(retLevNcOutFilePath, 'returnlevelFitErrorGPD', retLevErrFitGPD);

nccreate(retLevNcOutFilePath, 'returnlevelTransErrorGPD', ...
  'dimensions', {'pt', npt, 'return_period', nretper, 'rlyear', nRlTime});
ncwrite(retLevNcOutFilePath, 'returnlevelTransErrorGPD', retLevErrTransGPD);
```

#### 3.3.4 Parameters (GPD)

```matlab
nccreate(retLevNcOutFilePath, 'shapeParamGPD', 'dimensions', {'pt', npt});
ncwrite(retLevNcOutFilePath, 'shapeParamGPD', shapeGPD);

% NOTE (legacy scripts): shapeParamGPDErr is sometimes defined as (pt,rlyear),
% while shapeGPDErr is stored as (pt,1).
% In these guidelines we treat the shape-parameter uncertainty as time-invariant
% and store it as shapeParamGPDErr(pt). If you adopt a time-varying convention,
% allocate and write shapeGPDErr(pt,rlyear) consistently.
nccreate(retLevNcOutFilePath, 'shapeParamGPDErr', 'dimensions', {'pt', npt});
ncwrite(retLevNcOutFilePath, 'shapeParamGPDErr', shapeGPDErr);

nccreate(retLevNcOutFilePath, 'scaleParamGPD', 'dimensions', {'pt', npt, 'rlyear', nRlTime});
ncwrite(retLevNcOutFilePath, 'scaleParamGPD', scaleGPD);

nccreate(retLevNcOutFilePath, 'scaleParamGPDErr', 'dimensions', {'pt', npt, 'rlyear', nRlTime});
ncwrite(retLevNcOutFilePath, 'scaleParamGPDErr', scaleGPDErr);

nccreate(retLevNcOutFilePath, 'thresholdParamGPD', 'dimensions', {'pt', npt, 'rlyear', nRlTime});
ncwrite(retLevNcOutFilePath, 'thresholdParamGPD', thresholdGPD);

nccreate(retLevNcOutFilePath, 'thresholdParamGPDErr', 'dimensions', {'pt', npt, 'rlyear', nRlTime});
ncwrite(retLevNcOutFilePath, 'thresholdParamGPDErr', thresholdGPDErr);
```

**Important consistency check:** ensure that each `nccreate(...dimensions...)` matches the actual shape of the MATLAB array being written (for every variable), otherwise `ncwrite` will error or silently permute/truncate in unexpected ways.

---

### 3.4 Save coordinates and axes (required for self-description)

Finally, the file stores:
- point coordinates `lon(pt)`, `lat(pt)`
- point index `pt`
- return period axis `return_period(nretper)`
- evaluation-year axis `rlyear(nRlTime)`

```matlab
nccreate(retLevNcOutFilePath, 'lon', 'dimensions', {'pt', npt});
ncwrite(retLevNcOutFilePath, 'lon', lon);

nccreate(retLevNcOutFilePath, 'lat', 'dimensions', {'pt', npt});
ncwrite(retLevNcOutFilePath, 'lat', lat);

nccreate(retLevNcOutFilePath, 'pt', 'dimensions', {'pt', npt});
ncwrite(retLevNcOutFilePath, 'pt', 1:npt);

nccreate(retLevNcOutFilePath, 'return_period', 'dimensions', {'return_period', nretper});
ncwrite(retLevNcOutFilePath, 'return_period', returnPeriodsInYears);

nccreate(retLevNcOutFilePath, 'rlyear', 'dimensions', {'rlyear', nRlTime});
ncwrite(retLevNcOutFilePath, 'rlyear', retLevYears);
```

Guideline:
- Always include these axes so downstream users can map and subset results without external context.
- Prefer using the same names for coordinate variables and their dimensions (e.g., `pt`, `return_period`, `rlyear`) so the NetCDF structure is unambiguous.

---

### 3.5 Recommended additions (small but valuable)

For robustness and interpretability, consider also storing:
- `fit_ok(pt)` (boolean) and/or counts of peaks/block maxima used,
- key configuration metadata (e.g., `transfType`, `timeWindow`, `minPeakDistanceInDays`, POT percentile settings),
- a consistent missing value convention (e.g., `_FillValue`) if needed by downstream tools.

A minimal, practical approach is to write key run settings as global attributes:

```matlab
ncwriteatt(retLevNcOutFilePath,'/','transfType', transfType);
ncwriteatt(retLevNcOutFilePath,'/','timeWindow_days', timeWindow);
ncwriteatt(retLevNcOutFilePath,'/','minPeakDistanceInDays', minPeakDistanceInDays);
ncwriteatt(retLevNcOutFilePath,'/','returnPeriodsInYears', mat2str(returnPeriodsInYears));
ncwriteatt(retLevNcOutFilePath,'/','retLevYears', mat2str(retLevYears(:)'));
```

---

## 4. Running large-scale tsEVA on memory-limited machines (streaming and chunking)

tsEVA itself is typically light enough for large-scale runs, but memory can become a bottleneck when the input dataset is too large to load at once, or when parallel workers replicate data structures. This workflow supports memory-limited environments as long as the implementation ensures that **only the time series needed for the current location (or current chunk of locations) is resident in memory**.

---

### 4.1 Symptom-driven troubleshooting checklist

Common failure modes:
- **Out-of-memory** errors while loading the domain (input arrays too large).
- Memory blow-up when starting `parpool` / entering `parfor` (worker replication of variables).
- Gradual growth over time (accidental accumulation of large temporary variables, or saving too much per point).

Quick checks:
- Reduce `nParWorker`: each worker needs its own working memory.
- Ensure **no large arrays** (e.g., full `tas(time,lat,lon)`) are created outside the loop.
- Keep point-level variables local and small; do not store raw series for all points.

---

### 4.2 Recommended baseline: on-demand per-point loading inside the loop

The safest approach (lowest memory) is:
- in each iteration, **load only one point time series**,
- run `tsEvaNonStationary` (or `tsEvaStationary`),
- write results into preallocated outputs,
- discard temporary variables.

Key principle:
- Input time series is loaded in memory **only while processing that location**.

Implementation notes:
- Prefer I/O patterns that support reading a single point efficiently (e.g., NetCDF hyperslabs via `start/count`, per-point MAT files, or a lightweight accessor function).
- Avoid “broadcasting” large read-only arrays into the `parfor` workspace.

---

### 4.3 Chunked workflow: load a block of series, then `parfor` over the block

If per-point I/O is slow (many small reads), an alternative is to load a *chunk* of points (or a tile of the grid) into memory, then process that chunk with `parfor`.

High-level pattern:
1) Choose a chunk size that fits in RAM (e.g., `nChunkPts × nTime`).
2) Load the chunk series into a local array `Xchunk`.
3) Run `parfor` over points **within the chunk**, reading series from `Xchunk`.
4) Write results into the preallocated global outputs using the global indices.
5) Clear the chunk and move to the next one.

Benefits:
- Fewer I/O calls (faster on network filesystems).
- Memory remains bounded by the chunk size.

Rule of thumb:
- Chunk so that `(nChunkPts × nTime × 8 bytes)` plus overhead stays comfortably below available RAM,
  then reduce further if using many workers (each worker adds overhead).

---

### 4.4 Minimize per-worker memory in `parfor`

`parfor` can replicate variables across workers. To keep memory stable:
- Keep large metadata out of the `parfor` body unless needed.
- Preallocate outputs once (Section 1) and only write sliced pieces: `out(ipt,:,:) = ...`.
- Avoid collecting per-point diagnostics in growing arrays inside `parfor`.

If `saveAnalysisDetails` is enabled:
- Use output reduction (e.g., `tsEvaReduceOutputObjSize(...)`) before saving analysis objects.
- Save **one point per file** (or one chunk per file) to avoid holding many objects in memory.

---

### 4.5 Practical levers to solve memory issues quickly

If you hit memory limits, apply these in order (usually the first two are enough):

1) **Reduce parallelism**  
   Decrease `nParWorker`. Memory usage often scales roughly with the number of workers.

2) **Ensure streaming input**  
   Do not load the full domain array. Load only per-point (or per-chunk) series.

3) **Store fewer/heavier things**  
   - Disable `saveAnalysisDetails` unless needed.
   - If enabled, reduce objects before saving and avoid saving raw series.

4) **Use lighter numeric types for large outputs**  
   Consider `single` for large `retLev*` arrays if precision requirements allow (this is a product decision).

5) **Write incrementally for very large domains**  
   If outputs themselves become too large, write results in blocks (e.g., chunk-by-chunk NetCDF writes or per-chunk MAT files) and merge later.

---

### 4.6 How to diagnose what is causing the bottleneck

A simple diagnostic strategy is to run:
- a small subset of points with `for` (debug mode) and inspect memory usage,
- then enable `parfor` and vary `nParWorker`,
- then switch between per-point loading vs chunk loading.

Interpretation:
- If memory spikes at `parpool` start or at `parfor` entry → replication/worker overhead is dominant → reduce workers and avoid broadcasting large variables.
- If memory spikes during data reads → you are loading too much data at once → enforce per-point or per-chunk loading.
- If memory grows gradually over the loop → you are accumulating data unintentionally (e.g., saving full objects/series, or appending inside the loop).

---

## 5. Recommended practice for climate projection datasets (CMIP/CORDEX and derived products)

A common tsEVA application is the analysis of climate-model projections (e.g., CMIP, CORDEX, or downscaled/derived datasets). In these contexts, the goal is usually to quantify how extremes evolve under forcing scenarios, while properly representing both **internal variability** and **model uncertainty**.

---

### 5.1 Join historical + scenario into a single continuous time series

**Recommendation:** run tsEVA on a time series that concatenates the **baseline historical** period with the **future scenario** period (RCP/SSP), forming one continuous record.

Rationale:
- tsEVA’s non-stationary inference benefits from a longer record to constrain the transformation and the tail model.
- It avoids artificial discontinuities that may arise if historical and scenario periods are treated separately.
- It supports a consistent mapping of return levels across the full evaluation timeline (e.g., 1970–2100).

Practical note:
- Ensure time stamps are continuous and correctly ordered.
- Apply consistent QC/masking rules across the joined series.

---

### 5.2 Use multi-model ensembles (do not rely on a single model)

**Recommendation:** in climate-change applications, analyze **multiple models**, not just one.

Rationale:
- A single model rarely captures the full range of plausible future evolution of extremes.
- Structural differences among models (dynamics, parameterizations, resolution, downscaling) often dominate uncertainty in projected extremes.
- Multi-model analysis is typically needed to make defensible statements about changes in return levels or return periods.

Operational implication:
- The large-scale workflow should support looping over models (and scenarios) as outer loops, with the spatial `parfor` as the inner scalable computation.

---

### 5.3 Run tsEVA separately for each model (never on the ensemble mean)

**Recommendation:** if you have an ensemble, run tsEVA **independently for each model member**; do **not** run tsEVA on the ensemble mean time series.

Rationale:
- The ensemble mean is not a physically realizable trajectory and typically has **reduced variance**, especially in high-frequency variability.
- Extremes are **nonlinear** statistics; taking the mean first generally biases extremes low and distorts tail behavior.
- The correct workflow is:
  1) fit tsEVA per model → obtain per-model return levels through time,
  2) aggregate the resulting return levels across models (e.g., median, quantiles) to summarize ensemble uncertainty.

Practical note:
- Store per-model outputs in separate files (or separate dimensions) so that ensemble statistics can be computed transparently and reproducibly later.

---

## 6. Companion example code (optional)

A companion example implementation of the large-scale monovariate workflow is provided in:

- `tsevaMonovariateLargeScale.zip`  
  (scripts: `tsevaPerformModel.m`, `tsevaLoopModels.m`)

These scripts implement an end-to-end pipeline consistent with the guidelines in this document:

- output preallocation (fixed tensor sizes),
- `parfor` loop over points (and the option to switch to `for` for debugging),
- optional per-point saving of reduced analysis details (for diagnosis and reproducibility),
- NetCDF export of return levels, uncertainties, parameters, and axes.

**Important scope note.** The scripts also contain project-specific utilities (e.g., I/O helpers, file naming conventions, saving helpers, filtering rules) that may require adaptation for other datasets and computing environments. Treat them as a practical reference implementation and a starting point for customization, not as a mandatory API contract.
