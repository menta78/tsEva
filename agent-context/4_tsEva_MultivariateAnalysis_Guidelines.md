# Soft guidelines for multivariate sampling parameters in tsEVA

This section provides **soft (non-prescriptive) guidelines** for choosing the
parameters controlling univariate declustering and multivariate event pairing
in the tsEVA 2.0 framework.

These parameters define how *independent events* are identified in each margin
and how they are combined into *compound (multivariate) events*. Their choice
should reflect both **physical process timescales** and the **impact-based
definition of compounding** adopted in the analysis.

---


## Choice and limitations of POT–GPD and block-maxima–GEV in multivariate tsEVA

In tsEVA, the choice between a Peak-Over-Threshold (POT) approach with GPD margins and a block-maxima approach with GEV margins is not a purely statistical decision, but reflects how **extreme events are defined and interpreted** in the context of the scientific or impact-related question.

The POT–GPD approach is generally preferable when extremes are naturally interpreted as **events**, when multiple independent extremes may occur within a block (e.g. within a year), and when the **timing and compounding of events** is relevant. This is typically the case for storm-driven hazards such as storm surge, waves, runup, total water level, and river discharge. POT sampling allows a much larger number of extreme events to be retained, which is advantageous for dependence estimation and copula-based joint analysis. However, in multivariate settings—especially when the number of variables is large or when dependence is heterogeneous—it becomes increasingly unlikely that all variables exceed their thresholds simultaneously. As a result, joint POT samples often combine truly extreme values in some dimensions with non-extreme values in others. While the resulting joint distribution remains statistically coherent and asymptotically valid, its interpretation **close to the threshold** becomes less straightforward, as the implicit assumption that all components are simultaneously extreme is violated. In practice, POT–GPD joint models should therefore be interpreted as most reliable **away from the threshold**, in the extrapolated tail.

The block-maxima–GEV approach avoids this issue, as each joint observation corresponds by construction to a block maximum in all variables. This makes GEV-based joint distributions intrinsically multivariate-extreme and conceptually robust even in higher dimensions or under weak dependence. GEV is therefore well suited for phenomena characterized by **long intrinsic timescales** (e.g. heatwaves, droughts) or for impact frameworks in which extremes are assumed to **accumulate at the block scale** (e.g. annual or seasonal impacts). The main drawback of GEV is that it disregards the physical linkage between individual events within the block and implicitly assumes that extremes occurring in the same block contribute jointly to impact, even if they are not causally or temporally connected. In addition, for hazards that cluster seasonally—particularly winter storm hazards—annual block definitions may artificially split physically connected events across block boundaries (e.g. events occurring on December 31 and January 1). In such cases, redefining the block calendar (e.g. hydrological year or seasonally shifted year) as a preprocessing step is recommended to ensure that block maxima represent physically meaningful compound events.

In summary, POT–GPD and block-maxima–GEV represent complementary approaches with different assumptions and limitations. The choice should be guided by the **event definition**, the **timescale of the processes**, the **dimensionality of the joint analysis**, and the **intended interpretation of compound impacts**, rather than by goodness-of-fit considerations alone.





## Choice of copula families in tsEVA

In practice, tsEVA primarily supports **Gaussian** and **Gumbel** copulas for multivariate extreme value analysis. This choice is intentional and reflects both the objectives of joint hazard analysis and the properties of the sampling strategies adopted in tsEVA.

Gaussian and Gumbel copulas represent the two dominant dependence regimes relevant for environmental extremes. The Gaussian copula provides a parsimonious and interpretable representation of symmetric dependence, while the Gumbel copula is the simplest and most robust model capturing **upper-tail dependence**, which is the key feature of interest for compound extreme events and risk assessment.

Copulas such as Clayton or t, although theoretically admissible, tend to perform well in multivariate POT–GPD settings for reasons that are often unrelated to the physical dependence of extremes. In joint POT sampling, exceedances are defined marginally, and joint samples are typically dominated by events in which only a subset of variables is truly extreme. Copulas with greater flexibility in the bulk or symmetric tail behavior may therefore achieve better likelihood-based scores by fitting the **sampling structure induced by POT**, rather than the dependence of joint extremes themselves. Their apparent superiority is thus frequently an artifact of the sampling mechanism, and does not necessarily translate into more reliable or interpretable joint tail inference.

By restricting the copula choice to Gaussian and Gumbel families, tsEVA prioritizes **interpretability, robustness, and physical consistency** over marginal gains in statistical fit. This design choice reduces the risk of overfitting, limits the influence of sampling artifacts on dependence identification, and ensures that inferred joint behavior remains directly aligned with the upper-tail processes relevant for hazard and impact assessment.







## `minPeakDistanceInDaysMonovarSampling`
*(Univariate declustering / event independence)*

**Meaning**  
Minimum temporal separation (in days) between two peaks of the **same variable**
for them to be treated as independent extreme events.

**General principle**  
This parameter should represent the typical **duration of an event episode**
for the considered variable, or a longer timescale if justified by the impact
definition.

### Typical guidance by hazard type

- **Coastal hazards (storm surge, wave height, wave runup, total water level)**  
  A value of approximately **3 days** is often appropriate and corresponds to
  a typical storm duration.  
  This avoids counting multiple peaks generated by the same storm as independent
  events.

- **River discharge / runoff**  
  Strongly dependent on **catchment size and hydrologic response**:
  - Small or flashy catchments may justify shorter values.
  - Large catchments may require **several days**, as flood waves can be long-lasting.

- **Heatwaves and droughts**  
  These are intrinsically persistent phenomena. Much **longer values**
  (multi-day to multi-week) may be required to ensure true independence of events.

### Impact-based perspective

From an impact or vulnerability perspective, it can be acceptable—and sometimes
preferable—to set `minPeakDistanceInDaysMonovarSampling` **larger than the minimum
physical duration**, if independence is defined in terms of **recovery time**,
system fatigue, or societal response rather than pure process physics.

---

## `maxPeakDistanceInDaysMultivarSampling`
*(Multivariate pairing window / compound-event definition)*

**Meaning**  
Maximum allowed temporal separation (in days) between peaks of **different
variables** for them to be considered part of the same compound event.

This parameter defines what is meant by a *joint* or *compound* event.

---

### Storm-episode (process-driven) compound events

**Soft default guideline**

maxPeakDistanceInDaysMultivarSampling ≥ max(minPeakDistanceInDaysMonovarSampling)

**Rationale**

If the multivariate pairing window is shorter than the episode duration of at
least one variable, there is a risk of pairing **sub-peaks** rather than the
dominant peaks of each variable within the same event episode.

Ensuring that the pairing window spans at least the univariate episode length
reduces the likelihood of this failure mode and favors pairing of physically
meaningful event maxima.

---

### Impact-based compound events (recovery / cascading impacts)

**Soft guideline**

maxPeakDistanceInDaysMultivarSampling ≫ max(minPeakDistanceInDaysMonovarSampling)

**Rationale**

When compound events are defined from an **impact perspective**, the relevant
timescale may be the **recovery or reduced-capacity window** rather than the
physical simultaneity of drivers.

Examples include:
- erosion or infrastructure weakening followed by high water levels days later,
- storm surge followed by river flooding,
- cascading or sequential hazards affecting the same exposed system.

In these cases, a much larger multivariate pairing window may be appropriate,
provided that the definition is clearly stated and justified.

---

### Final remarks

- There is **no single universally correct choice** for these parameters.
- Their selection must be consistent with the **scientific question**, the
  **physical processes**, and the **impact framing** of the analysis.
- Sensitivity tests using multiple values of
  `maxPeakDistanceInDaysMultivarSampling` are strongly recommended to verify the
  robustness of inferred dependence structures and joint return levels.