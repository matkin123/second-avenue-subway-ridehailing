# Did the Second Avenue Subway reduce ride-hailing demand?

**A difference-in-differences study of the Q line extension's effect on Uber/Lyft trips in Manhattan, 2016–2017.**

<p align="center">
  <img src="analysis/fig_subway_ridehailing.png" width="820" alt="Subway ridership increased after the Q line opened, while ride-hailing in the treated corridor fell below its counterfactual">
</p>

When the Second Avenue Subway (Q line, Phase 1) opened on **1 January 2017**, it gave the
Upper East Side its first new rapid-transit spine in decades. This project asks whether that
new transit capacity pulled riders away from for-hire ride-hailing (Uber/Lyft) in the
corridor it serves — and estimates the size of the substitution.

**Headline result:** ride-hailing demand in the treated corridor fell by **about 10%**
relative to the rest of Manhattan after the line opened. The preferred estimate is
**−10.6%** (a pre-trend-matched synthetic control); estimators range from **−7% to −13%**.

---

## Design in one paragraph

Monthly ride-hailing trip counts are observed for Manhattan taxi zones across 2016–2017.
Seven zones along the new line (the Upper East Side / Yorkville corridor) are the **treatment**
group; the remaining Manhattan zones are **controls**. A difference-in-differences comparison
isolates the change in treated-zone demand after January 2017 that is *over and above* the
city-wide trend — the latter being substantial, since ride-hailing was growing everywhere at
the time. Because treatment timing is common to all treated zones, the design is immune to the
staggered-adoption pathologies (Goodman-Bacon, de Chaisemartin–D'Haultfœuille) that complicate
many modern DiD studies.

<p align="center">
  <img src="analysis/map_treatment.png" width="340" alt="Treated Upper East Side taxi zones vs. control zones across Manhattan">
</p>

## Headline results

| Estimator | Effect | Notes |
|---|---:|---|
| **Synthetic control, per-zone average (preferred)** | **−10.6%** | near-perfect pre-treatment fit (RMSPE ≈ 0.003) |
| Two-way fixed-effects DiD (61-zone main) | −13.2% | cluster-robust SE 0.055, p < 0.01, ≈ −10k trips/zone/mo |
| TWFE + group-specific linear trend | −9.0% | shrinks but survives |
| Drop month of opening (Dec '16 / Jan '17) | −13.2% / −14.0% | insensitive to the transition month |
| Synthetic control, pooled | −7.0% | conservative floor (plotted below) |
| HonestDiD (relative-magnitudes) | CI stays negative | strict Mbar = 0 bound ≈ [−21%, −6.7%] |

**Honest range: −7% to −13%, central estimate ≈ −10%.**

### The story in one chart

The treated-zone mean and its synthetic counterfactual track almost perfectly through 2016,
then separate after the line opens in January 2017 — the treated corridor grows more slowly
than its synthetic twin.

<p align="center">
  <img src="analysis/scm_aggregate_61.png" width="760" alt="Synthetic control: treated mean vs. synthetic counterfactual">
</p>

### Parallel trends and the event study

Pre-period event-study coefficients sit on zero (pre-trend slope +0.002 log-pts/mo, p = 0.49),
and every post-period month is negative.

<p align="center">
  <img src="analysis/event_study_61.png" width="760" alt="Event study, 61-zone sample">
</p>

Raw ridership by group shows the context: both groups rise together through 2016 (the
ride-hailing boom), with the treated corridor's growth bending down after the opening.

<p align="center">
  <img src="analysis/trajectories_61_raw.png" width="760" alt="Mean ride-hailing trips per zone per month, by group">
</p>

## Robustness suite

Every check lives in a numbered, runnable script under [`analysis/`](analysis/).

| Check | Script | What it probes |
|---|---|---|
| Main DiD (two-way FE, 61 zones) | `05_did_main_61zones.R` | the headline effect |
| Full model suite (OLS / FE / event study / spatial error models) | `06_models_full_suite.R` | specification and spatial-error sensitivity |
| Pre-trends & parallel-trends diagnostics | `07_pretrends.R` | event study, group trends, month drops |
| HonestDiD (Rambachan–Roth 2023) | `08_honest_did.R` | sensitivity to parallel-trends violations (relative-magnitudes & smoothness) |
| Synthetic control (Abadie, via CVXR) | `09_synthetic_control.R` | design-based counterfactual, pooled and per-zone |
| Spatial spillover (donut / Butts 2023 ring / SLX) | `10_spatial_spillover.R` | SUTVA / interference from neighbouring control zones |
| Hausman / FE-vs-RE (Mundlak-Wu) | `11_hausman.R` | estimator choice, with the degenerate-Hausman explained |
| 52-zone contiguous subset | `12_did_robustness_52zones.R` | a tighter, fully-contiguous spatial subsample |

A note on identification: the fixed-effects estimator is correct **by design** here — subway
proximity is correlated with time-invariant zone traits (income, density, land use), so a
random-effects assumption fails. `11_hausman.R` documents why the textbook Hausman test
returns a degenerate p ≈ 1 (all covariates are a single cross-section, so the within estimator
drops them) and replaces it with a balance table, a cluster-robust Mundlak-Wu test (p ≈ 0.07,
leaning FE), and a side-by-side FE/RE comparison (they agree to within 0.5%).

## Repository layout

```
.
├── data/                     raw inputs
│   ├── 2016sum.csv           monthly ride-hailing counts by taxi zone, 2016
│   ├── 2017sum.csv           monthly ride-hailing counts by taxi zone, 2017
│   └── MN.csv                NYC PLUTO land-use records (Manhattan)
├── analysis/
│   ├── 00_map_treatment.R    treatment-vs-control map
│   ├── 01_build_ridehailing.R … 12_did_robustness_52zones.R
│   ├── *.rds                 committed checkpoints (lets you skip the rebuild)
│   ├── *.png                 figures
│   ├── *.csv                 result tables
│   └── geo/                  TLC taxi-zone shapefile
└── docs/
    └── summary.pdf           one-page plain-language summary
```

## Reproducing the results

All scripts are run **from the repository root** and read the committed checkpoints in
`analysis/`, so you can reproduce every model and figure without re-downloading any external
data:

```bash
Rscript analysis/05_did_main_61zones.R      # headline DiD
Rscript analysis/07_pretrends.R             # event study + parallel-trends figures
Rscript analysis/09_synthetic_control.R     # synthetic control
Rscript analysis/10_spatial_spillover.R     # spillover-robust DiD
Rscript analysis/08_honest_did.R            # HonestDiD bounds (needs the HonestDiD package)
```

To rebuild the data pipeline from the raw inputs instead of the checkpoints, set a Census API
key and run the build scripts in order:

```bash
export CENSUS_API_KEY=your_key_here         # census.gov/data/developers
Rscript analysis/01_build_ridehailing.R     # clean trip counts (corrects a TLC data defect)
Rscript analysis/02_build_census.R          # ACS demographics -> taxi zones
Rscript analysis/03_build_landuse.R         # PLUTO land use -> taxi zones
Rscript analysis/04_build_sample_and_weights.R  # 61-zone sample + spatial weights
```

### A data-quality note

The published NYC TLC November-2017 for-hire-vehicle trip file duplicated whole daily batches,
inflating every taxi zone's count by roughly 1.45×. `01_build_ridehailing.R` de-duplicates
November 2017 against the independently-reported dispatch totals and retains the (clean)
December 2017 holiday peak. This correction is applied before any modelling.

## Software

R, with `plm`, `spdep`, `splm`, `sandwich`, `lmtest`, `HonestDiD`, `CVXR`, `sf`, `tidycensus`,
`tigris`, `dplyr`, and `ggplot2`.

---

*A one-page, plain-language summary of this project is in [`docs/summary.pdf`](docs/summary.pdf).*
