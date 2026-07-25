# DAP-2024-SurvivalAnalysis

Analysis code for the manuscript **"Factors Associated with Longevity in
Companion Dogs: Initial Findings from the Dog Aging Project"** (Chiu et al.).
The manuscript, supplemental figures/tables, and supplemental separate files
are in `Manuscript Files/`.

The full pipeline is orchestrated by **`00.RunAll.R`**, which sources the
scripts below in order. Each stage runs from a clean workspace
(`rm(list=ls())` between calls). Inputs live in `data/`, derived outputs in
`results/`, and plots in `figures/`.

## Requirements

- R 4.5.1 (as reported in the manuscript)
- CRAN packages: `survival` (v3.8-3), `flexsurv`, `survminer`, `lubridate`,
  `dplyr`, `tidyverse`, `ggpubr`, `viridisLite`, `stringr`, `haven`,
  `usdata`, `wCorr`, `emmeans`

## Input data (in `data/`)

- `SurvivalAnalysisCuratedDogs_thru_2024_ran_2025_12_16.csv` — DAP curated
  survival release with mortality follow-up through 2024
- `DAP_2024_DogOverview_v1.0.csv` — dog-level demographic overview
- `DAP_2024_CODEBOOK_v1.0.csv` — HLES codebook used to drive the SWAS
- HLES data files loaded by the SWAS scripts
- State-level human mortality / life expectancy files (HDPulse) used by
  `4.2.Geography_DAP_vs_Human.R`

## Pipeline (as sourced by `00.RunAll.R`)

| # | Script | Purpose | Manuscript element |
|---|---|---|---|
| 1 | `0.process_DAP_datafiles.R` | Load and clean the curated release; parse dates; compute `first.age`, `last.age`, `event`; factorize Size × Breed_Class × Sex; save `SurvivalData.RData`. | Cohort construction (N = 41,047) |
| 2 | `1.Cohort_Descriptive_Stats.R` | Crude mortality rates and entry/follow-up/death-age quantiles by Breed_Class × Size × Sex strata; prop tests and stratum ANOVA. | **Table 1**; supplemental AOV tables |
| 3 | `2.1.Survival_DAP_demographics.R` | Primary Cox PH: `Surv(first.age, last.age, event) ~ Size + Breed_Class + Sex`; K-M by 20 strata; median/IQR lifespans; stratified size/breed/sex effect p-values. | **Figure 1** (K-M + forest); **Table 2** |
| 4 | `2.2.Survival_DAP_demographics.Alt.R` | Sensitivity analysis — alternative parameterization (10 kg weight bins). | Supplemental sensitivity results |
| 5 | `2.3.Survival_DAP_demographics.Alt2.R` | Sensitivity analysis — alternative reference strata / model form. | Supplemental sensitivity results |
| 6 | `3.Survival_CommonBreeds.R` | Repeats the demographic survival analysis on the 16 most common single AKC breeds (n > 250). | **Table 3**; breed-specific K-M plots |
| 7 | `4.1.Geography_DAP.R` | Cox stratified by Size × Breed × Sex, adjusting for rural/suburban/urban and state; `emmeans` with `method = "eff"` to obtain each state's deviation from the count-weighted grand mean (no reference state pinned to zero). Saves `Geoeffect-Cox-results.Rdata`. | Geographic results text |
| 8 | `4.2.Geography_DAP_vs_Human.R` | Joins state-level dog HRs to human age-adjusted mortality and life expectancy (HDPulse); inverse-variance-weighted regression plus `wCorr` weighted Pearson/Spearman. | **Figure 2** (r = 0.43, ρ = 0.63, p = 0.0017) |
| 9 | `5.1.HLES_SWAS_analysis.R` | Survey-wide association study over the 804 HLES codebook variables (grouped by module: `dd, oc, pa, de, db, df, dt, mp, hs`, …). Univariate Cox stratified by Size × Breed × Sex; modal level as reference for categorical variables; Benjamini–Yekutieli FDR control. | Models underlying Figs 3–4 |
| 10 | `5.2.0.HLES_SWAS_figures.R` | Manhattan-style −log10(q) plot by module and forest plots of selected hits. | **Figure 3**, **Figure 4**, **Table 4** |
| 11 | `5.2.1.HLES_SWAS_figures_SensAn_MatureAdult.R` | SWAS sensitivity analysis restricted to the "Mature Adult" lifestage (n = 20,505). | Supplemental SWAS sensitivity |
| 12 | `5.2.2.HLES_SWAS_figures_SensAn_CommonBreeds.R` | SWAS sensitivity analysis restricted to the 16 most common single breeds (reduced genetic heterogeneity). | Supplemental SWAS sensitivity |

### Helper

- `ggforest2.R` — custom `survminer::ggforest` variant, `source()`-ed by
  `4.1.Geography_DAP.R`, `5.1.HLES_SWAS_analysis.R`, and the SWAS figure
  scripts. Required by the pipeline; not sourced directly from `00.RunAll.R`.

## Outputs

- `results/` — every CSV/RData artifact named after the script that produced
  it, so any manuscript number can be traced back to a specific script
  (e.g. `Survival_DAP_demographics.SumStats.csv`,
  `Geoeffect-Cox-Human.MR.csv`, `Supp.HLES_Cox_signif_results.csv`).
- `figures/` — PDF/SVG panels used in the main and supplemental figures.

## To reproduce

```r
setwd("<repo root>")
source("00.RunAll.R")
```

## Files NOT part of the manuscript pipeline

The following files exist in the repository but are **not** sourced by
`00.RunAll.R` and are not required to reproduce the manuscript:

- `PrelimDataFigures.R` and the associated `PrelimDataFig.*.pdf` / `.svg`
  (ad-hoc preview figures)
- `Crude-Rates-By-Cohort-For-Daniel.R` (one-off external request)
- `X.HR_interp.R` (interactive hazard-ratio interpretation utility)
