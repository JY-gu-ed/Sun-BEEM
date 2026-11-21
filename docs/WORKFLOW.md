# Analysis workflow

The scripts document Sun-BEEM exposure definitions, statistical models and result selection. Paths below are relative to `SUN_BEEM_PROJECT_ROOT`, which holds project inputs and outputs.

## Preparation and main analyses

| Stage | Code | Inputs | Main outputs |
|---|---|---|---|
| Exposure construction | `analysis/step1/1-merge_uv_data.R` to `3-dichotomization.R` | `raw_data_3_8/uv_3_8.csv`, `Dataset_neid_SWR.csv` | Recoded components and Sun-BEEM score in `washed_data/` |
| Supporting variables | Remaining `analysis/step1` scripts | Assessment fields, ancestry, vitamin D and NMR exports | Supporting variables; log(1 + x), standardised NMR measures |
| Covariates | `analysis/step2/var_wash_imputa.R` | Social and alcohol field exports | Recoded and imputed covariates |
| Record linkage | `analysis/step3` | Hospital, cancer and death registration exports; assessment dates | Outcome definitions, dates and prevalent diagnoses |
| Analysis cohorts | `analysis/step4` | Prepared exposures, covariates and linked outcomes | Mortality and incidence analysis tables |
| Main results | `analysis/step5` to `step7` | Analysis tables in `washed_data/` | Cox estimates, trends, sensitivity/stratified models and scenario calculations |
| Protein preparation | `analysis/step8/1-datacombine.R` | Main cohort, `raw_data_3_8/olink.csv`, `step8/daily_uva_temp.csv`, `step8/UVA_UKB_allethn.csv` | `step8/mea_set_merged.rds`, `biomarkers_olink.rds` |
| Protein associations | `analysis/step8/2-lm.R` | Prepared protein object and biomarker names | Adjusted linear associations and FDR in `step8/results/` |
| Protein mediation | `analysis/step8/3-cancer-meadiation.R`, `4-cvd-meadiation.R` | Prepared protein object and biomarker names | Coefficient-change estimates, bootstrap intervals and run settings |
| Protein integration | `analysis/step8/5-cancer-combine.R`, `6-cvd-combine.R` | Protein associations and selected bootstrap run | Joined estimates and original significance filters |
| Protein pathways | `analysis/step8/11-pathwaybackground.R` | Protein associations and Olink annotations used in the script | Enrichment results and displays |
| NMR associations | `analysis/step10/1-combine.R`, `2-lm.R` | Prepared NMR variables and clinical cohort | NMR analysis object and adjusted associations |
| NMR annotation | `analysis/step10/3-washratio.R`, `4-removeratio.R` | Nightingale field/group dictionary and association estimates | Non-ratio annotated results |
| NMR mediator selection | `analysis/step10/10-supplement.R` | Annotated NMR mediator results named in the script | Protective-mediator supplementary tables |

`analysis/step8/9-overlap.R` compares previously saved CVD and cancer Top15 lists. It follows preparation of those lists in the figure workflow.

Input exports use UK Biobank field names as shown in each script. Outdoor time uses fields 1050/1060, protection uses 2267, solarium use uses 2277, and expanded adjustment uses 738, 6142, 20277 and 816. Residential SWR is supplied in `swr2008_kJ`. The scripts retain the original outcome codes, cutpoints, covariates and transformations.

## Component and revision analyses

Run these stages in order when their inputs are available:

1. `analysis/revision/00_component_models.R`: saved cohort preparation and Cox loops for single components, mutually adjusted components, leave-one-component-out scores and score specifications. Writes `elife_revision/results/targeted_cox_results.csv` and the corrected follow-up object used by the checking script.
2. `analysis/revision/01_revision_analysis.R`: checks the baseline against that output, reuses component estimates and adds continuous/quartile, incidence, selected-biomarker and extended-adjustment analyses. Writes to `JID_R1_20261104/results/`.
3. `analysis/revision/01b_extra_confounder_recheck.R`: subsequent employment-coding and extended-adjustment check, with its original outputs retained.
4. `analysis/revision/02_omics_point_checks.R`: saved checking script for protein and metabolite coefficient-change point estimates, using the omic objects and corrected follow-up.

`00_component_models.R` contains the component-model portion of the saved script through the `targeted_cox_results.csv` export. Its later exploratory scenario blocks are separate from the publication's scenario calculations in `step5`.

The Methods-oriented NMR implementation uses `analysis/step10/5-cancer-mediation.R` and `6-cvd-mediation.R`, followed by `7-cancer-combine.R` and `8-cvd-combine.R`. These implement the stated Cox coefficient-change models and participant bootstrap with an unadjusted Sun-BEEM–metabolite association P<0.01 filter. Non-ratio selection uses the existing Nightingale dictionary. `10-supplement.R` exports the corresponding protective-mediator tables.

The added NMR runners and adjusted construct-validation model implement the final Methods.

## Bootstrap runs

Set `SUN_BEEM_N_BOOT` to the requested count (integer, at least 10) and `SUN_BEEM_SEED` to the run seed. The default count is 1,000 for protein and metabolite mediation, matching the final Methods. Seed 123 is the default for this release's execution controls. Both layers use the same covariates, coefficient-change estimator and percentile intervals specified in the Methods.

Protein outputs use `step8/results/bootstrap_<count>_seed_<seed>/`; NMR outputs use the same naming under `step10/results/`. Each run saves requested and successful resampling counts, input fingerprints, settings and `sessionInfo()`. Integration scripts select the same count/seed directory. Fixed per-biomarker random streams keep resampling consistent when a run is resumed. Different counts and seeds use separate directories.

## Environment

`requirements-R.txt` lists the R packages used. `clusterProfiler`, `enrichplot` and `org.Hs.eg.db` are Bioconductor packages; `grid` is included with R. The installation helper installs missing CRAN/Bioconductor packages. Record the environment used with `sessionInfo()` in the project output directory. Source syntax was checked with R 4.5.3.

## Methods alignment and required metadata

The primary mortality preparation defines all-cause deaths from registered death dates, including deaths without a primary cause, and propagates this definition to clinical and molecular cohorts. Administrative mortality censoring remains 4 August 2025.

Incident models end follow-up at the first qualifying event, death, recorded loss to follow-up or the relevant registry cutoff. Set `SUN_BEEM_HOSPITAL_END_DATE` and `SUN_BEEM_CANCER_END_DATE` from the project linkage metadata. Incidence cohorts carry `death_date` and `loss_to_followup_date` (field 191) from preparation. Run the step3/step4 stages before the incidence models.

`analysis/step5/13-construct-validation.R` fits the category comparisons adjusted for age group, sex, BMI category, season and assessment centre, and the unadjusted component comparisons. It reads baseline centre field `p54_i0` from the request/social export, or from the CSV specified by `SUN_BEEM_BASELINE_METADATA`. Season is derived from baseline assessment month: winter December–February, spring March–May, summer June–August and autumn September–November. Missing adjustment values are excluded per response. These covariates are retained in the saved model description.

`analysis/step5/14-complete-case-sensitivity.R` uses pre-imputation core-covariate membership saved by step2 and repeats adjusted mortality models. The PWAS direction/enrichment input threshold is FDR<0.01; the enrichment-term threshold remains FDR<0.05.

Use a new `SUN_BEEM_PROJECT_ROOT` for executing this Methods-oriented release. Source inputs and generated outputs remain in that project directory. This keeps existing project results intact.
