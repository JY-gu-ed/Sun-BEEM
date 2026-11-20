# Sun-BEEM

Analysis code for *Long-term sun exposure is associated with reduced all-cause mortality and distinct immune–metabolic signatures* (JID-2026-1128).

This repository presents the Sun-BEEM analysis workflow and statistical logic: exposure construction, clinical outcome models, molecular associations, mediation and sensitivity analyses. It contains code and documentation. UK Biobank inputs are accessed under authorisation and kept in the research environment.

## Code organisation

- `analysis/step1`–`step4`: exposures, covariates, record linkage and cohort preparation.
- `analysis/step5`–`step7`: mortality/incidence models, scenario calculations and sensitivity analyses.
- `analysis/step8`: protein associations, mediation, integration and pathways.
- `analysis/step10`: NMR preparation, associations, annotation and mediator selection.
- `analysis/revision`: component models, score robustness, expanded adjustment and selected-biomarker analyses.
- `tools`: dependency installation, project-directory preparation and the code-only check.

See [the workflow guide](docs/WORKFLOW.md) for execution order, input/output dependencies and bootstrap settings.

## Setup

Install missing R dependencies from the code directory:

```r
source("tools/install_dependencies.R")
```

Set a project directory for inputs and generated outputs, then prepare its folders:

```r
Sys.setenv(SUN_BEEM_PROJECT_ROOT = "/path/to/project")
Sys.setenv(SUN_BEEM_EXTERNAL_ROOT = "/path/to/external/inputs")
source("tools/prepare_project.R")
```

Run the selected stage using its code path, for example:

```r
source("/path/to/Sun-BEEM/analysis/step1/1-merge_uv_data.R")
```

Bootstrap count and seed are configurable:

```r
Sys.setenv(SUN_BEEM_N_BOOT = "1000", SUN_BEEM_SEED = "123")
source("/path/to/Sun-BEEM/analysis/step8/3-cancer-meadiation.R")
```

Each count/seed combination has a separate output directory and saved run settings. The example sets 1,000 resamples; the default count is 50.

## Data access and repository checks

UK Biobank access information: https://www.ukbiobank.ac.uk/use-our-data/

Keep project data, generated results and local run records outside the code directory. Before publishing code changes, run the [UK Biobank Git Audit Tool](https://github.com/UK-Biobank/UKB-Git-Audit-Tool) and the source-file check:

```sh
python tools/check_code_only.py
```
