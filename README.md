# Sun-BEEM

Analysis code for *Long-term sun exposure is associated with reduced all-cause mortality and distinct immune–metabolic signatures* (JID-2026-1128).

## Repository contents

This repository contains source code and documentation only. No participant-level data, derived datasets, model objects, numerical result files, figures or execution logs are included.

- `analysis/step1`–`step4`: exposure construction, covariates, record linkage and analysis-cohort preparation.
- `analysis/step5`–`step7`: clinical outcomes, Cox models, scenario calculations and sensitivity analyses.
- `analysis/step8`: proteomic association, mediation and pathway analysis.
- `analysis/step10`: NMR preprocessing and association analysis.
- `analysis/revision`: component, score-sensitivity, selected-biomarker and expanded-adjustment analyses.
- `tools/check_code_only.py`: checks the repository for prohibited data files and sensitive content.

The numbered directories follow the original project workflow. Input tables and intermediate objects referenced by the scripts must be generated or obtained within an authorised analysis environment; they are not bundled here. No automatic end-to-end execution is enabled.

## Data access

UK Biobank data require approval from UK Biobank and are not redistributed by the authors. See https://www.ukbiobank.ac.uk/use-our-data/ for access information. Field names, outcome definitions, thresholds and model parameters in the code describe the analytical methods; they are not participant records.

## Environment and execution

Install the R packages listed in `requirements-R.txt` in the authorised environment. This package list is not a version-locked environment specification.

Set `SUN_BEEM_PROJECT_ROOT` to the project directory containing the authorised input and working directories. Set `SUN_BEEM_EXTERNAL_ROOT` if a script requires separately managed inputs. `.env.example` shows generic placeholders only; never commit a populated `.env` or `.Renviron` file.

Scripts use the original project-relative input and output conventions. Run the required preparation stages before downstream analyses and inspect the input-file references in each script. Some analyses depend on saved intermediate objects from earlier stages. Generated outputs must remain outside this repository.

Start each R process from the authorised project directory so that `here()` resolves the same project-relative inputs:

```sh
cd /path/to/authorised/project
Rscript /path/to/Sun-BEEM/analysis/step1/1-merge_uv_data.R
```

### Bootstrap configuration

The number of bootstrap resamples in the protein mediation scripts can be set as needed with `SUN_BEEM_N_BOOT`. For example, set `Sys.setenv(SUN_BEEM_N_BOOT = "1000")` in R before sourcing a mediation script. The value must be an integer of at least 2. If unset, the scripts retain the archived source setting of 50. The requested count is printed when the script runs; retain it with the run record. This configurable release does not establish the count used for previously generated manuscript results.

Before committing, run:

```sh
python tools/check_code_only.py
```

The check reads source files only; it does not run the UK Biobank analyses.
