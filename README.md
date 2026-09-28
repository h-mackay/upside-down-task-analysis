# upside-down-task-analysis

Eye-tracking and behavioural analysis pipeline for a study investigating object affordances and canonicality in visual feature processing.

## Structure

```
analysis/           numbered stage scripts (01-04)
R/                  helper functions
data/raw/           raw data (not included; see Data)
data/derivatives/   built data (git ignored)
outputs/            tables, figures, diagnostics (git ignored)
config.yml          paths, factors, exclusions, seed
```

## Requirements

- R 4.6.1
- Package versions are recorded in `renv.lock`; `renv::restore()` installs them

## Data

Participant data are not included in this repository.

The scripts require two folders in `data/raw/` (paths set in `config.yml`):

```
data/raw/et_processed/   Tobii exports; P<id>_<Orientation|Saturation>_TB_processed.xlsx
data/raw/eprime_clean/   E-Prime files; P<id>_<Orientation|Saturation>_TB_cleaned.xlsx
```

## Reproducing

1. Clone the repository
2. Open `upside-down-task-analysis.Rproj` in RStudio, or start R in the project folder (activates project environment)
3. Run `renv::restore()` to install the package versions from the lockfile
4. Put the data in `data/raw/` (see Data)
5. Run the scripts in `analysis/` in numbered order:

| Script | Writes |
|---|---|
| `01_data.R` | `data/derivatives/` (`trial_data.rds`, `et_data.rds`, `behav_data.rds`); starts `outputs/run_log.txt` |
| `02_anova.R` | `outputs/tables/et/`, `outputs/tables/behavioural/`, `outputs/diagnostics/` |
| `03_figures.R` | `outputs/figures/et/`, `outputs/figures/behavioural/` |
| `04_correlations.R` | `outputs/tables/rmcorr/`, `outputs/figures/rmcorr/` |
