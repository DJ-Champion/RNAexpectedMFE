# RNAexpected

Expected RNA minimum free energy (MFE) from GC content and length alone.

RNAfold MFE scales strongly with sequence length and GC content, so a raw MFE
says little about whether a particular mRNA is unusually structured. This
project fits a simple model of the MFE that RNAfold gives, on average, for a
*random* sequence of a given GC content and length. An observed MFE can then be
compared against that baseline (for example `observed - expected`).

The recommended model (C) is

```
expected MFE = (a · GC^b + c) · length + d,   capped at 0 kcal/mol
a = -0.8392021   b = 2.3535788   c = -0.1517445   d = 9.3178327
```

with GC as a fraction (0–1) and length in nt. It was fitted to 406,470
synthetic sequences spanning 20–4000 nt and 0–100% GC, by nonlinear least
squares weighted by 1/length (MFE noise grows with length). Its systematic
error is about 1–3.5% for sequences of 100–4000 nt, but it predicts too little
structure for very short sequences (about −13% at 50–100 nt), and lengths
above 4000 nt are extrapolated. See [docs/methods.md](docs/methods.md) for the
full method, model comparison and limitations.

**Coefficient change (November 2026):** switching from an unweighted to a
1/length-weighted fit changed model C's coefficients (mainly d: 13.56 → 9.32).
Projects that hard-code the older values (`mRNA_stability_R`, `RNAstab`)
need updating.

## Using the R package

Install from GitHub (the repository is `RNAexpectedMFE`; the package is
`RNAexpected`):

```r
# install.packages("remotes")
remotes::install_github("DJ-Champion/RNAexpectedMFE")
```

or from a local clone with
`install.packages("path/to/RNAexpectedMFE", repos = NULL, type = "source")`.

```r
library(RNAexpected)

expected_mfe(gc = 0.45, length = 1500)
#> [1] -410.5045
expected_mfe(gc = c(0.3, 0.5, 0.7), length = 1000)
#> [1] -191.7702 -306.6250 -504.9138
expected_mfe(0.45, 1500, model = "A")   # alternative models A, B, C
expected_mfe_coefficients("C")
```

- `gc` must be a **fraction**. Values above 1 are an error, because they
  almost always mean percentages were passed. Divide by 100 first.
- Lengths outside 20–4000 nt still return a value, with a message that it is
  an extrapolation.
- Predictions are capped at 0 kcal/mol: very short or GC-poor sequences can't
  have a positive MFE.

The package has no dependencies beyond base R.

## Reproducing the analysis

Everything in `results/` and `R/coefficients.R` is rebuilt from the raw
RNAfold output in `data-raw/` by:

```bash
Rscript analysis/run_all.R
```

This takes about 2 minutes and needs `dplyr`, `ggplot2`, `patchwork`, `purrr`,
`readr`, `stringr`, `tibble` and `tidyr`. Run it from the project root, or open
`RNAexpected.Rproj` and source the numbered scripts in order.

| Script | Does | Writes |
|---|---|---|
| `analysis/01_fit_models.R` | Fits models A, B and C (weights 1/length) to all sequences, and to replicates 1–8 for a held-out check | `R/coefficients.R`, `results/tables/model_coefficients.csv`, `results/models/fits.rds` |
| `analysis/02_compare_models.R` | Fit statistics, held-out error, error by length band, residual diagnostics | `results/tables/model_comparison.csv`, `results/tables/rmse_by_length.csv`, `results/figures/rmse_by_length.*`, `results/figures/diagnostics/model_*.*` |
| `analysis/03_figures.R` | Key figures: mean MFE vs length (by GC) and mean MFE vs GC (by length), one version per model in `figure_models` (default C); model error vs length with human region lengths | `results/figures/mfe_vs_length_model_C.*`, `results/figures/mfe_vs_gc_model_C.*`, `results/figures/error_vs_length.*`, `results/tables/rel_error_by_length.csv` |

`R/coefficients.R` is generated: don't edit it by hand. If a refit changes the
coefficients, the package test in `tests/` fails on purpose, because downstream
projects (`mRNA_stability_R`, `RNAstab`) hard-code the current model C values.

## Layout

```
R/                 package code: model definitions, generated coefficients, expected_mfe()
man/, tests/       package help pages and tests (R CMD check)
analysis/          the fitting/comparison pipeline; analysis/R/ holds its helpers
data-raw/          frozen inputs: RNAfold MFEs (*_rnafold_results.csv) and measured
                   GC (*_actualGC.txt), split into length bands a–h; plus
                   human_region_summary.csv (real transcript lengths/GC, for context)
results/           tables and figures produced by the pipeline
docs/methods.md    methods, model comparison, limitations
```

To add a candidate model, add an entry to `mfe_models` in `R/models.R`
(formula plus start values) and rerun the pipeline. It is fitted, compared,
and exposed through `expected_mfe(model = ...)` automatically.

## Checking the package

```bash
R CMD build . && R CMD check RNAexpected_0.1.0.tar.gz
```

## Data provenance

The synthetic sequences were folded with RNAfold (ViennaRNA) using default
settings. The script that generated the sequences, and the exact RNAfold
version, were not kept. The RNAfold output in `data-raw/` is therefore treated
as a frozen input. See the provenance notes in [docs/methods.md](docs/methods.md).
