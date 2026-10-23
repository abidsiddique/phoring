## Resubmission Comments for phoring 0.1.0

This is a resubmission of package 'phoring' following incoming pre-test feedback from CRAN auto-check service (Prof. Uwe Ligges).

### Pre-test feedback resolution

1. **LaTeX & PDF Manual Build Errors (Unicode ^^H and ^^G)**:
   - Root cause: Stray unescaped control characters in roxygen math formulas (`\eqn{\hat{\alpha}_i}`, `\eqn{\exp(\beta) \pm 95\%}`, and `\eqn{\hat{\beta}(\tau)}`) generated non-printable ASCII bytes (0x07 and 0x08).
   - Fix: Corrected escaping in `R/plots.R` and regenerated all Rd files (`plot.econ_panel.Rd`, `plot.econ_count.Rd`, `plot.econ_rdd.Rd`, and `plot.econ_qreg.Rd`). Verified that all Rd files convert cleanly to LaTeX and HTML math without errors or warnings.

2. **Spell-check NOTE on Econometric Terms in DESCRIPTION**:
   - Added `inst/WORDLIST` containing domain-specific econometric terms and author names (`cointegration`, `heteroskedasticity`, `tobit`, `Engle`, `GLS`, `Granger`, `HAC`, `KPSS`, `Koyck`, `Perron`, `Prais`, `SLS`, `Winsten`).

### Test environments
* local macOS Sequoia 15.6 (aarch64), R 4.5.2
* win-builder (devel and release)
* Debian (r-devel-linux-x86_64-debian-gcc)

### R CMD check results
There were 0 ERRORs, 0 WARNINGs, 0 unexpected NOTEs.

### Notes on Package Architecture
* 'phoring' is implemented strictly using base R ('stats', 'graphics', 'grDevices', 'utils') with zero external package dependencies, ensuring lightweight, high-performance portability.
* Complete academic publication-grade S3 visualization suite with serif typography and multi-panel diagnostic support for all econometric models and tests.
* All statistical algorithms and estimation methods are aligned with established econometric theory and peer-reviewed literature.
