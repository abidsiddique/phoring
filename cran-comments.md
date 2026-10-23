## Resubmission Comments for phoring 0.1.0

This is a resubmission of package 'phoring' following incoming pre-test feedback from CRAN auto-check service (Prof. Uwe Ligges).

### Pre-test feedback resolution

1. **LaTeX & PDF Manual Build Errors (Unicode ^^H and ^^G), HTML Math Rendering Note, and Left-over `phoring-manual.tex`**:
   - Root cause: Stray unescaped control characters (`0x07`, `0x08`, `0x09`) in math formulas (`\eqn{\hat{\alpha}_i}`, `\eqn{\exp(\beta) \pm 95\%}`, and `\eqn{\hat{\beta}(\tau)}`) and unescaped `%` characters in `.Rd` descriptions caused `Rd2pdf` and HTML math rendering to fail, leaving `phoring-manual.tex` in the check directory.
   - Fix: Corrected all LaTeX math and `%` escaping in `R/plots.R` and regenerated all `.Rd` documentation files in `man/`. Verified that all 67 `.Rd` files parse and compile cleanly to PDF and HTML with zero errors, warnings, or left-over `.tex` files.

2. **Spell-check NOTE on Econometric Terms in DESCRIPTION**:
   - Flagged terms (`Engle`, `GLS`, `Granger`, `HAC`, `KPSS`, `Koyck`, `Perron`, `Prais`, `SLS`, `Winsten`, `cointegration`, `heteroskedasticity`, `tobit`) are standard econometric method names, statistical acronyms, and author surnames (Engle, Granger, Koyck, Perron, Prais, Winsten).
   - Added `.aspell/defaults.R` and `.aspell/phoring.rds` (along with `inst/WORDLIST`) registering these domain-specific econometric terms.

### Test environments
* local macOS Sequoia 15.6 (aarch64-apple-darwin20), R 4.5.2
* win-builder (devel and release)
* Debian (r-devel-linux-x86_64-debian-gcc)

### R CMD check results
0 ERRORs | 0 WARNINGs | 0 unexpected NOTEs.

### Notes on Package Architecture
* 'phoring' is implemented strictly using base R ('stats', 'graphics', 'grDevices') with zero external package dependencies, ensuring lightweight, high-performance portability.
* Complete academic publication-grade S3 visualization suite with serif typography and multi-panel diagnostic support for all econometric models and tests.
