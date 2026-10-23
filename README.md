# phoring

> **A Comprehensive Econometric Estimation and Diagnostic Suite for R**  
> *Developed by **Md. Abid Siddique Emon**.*  
> *Unified formula-based estimation, robust covariance structures, automated specification testing, economic interpretation, and publication-ready academic visualizations.*

---

## Overview

In conventional R workflows, estimating an econometric model and conducting comprehensive diagnostic testing requires loading and coordinating multiple separate packages (`sandwich`, `lmtest`, `car`, `plm`, `AER`, `tseries`, `quantreg`):

```r
# Conventional fragmented workflow (requiring multiple packages):
library(sandwich); library(lmtest); library(car); library(tseries); library(plm)
mod <- lm(wage ~ educ + exper + tenure, data = df)
coeftest(mod, vcov = vcovHC(mod, type = "HC1"))
bptest(mod)
vif(mod)
dwtest(mod)
resettest(mod)
```

With **`phoring`**, standard formula syntax delivers full coefficient estimates, robust and clustered standard errors, comprehensive diagnostic test batteries, economic interpretation guidance, and publication-ready tables in a single call:

```r
library(phoring)

# Estimate regression with robust standard errors and automated diagnostics:
reg(wage ~ educ + exper + tenure, data = wage_data)
```

---

## 1. Plain-English Interpretation & Diagnostic Tools

### A. Plain-English Model Interpretation (`interpret`)
Auto-detects functional forms (Log-Log elasticity, Log-Level semi-elasticity, Level-Log, Level-Level) and solves for quadratic turning points $x^* = -\beta_1 / (2\beta_2)$:

```r
data(wage_data)

# Quadratic Mincer equation:
mod <- reg(wage ~ educ + exper + I(exper^2), data = wage_data)
interpret(mod)
# Automatically reports:
# * Shape: Inverted U-shaped (Maximum)
# * Turning point: exper* = 28.4 years (within sample range [1, 45])
# * Insight: Experience increases wage up to 28.4 years, after which diminishing returns dominate.

# Log-Log Constant Elasticity Model:
data(house_prices)
mod_log <- reg(log(price) ~ log(rooms) + crime, data = house_prices)
interpret(mod_log)
# Reports: "A 1% increase in rooms is associated with a 1.25% change in price."
```

### B. Outlier & Influential Observation Detection (`outliers`)
Computes Cook's Distance ($D > 4/n$), Studentized Residuals ($|r| > 2.5$), and Leverage ($h > 2k/n$) to flag points distorting estimates:

```r
outliers(mod, top_n = 5)
```

### C. Structural Break / Parameter Stability Test (`chow_test`)
Tests whether coefficients differ systematically across two subsamples (e.g. before/after policy, male vs. female):

```r
chow_test(wage ~ educ + exper, split_by = "female", data = wage_data)
```

### D. Publication-Grade Academic Visualization Suite (`plot` & `econ_plot`)
Zero external graphics dependencies — generates publication-ready figures with serif typography, L-frame axes (`bty = "l"`), and Oxford navy palette. Every model in the package has dedicated diagnostic visualization supporting individual panel selection via `which =`:

```r
# Standard regression (13 diagnostic views available):
plot(mod)                           # Academic 2x2 grid (Resid vs Fitted, Q-Q, Scale-Loc, Cook's)
plot(mod, which = "leverage")       # Individual panel: Residuals vs Leverage
plot(mod, which = "added_var")      # Partial regression (added-variable) plots
plot(mod, which = "forest")         # Coefficient forest plot with 95% CIs
plot(mod, which = "cusum")          # Parameter stability CUSUM test

# Specialized model-specific plots:
plot(mod_iv)                        # IV/2SLS: First-stage, instrument strength, Hausman
plot(mod_panel)                     # Panel: Entity FE, within-demeaned, residuals
plot(mod_logit)                     # Logit: ROC curve, calibration, confusion matrix
plot(mod_qreg)                      # Quantile regression: Interquantile coefficient trajectories
plot(mod_did)                       # DiD: Parallel trends, 2x2 treatment trajectory, event study
plot(mod_rdd)                       # RDD: Discontinuity scatter, McCrary density sorting
plot(mod_coint)                     # Cointegration: Long-run comovement & ECM equilibrium
plot(mod_garch)                     # GARCH: Volatility clustering & conditional sigma
```


### E. Non-Nested Model Selection: Davidson-MacKinnon J-Test (`j_test`)
Tests competing theoretical specifications against each other when neither model is a subset of the other (Davidson and MacKinnon, 1981):

```r
m1 <- reg(wage ~ educ + exper, data = wage_data)
m2 <- reg(wage ~ tenure + mother_educ, data = wage_data)
j_test(m1, m2)
```

### F. Goldfeld-Quandt Heteroskedasticity Test (`gq_test`)
Orders the sample by a regressor, omits central observations, and tests whether error variance differs across subsamples (Goldfeld and Quandt, 1965):

```r
gq_test(wage ~ educ + exper, split_by = "educ", data = wage_data)
```

---

## 2. Causal Inference & Policy Evaluation

### A. Sharp Regression Discontinuity Design (`rdd`)
Estimates Local Average Treatment Effects (LATE) at policy thresholds using local linear regression with triangular/uniform kernels and McCrary density continuity checks:

```r
# Sharp RDD at cutoff score = 0:
rdd(outcome ~ test_score, data = df, cutoff = 0)
```

### B. Difference-in-Differences (`did`)
2x2 means breakdown and treatment interaction with clustered standard errors:

```r
data(min_wage)
did(fte ~ treated * post, data = min_wage)
```

### C. Instrumental Variables 2SLS (`iv`)
Estimates 2SLS with Weak Instrument First-Stage $F$-statistic ($F > 10$), Wu-Hausman endogeneity test, and Sargan overidentification test:

```r
data(college_distance)
iv(wage ~ exper + tenure + south | educ ~ near_college + mother_educ, data = college_distance)
```

### D. Table 1 Covariate Balance Testing (`desc`)
Generates publication-ready "Table 1" balance checks with normalized differences and two-sample $t$-tests:

```r
desc(wage_data, vars = c("wage", "educ", "exper"), by = "female")
```

---

## 3. Distributional & Quantile Regression (`qreg`)

Estimate effects across the entire conditional distribution (10th, 25th, median 50th, 75th, 90th percentiles) with robust Huber sandwich standard errors:

```r
data(wage_data)

# Median regression (tau = 0.5):
qreg(wage ~ educ + exper, data = wage_data, tau = 0.5)

# Multi-quantile distribution profile:
qreg(wage ~ educ + exper, data = wage_data, tau = c(0.1, 0.25, 0.5, 0.75, 0.9))
```

---

## 4. Complete Panel Data Econometrics (`panel`)

Supports all major panel estimators with entity-clustered standard errors by default:

```r
data(firm_investment)

# 1. One-Way Fixed Effects (Within Estimator)
panel(inv ~ value + capital, data = firm_investment, id = "firm_id", time = "year", model = "fe")

# 2. Two-Way Fixed Effects (Entity + Time Shock FE)
panel(inv ~ value + capital, data = firm_investment, id = "firm_id", time = "year", model = "twfe")

# 3. Random Effects (Swamy-Arora GLS)
panel(inv ~ value + capital, data = firm_investment, id = "firm_id", time = "year", model = "re")

# 4. First Differences (FD)
panel(inv ~ value + capital, data = firm_investment, id = "firm_id", time = "year", model = "fd")

# 5. Between Estimator (Long-Run Group Means)
panel(inv ~ value + capital, data = firm_investment, id = "firm_id", time = "year", model = "between")
```

**Every panel run automatically includes:**
* **Hausman Specification Test:** FE vs. RE test ($p < 0.05 \rightarrow$ Fixed Effects required).
* **Breusch-Pagan LM Test for Panel:** Tests Pooled OLS vs. Random Effects ($p < 0.05 \rightarrow$ Panel effects present).
* **Wooldridge Panel AR(1) Test:** Tests for first-order serial correlation in errors.

---

## 5. Dynamic Models & Time Series Econometrics

```r
data(macro_data)

# 1. Koyck Distributed Lag Model:
# Computes Short-run multiplier, Long-run multiplier, Median lag, and Mean lag:
koyck(inflation ~ unemployment, data = macro_data)

# 2. Durbin's h Test for Dynamic / Autoregressive Models:
# (Standard Durbin-Watson is invalid when lagged Y_{t-1} is a regressor):
m_dyn <- reg(inflation ~ inf_lag + unemployment, data = macro_data)
durbin_h(m_dyn, lag_var = "inf_lag")

# 3. Prais-Winsten & Cochrane-Orcutt FGLS AR(1) Correction:
# Transforms data via rho-differencing to eliminate serial correlation:
prais_winsten(inflation ~ fed_funds + unemployment, data = macro_data)

# 4. Error Correction Model (ECM):
# Estimates short-run adjustments and speed of adjustment (gamma) to equilibrium:
ecm(inflation ~ unemployment + fed_funds, data = macro_data, lags = 1)

# 5. Vector Autoregression (VAR):
# System estimation with AIC/BIC selection and pairwise Granger causality matrix:
var_model(macro_data, vars = c("inflation", "unemployment", "fed_funds"), p = 1)

# 6. Multi-Test Stationarity Battery (ADF + Phillips-Perron + KPSS in 1 table):
unit_root(macro_data$inflation)

# 7. Engle-Granger Two-Step Cointegration Test:
coint(inflation ~ unemployment + fed_funds, data = macro_data)

# 8. Bidirectional Granger Causality Test:
granger(inflation ~ fed_funds, data = macro_data, lags = 2)

# 9. Engle's ARCH LM Test (Volatility Clustering):
m_ts <- reg(inflation ~ fed_funds, data = macro_data, diagnostics = FALSE)
arch_test(m_ts, lags = 2)
```

---

## 6. Complete Microeconometrics & Binary Choice

```r
data(labor_supply)

# 1. Logit with Odds Ratios, ROC-AUC, and Average Marginal Effects (AME):
logit(inlf ~ educ + exper + age + nwifeinc + kidslt6, data = labor_supply)

# 2. Probit with Average Marginal Effects:
probit(inlf ~ educ + exper + age + nwifeinc + kidslt6, data = labor_supply)

# 3. Tobit Censored Regression (with McDonald-Moffitt Decomposition):
tobit(hours ~ educ + exper + age + nwifeinc, data = labor_supply, left = 0)

# 4. Count Data (Poisson Regression with Cameron & Trivedi Overdispersion Test):
count_reg(docvis ~ educ + age + kidslt6, data = labor_supply)
```

---

## 7. Model Comparison & Word Export (`word_table`)

Compile multiple specifications side-by-side with custom labels and export directly to Microsoft Word (`.docx`) or HTML:

```r
m1 <- reg(wage ~ educ, data = wage_data, diagnostics = FALSE)
m2 <- reg(wage ~ educ + exper, data = wage_data, diagnostics = FALSE)
m3 <- reg(wage ~ educ + exper + tenure, data = wage_data, diagnostics = FALSE)

# Export to Word with custom column headers:
word_table(m1, m2, m3, 
           labels = c("Bivariate", "+ Experience", "+ Tenure"),
           file = "wage_models.docx", 
           title = "Table 1: Determinants of Wages")
```

---

## Curated Practice Datasets Included in `phoring`

| Dataset | Category | Observations | What It's Used For |
| :--- | :--- | :---: | :--- |
| **`house_prices`** | Hedonic Pricing / OLS | 506 | Harrison-Rubinfeld housing values, crime rates, pollution (`nox`), VIF & RESET tests. |
| **`wage_data`** | Labor Economics / OLS | 600 | Mincer wage equations, experience curvature, gender wage differential, firm clustering. |
| **`college_distance`** | Instrumental Variables (2SLS) | 800 | Card (1995) returns to schooling using college proximity (`near_college`) as IV. Weak IV & Hausman tests. |
| **`firm_investment`** | Panel Data (Micro) | 200 (40 $\times$ 5) | Grunfeld corporate investment: FE, TWFE, RE, FD, Hausman test, and BP panel test. |
| **`crime_panel`** | Panel Data (Macro/Policy) | 200 (50 $\times$ 4) | County crime rates, police presence, and unemployment with entity-clustered SEs. |
| **`labor_supply`** | Binary & Censored Micro | 753 | Mroz (1987) female labor supply (`inlf`), annual hours (`hours` for Tobit), and doctor visits (`docvis` for Poisson). |
| **`min_wage`** | Causal Inference / DiD | 400 | Card & Krueger (1994) NJ & PA fast food minimum wage hike, 2x2 means breakdown. |
| **`macro_data`** | Time Series Econometrics | 120 quarters | US inflation, unemployment, fed funds rate. `unit_root()`, `coint()`, `granger()`, `arch_test()`, and `ardl()`. |

---

## Installation

You can install the development version of **phoring** directly from GitHub:

```r
# install.packages("remotes")
remotes::install_github("abidsiddique/phoring")
```

Once published on CRAN, install the official release version using:

```r
install.packages("phoring")
```

---

## Author & Citation
* **Author:** Md. Abid Siddique Emon  
* **Email:** [contact.abidsiddique@gmail.com](mailto:contact.abidsiddique@gmail.com)  
* **ORCID:** [0009-0007-6576-027X](https://orcid.org/0009-0007-6576-027X)  
* **GitHub:** [https://github.com/abidsiddique](https://github.com/abidsiddique)  
* **License:** MIT
