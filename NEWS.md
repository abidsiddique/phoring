# phoring 0.1.0

* Initial release of `phoring`: A unified toolkit for econometric estimation, automated specification testing, and publication-ready diagnostics.
* **Core OLS & Diagnostics:** Automated heteroskedasticity (Breusch-Pagan, White), autocorrelation (Durbin-Watson, Breusch-Godfrey), multicollinearity (VIF), functional form (Ramsey RESET), and normality (Jarque-Bera) tests.
* **Standard Errors:** Robust (HC0, HC1, HC2, HC3), cluster-robust ($G-1$ degrees of freedom), and Newey-West HAC standard errors, with full support for Weighted Least Squares (WLS).
* **Rank Deficiency & Collinearity Protection:** Informative diagnostic guidance identifying exact linear dependencies and collinear predictors in dummy variable traps.
* **Instrumental Variables (2SLS):** Automated weak instruments first-stage F test, Wu-Hausman endogeneity test, and Sargan overidentification test with cluster-robust SEs and multi-endogenous regressor support.
* **Panel Data:** Within Fixed Effects, Two-Way Fixed Effects for both balanced and unbalanced panels, Swamy-Arora Random Effects GLS, First Differences, Between Estimator, Hausman specification test, and Wooldridge panel autocorrelation test.
* **Microeconometrics:** Logit and Probit models with Average Marginal Effects (AME) and delta-method standard errors; Tobit censored regression with McDonald-Moffitt decomposition; Poisson and true Negative Binomial regression with Cameron-Trivedi overdispersion test and Likelihood Ratio validation.
* **Causal Inference:** Difference-in-Differences (DiD) with 2x2 means breakdown and covariate adjustments; Sharp Regression Discontinuity Design (RDD) with triangular/uniform kernels and McCrary density continuity check.
* **Quantile Regression:** Koenker-Bassett quantile estimation across single or multiple quantiles with Huber sandwich covariance.
* **Time Series & Dynamics:** Augmented Dickey-Fuller (trend, drift, none), Phillips-Perron, and KPSS stationarity tests with continuous MacKinnon response surface p-values; Engle-Granger cointegration test; bidirectional Granger causality; Engle's ARCH test; ARDL(p, q); Koyck distributed lag with Liviatan IV option; Prais-Winsten AR(1) correction; Error Correction Models (ECM); and Vector Autoregression (VAR).
* **Reporting & Export:** `word_table()` multi-model side-by-side comparison tables exported to Word and HTML, `interpret()` plain-English economic guidance, and a comprehensive academic `plot()` suite with serif typography, L-frame axes, and model-specific diagnostic graphs across all 28 econometric model and test classes.
* **S3 Ecosystem:** Universal S3 method coverage (`coef`, `vcov`, `residuals`, `fitted`, `nobs`, `confint`, `predict`, `logLik`, `AIC`, `BIC`, `plot`).
