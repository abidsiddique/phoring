#' Cross-Sectional and Panel Wage Determination Dataset
#'
#' A simulated labor economics dataset modeled after the Mincerian wage equation
#' and empirical labor economics literature. Includes individual characteristics,
#' parental education instruments, firm identifiers for clustering, and panel time periods.
#'
#' @format A data frame with 600 observations and 11 variables:
#' \describe{
#'   \item{wage}{Hourly wage in US dollars}
#'   \item{educ}{Years of education}
#'   \item{exper}{Years of labor market experience}
#'   \item{tenure}{Years with current employer}
#'   \item{female}{Indicator variable: 1 if female, 0 if male}
#'   \item{married}{Indicator variable: 1 if married, 0 otherwise}
#'   \item{mother_educ}{Mother's years of education (instrument for educ)}
#'   \item{father_educ}{Father's years of education (instrument for educ)}
#'   \item{firm_id}{Firm identifier for clustering and panel analysis}
#'   \item{year}{Year of observation (2020-2023)}
#'   \item{inlf}{Indicator variable: 1 if in labor force / high wage, 0 otherwise}
#' }
#' @source Simulated based on standard labor econometrics specifications.
"wage_data"

#' Hedonic Housing Prices and Urban Environmental Quality
#'
#' A cross-sectional dataset of 506 housing tracts modeled after the classic
#' Boston housing study (Harrison & Rubinfeld, 1978). Ideal for multiple regression,
#' heteroskedasticity testing (Breusch-Pagan, White), and Ramsey RESET functional form tests.
#'
#' @format A data frame with 506 rows and 7 variables:
#' \describe{
#'   \item{price}{Median value of owner-occupied homes in $1,000s}
#'   \item{rooms}{Average number of rooms per dwelling}
#'   \item{crime}{Per capita crime rate by town}
#'   \item{nox}{Nitric oxides concentration (parts per 10 million)}
#'   \item{dist}{Weighted distances to five employment centres}
#'   \item{tax}{Full-value property-tax rate per $10,000}
#'   \item{ptratio}{Pupil-teacher ratio by town school district}
#' }
#' @source Modeled after Harrison and Rubinfeld (1978) hedonic housing study.
"house_prices"

#' Returns to Schooling and Proximity to College Dataset
#'
#' A dataset of 800 individuals designed for Instrumental Variables (2SLS) estimation,
#' based on David Card's (1995) landmark study using geographic proximity to a 4-year
#' college as an exogenous instrument for educational attainment.
#'
#' @format A data frame with 800 rows and 9 variables:
#' \describe{
#'   \item{wage}{Hourly wage in US dollars}
#'   \item{educ}{Years of completed education (endogenous regressor)}
#'   \item{exper}{Years of potential labor market experience}
#'   \item{tenure}{Years with current employer}
#'   \item{near_college}{Indicator: 1 if grew up near a 4-year college (excluded instrument)}
#'   \item{mother_educ}{Mother's years of education (excluded instrument)}
#'   \item{father_educ}{Father's years of education (excluded instrument)}
#'   \item{south}{Indicator: 1 if residing in the South}
#'   \item{black}{Indicator: 1 if Black}
#' }
#' @source Modeled after Card (1995), "Using Geographic Variation in College Proximity".
"college_distance"

#' Corporate Investment Panel Dataset (Grunfeld Style)
#'
#' A balanced panel dataset of 40 manufacturing firms observed annually over 5 years (2019-2023),
#' totaling 200 observations. Modeled after the classic Grunfeld investment panel. Ideal
#' for Fixed Effects (Within), Random Effects (GLS), and the Hausman specification test.
#'
#' @format A data frame with 200 rows and 5 variables:
#' \describe{
#'   \item{firm_id}{Firm identifier (1 to 40)}
#'   \item{year}{Year of observation (2019 to 2023)}
#'   \item{inv}{Gross corporate capital investment}
#'   \item{value}{Market value of the firm}
#'   \item{capital}{Real value of plant and equipment capital stock}
#' }
#' @source Modeled after Grunfeld (1958) corporate investment panel.
"firm_investment"

#' County-Level Crime Rates Panel Dataset
#'
#' A balanced panel dataset of 50 counties observed over 4 years (2020-2023), totaling
#' 200 observations. Useful for testing the economic model of crime with Fixed Effects,
#' Random Effects, and entity-clustered standard errors.
#'
#' @format A data frame with 200 rows and 6 variables:
#' \describe{
#'   \item{county_id}{County identifier (1 to 50)}
#'   \item{year}{Year of observation (2020 to 2023)}
#'   \item{crime_rate}{Reported crimes per 1,000 residents}
#'   \item{prbarr}{Probability of arrest (ratio of arrests to reported offenses)}
#'   \item{police}{Police officers per 1,000 residents}
#'   \item{unemp}{County unemployment rate (percentage)}
#' }
#' @source Modeled after Cornwell and Trumbull (1994) panel crime study.
"crime_panel"

#' Female Labor Force Participation Dataset (Mroz Style)
#'
#' A microeconometric cross-section of 753 married women modeled after Mroz (1987).
#' Ideal for binary choice estimation via Logit, Probit, and computing Average Marginal Effects (AME).
#'
#' @format A data frame with 753 rows and 7 variables:
#' \describe{
#'   \item{inlf}{Binary indicator: 1 if woman worked in the labor force, 0 otherwise}
#'   \item{nwifeinc}{Non-wife family household income (in $1,000s)}
#'   \item{educ}{Years of schooling}
#'   \item{exper}{Years of past work experience}
#'   \item{age}{Age in years}
#'   \item{kidslt6}{Number of young children under 6 years old}
#'   \item{kidsge6}{Number of older children between 6 and 18 years old}
#'   \item{hours}{Annual hours worked (0 for non-working women, ideal for Tobit)}
#'   \item{docvis}{Annual doctor visits (integer count, ideal for Poisson regression)}
#' }
#' @source Modeled after Mroz (1987), "The Sensitivity of an Empirical Model of Married Women's Hours of Work".
"labor_supply"

#' State Minimum Wage Policy Evaluation (Card & Krueger DiD)
#'
#' A Difference-in-Differences evaluation dataset of 400 fast-food restaurants in New Jersey
#' (treated state) and eastern Pennsylvania (control state) observed before and after New Jersey's
#' minimum wage increase.
#'
#' @format A data frame with 400 rows and 5 variables:
#' \describe{
#'   \item{store_id}{Store identifier (1 to 400)}
#'   \item{chain}{Fast-food chain name (Burger King, KFC, Roy Rogers, Wendys)}
#'   \item{treated}{State indicator: 1 if New Jersey (treated), 0 if Pennsylvania (control)}
#'   \item{post}{Time period: 1 if after the wage hike, 0 if before}
#'   \item{fte}{Full-time equivalent employment at the restaurant}
#' }
#' @source Modeled after Card and Krueger (1994) minimum wage study.
"min_wage"

#' Quarterly Macroeconomic Time Series (Phillips Curve)
#'
#' A quarterly macroeconomic time series dataset spanning 120 quarters (1994Q1 to 2023Q4).
#' Designed for time series regressions, serial correlation testing (Durbin-Watson, Breusch-Godfrey),
#' and Newey-West HAC standard errors.
#'
#' @format A data frame with 120 rows and 5 variables:
#' \describe{
#'   \item{quarter}{Calendar quarter (e.g., "1994Q1" to "2023Q4")}
#'   \item{inflation}{Annualized CPI inflation rate (percentage)}
#'   \item{unemployment}{Civilian unemployment rate (percentage)}
#'   \item{fed_funds}{Federal funds effective interest rate (percentage)}
#'   \item{gdp_growth}{Annualized real GDP growth rate (percentage)}
#' }
#' @source Simulated based on US macroeconomic time series dynamics.
"macro_data"
