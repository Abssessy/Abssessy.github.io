// =============================================================================
// This code sample demonstrates proficiency in Stata for econometric analysis,
// covering data management, time-series operations, and panel data techniques.
// All identifying information has been removed using generic variables and
// sample operations to maintain anonymity.
// =============================================================================

// ******************************
// 1. DATA MANAGEMENT AND CLEANING
// ******************************

// Load and examine built-in dataset
sysuse uslifeexp, clear
describe
summarize

// Data import and cleaning demonstration
import delimited using "sample_data.csv", clear  // Generic filename for anonymity

// Date variable processing
gen date = date(trddt, "YMD")
format date %td
drop trddt
rename date trddt

// Convert string variables to numeric
destring price volume, replace

// ******************************
// 2. TIME-SERIES ANALYSIS
// ******************************

// Set time-series structure
tsset trddt

// Calculate returns and moving average
gen ret = D.price / L.price * 100
gen ma_price = (L.price + L2.price + L3.price) / 3  // 3-period moving average

// Statistical analysis
summarize price ret, detail
histogram ret, fraction normal title("Return Distribution")
corrgram ret, lags(12)  // Autocorrelation analysis

// ******************************
// 3. PANEL DATA ANALYSIS
// ******************************

// Simulated panel data example
xtset firm_id year

// Fixed effects model
xtreg investment cashflow growth, fe
estimates store fe

// Random effects model
xtreg investment cashflow growth, re

// Model specification test
hausman fe . 

// Output formatted results
esttab fe, cells(b(star fmt(3)) se(par fmt(2))) ///
    starlevels(* 0.1 ** 0.05 *** 0.01) ///
    stats(N r2, fmt(%9.0g %9.3f)) title("Regression Results") 

// ******************************
// 4. RANDOM WALK ANALYSIS
// ******************************

// Set working directory and load data
cd "D:\data"
use stock_data.dta, clear

// Prepare time series data
sort trddt
gen day = _n
tsset day

// Log transformation of stock prices
gen lp = ln(price)
sktest price lp  // Test for normality

// Autocorrelation analysis
corr lp L.lp L2.lp L3.lp  // Price autocorrelation

// Calculate returns
gen r = lp - L.lp
corr r L.r L2.r L3.r  // Return autocorrelation

// Unit root tests
dfgls price  // DF-GLS test for unit root
reg r L(1/3).r  // AR model for returns

// Partial autocorrelation function
pac r, lags(10)

// Variance ratio test demonstration
egen sd1 = sd(r)  // Daily return volatility

// Five-day returns (assuming 5 trading days per week)
gen r5 = lp - L4.lp
egen sd5 = sd(r5)

// Variance ratio calculation
gen vr = sd5 / (sqrt(5) * sd1)  // Corrected formula

// ******************************
// 5. EVENT STUDY METHODOLOGY
// ******************************

// Event study example using estudy package
ssc install estudy
use event_study_data.dta, clear

// Market model event study with 3-day event window
estudy company1 company2 company3, ///
    datevar(date) evdate(20210615) dateformat(YMD) ///
    lb1(-1) ub1(1) indexlist(mkt) decimal(4)

// Fama-French three-factor model example
estudy company1 company2 company3, ///
    datevar(date) evdate(20210615) dateformat(YMD) ///
    lb1(-1) ub1(1) modtype(MFM) indexlist(mkt smb hml) ///
    diagnosticsstat(KP)

// Alternative event study implementation
eventstudy2 firm_id date using returns_data, ///
    ret(return) car1LB(-1) car1UB(1) mod(FM) ///
    marketfile(market_returns) mar(MKT)

// ******************************
// 6. FAMA-MACBETH REGRESSION
// ******************************

clear all
set more off

// Set working directory and load data
cd "D:\data"
import delimited "fm.csv", clear

// Rename variables for clarity
rename ÿþstkcd stkid
rename trdwnt trading_week
rename wretwd weekly_return
rename wsmvttl market_value
rename wclsprc closing_price

// Generate year variable from date string
gen year = real(substr(trading_week, 1, 4))

// Create market portfolio returns (equal-weighted)
bysort trading_week: egen market_return = mean(weekly_return)

// Save temporary dataset
save temp_data, replace

// Estimate initial betas using 2001 data
statsby beta_coef = _b[market_return], by(stkid) saving(beta_2001, replace): ///
    reg weekly_return market_return if year == 2001

// Portfolio formation based on estimated betas
use beta_2001, clear
drop if missing(beta_coef) | beta_coef < 0  // Remove problematic estimates
xtile portfolio = beta_coef, nq(30)        // Create 30 beta-sorted portfolios
keep stkid portfolio
save portfolios, replace

// Merge portfolio assignments with main dataset
use temp_data, clear
merge m:1 stkid using portfolios, keep(matched) nogenerate

// Calculate portfolio returns
bysort portfolio trading_week: egen portfolio_return = mean(weekly_return)

// Fama-MacBeth cross-sectional regression
preserve
    keep portfolio trading_week portfolio_return beta_coef
    duplicates drop  // Collapse to portfolio-level data

    // Stage 2: Cross-sectional regression for each period
    statsby gamma = _b[beta_coef], by(trading_week) saving(gamma_estimates, replace): ///
        reg portfolio_return beta_coef

// Statistical inference on risk premium
use gamma_estimates, clear
summarize gamma
local t_stat = r(mean) / (r(sd) / sqrt(r(N)))
display "Average risk premium: " r(mean)
display "T-statistic: `t_stat'"

// ******************************
// 7. FAMA-FRENCH THREE-FACTOR MODEL
// ******************************

// Load Fama-French factors
import delimited "FF3_factors.csv", clear
rename date trading_week
rename mktrf market_excess_return
rename smb smb_factor
rename hml hml_factor

// Merge with portfolio returns
merge 1:1 trading_week portfolio using portfolio_returns, nogenerate

// Run time-series regression for each portfolio
statsby alpha = _b[_cons] market_beta = _b[market_excess_return] ///
         smb_beta = _b[smb_factor] hml_beta = _b[hml_factor], ///
         by(portfolio) saving(ff3_results, replace): ///
         reg portfolio_excess_return market_excess_return smb_factor hml_factor

// Display results
use ff3_results, clear
list portfolio alpha market_beta smb_beta hml_beta, clean noobs

// =============================================================================
// Advanced techniques: Heckman selection, treatment effects, 
// PSM, DID, synthetic control, RDD
// =============================================================================

// ******************************
// 8. SELECTION BIAS CORRECTION
// ******************************

// Heckman Selection Model
// Using simulated labor market data for women
webuse womenwk, clear

// Generate selection indicator: 1 if wage observed, 0 otherwise
gen selected = (wage != .)

// First stage: Probit selection equation
probit selected age education married children
predict selection_index, xb
gen inverse_mills = normalden(selection_index) / normal(selection_index)

// Second stage: Wage equation with correction
reg wage education age inverse_mills
estimates store heckman_corrected

// Naive model without selection correction
reg wage education age
estimates store naive

// Compare results
esttab heckman_corrected naive, se star(* 0.1 ** 0.05 *** 0.01)

// ******************************
// 9. TREATMENT EFFECT ESTIMATION
// ******************************

// Treatment Effects Model
// Estimating union membership effect on wages
use http://www.stata-press.com/data/r15/union3, clear

// Maximum likelihood estimation
etregress wage age grade smsa black tenure, ///
    treat(union = south black tenure) 
estimates store treatment_ml

// Two-step estimator
etregress wage age grade smsa black tenure, ///
    treat(union = south black tenure) twostep
estimates store treatment_twostep

// ******************************
// 10. PROPENSITY SCORE MATCHING
// ******************************

// Evaluating job training program effects
use https://users.nber.org/~rdehejia/data/nsw_dw.dta, clear

// Estimate propensity scores
probit treat age educ black hisp married nodegree re74 re75
predict pscore, pr

// Nearest neighbor matching
psmatch2 treat, pscore(pscore) outcome(re78) neighbor(1) common caliper(0.05)

// Balance checks
pstest age educ black hisp married nodegree re74 re75, both

// ******************************
// 11. DIFFERENCE-IN-DIFFERENCES
// ******************************

// Policy implementation analysis
use "https://dss.princeton.edu/training/Panel101.dta", clear

// Basic DID specification
gen treated = (country > 4) & !missing(country)
gen time = (year >= 1994) & !missing(year)
gen did = treated * time

reg y treated time did
estimates store did_basic

// Event study approach
forvalues i = 3(-1)1 {
    gen pre_`i' = (year == (1994 - `i') & treated == 1)
}
gen current = (year == 1994 & treated == 1)
forvalues j = 1/3 {
    gen post_`j' = (year == (1994 + `j') & treated == 1)
}

reg y treated time pre_* current post_* 
estimates store did_event

// ******************************
// 12. SYNTHETIC CONTROL METHOD
// ******************************

// Tobacco control policy evaluation
use http://fmwww.bc.edu/repec/bocode/s/smoking.dta, clear
xtset state year

// Create synthetic California
synth cigsale retprice lnincome age15to24 beer ///
    cigsale(1975) cigsale(1980) cigsale(1988), ///
    trunit(3) trperiod(1989) xperiod(1980(1)1988) ///
    figure keep(synth_result, replace)

// Placebo tests
use synth_result, clear
gen effect = _Y_treated - _Y_synthetic
line effect _time, xline(1989) yline(0) ///
    title("Treatment Effect Over Time") ///
    ytitle("Effect on Cigarette Sales") ///
    xtitle("Year")

// ******************************
// 13. REGRESSION DISCONTINUITY
// ******************************

// Election analysis using vote share data
net get rd 
use votex, clear

// Basic RD specification
rd lne d, mbw(100)
estimates store rd_basic

// Covariate-adjusted RD
rd lne d, mbw(100) cov(i votpop black blucllr farmer fedwrkr ///
    forborn manuf unemplyd union urban veterans)
estimates store rd_adjusted

// Density test for manipulation
DCdensity d, breakpoint(0.5)

// Output results
esttab did_basic did_event, se star(* 0.1 ** 0.05 *** 0.01)
esttab rd_basic rd_adjusted, se star(* 0.1 ** 0.05 *** 0.01)

// ******************************
// 14. GARCH MODELING
// ******************************

// Set working directory and load oil price data
cd "D:\data\chapter6"
import excel using OIL.XLS, clear first case(lower)

// Generate log returns and set time series structure
sort date
gen day = _n
tsset day
gen p = 100 * (ln(spot) - ln(L1.spot)) // Percentage log returns

// Identify ARMA structure for mean equation
arima p, ma(1 3)
predict rres, residuals
corrgram rres, lags(12) // Check residual autocorrelation

// Test for ARCH effects (McLeod-Li test)
gen res2 = rres^2
regress res2 L1.res2 L2.res2 L3.res2 L4.res2
test L1.res2 L2.res2 L3.res2 L4.res2 // Significant F-stat indicates ARCH effects

// Estimate GARCH(1,1) model
arch p, arch(1) garch(1) ma(1 3)

// Diagnostic checks on standardized residuals
predict hres, residuals
predict hv, variance 
gen sdres = hres / sqrt(hv)
gen sdres2 = sdres^2
corrgram sdres, lags(12)  // Check autocorrelation
corrgram sdres2, lags(12) // Check remaining ARCH effects

// ******************************
// 15. TERM SPREAD ANALYSIS
// ******************************

import excel using QUARTERLY.XLS, clear first case(lower)
gen q = _n
tsset q
gen s = r5 - tbill // Term spread calculation

// ARIMA model identification
arima s, ar(1/2) ma(1 7)
predict res, residuals
gen res2 = res^2
corrgram res2, lags(12) // Clear evidence of volatility clustering

// GARCH model estimation and comparison
arch s, arch(1) garch(1) ar(1/2) ma(1 7)
est store garch11
estat ic // Obtain information criteria

arch s, arch(3) ar(1/2) ma(1 7)
est store arch3
estat ic

// Model comparison using IC
estimates stats garch11 arch3

// ******************************
// 16. EQUITY RETURNS ANALYSIS
// ******************************

import excel using NYSEReturns.xlsx, clear first case(lower)
destring return, replace force
sort entry
gen day = _n
tsset day
rename return r

// Volatility modeling with asymmetric effects
arch r, earch(1) egarch(1) ar(1/2) distribution(t)
estat ic

// Leverage effect test
predict s, residuals
gen s2 = s^2
regress s2 L1.s L2.s // Test for asymmetric impact

// ******************************
// 17. VAR ESTIMATION AND BACKTESTING
// ******************************

use http://www.stata-press.com/data/feus/index, clear

// Parametric VaR with normal distribution
su SP500
scalar mean = r(mean)
scalar sd = r(sd)
scalar VaR95 = mean + sd * invnormal(0.05)

// GARCH-VaR with t-distribution
arch SP500, arch(1) garch(1) distribution(t)
predict variance, variance
gen VaR_garch = mean + sqrt(variance) * invt(e(df_r), 0.05)

// Backtesting procedure
gen violations = SP500 < VaR_garch
su violations
scalar T1 = r(sum) // Number of violations
scalar T = r(N)    // Total observations
scalar empirical_coverage = T1/T

// Unconditional coverage test
scalar expected_violations = T * 0.05
scalar LR_uc = 2*(T1*log(T1/expected_violations) + (T-T1)*log((T-T1)/(T-expected_violations)))
di "LR UC statistic: " LR_uc " p-value: " chi2tail(1, LR_uc)

// ******************************
// 18. PORTFOLIO VAR WITH DCC-GARCH
// ******************************

mgarch dcc SP500 FTSE100, arch(1) garch(1)
predict variance*, variance
gen VaR_portfolio = 0.5*sqrt(variance_SP500) + 0.5*sqrt(variance_FTSE100) * invnormal(0.05)

// ******************************
// 19. MONTE CARLO SIMULATION FOR VAR
// ******************************

set obs 10000
set seed 1
drawnorm SP_sim FTSE_sim, n(10000) means(mean_V) sds(sd_V) corr(correlation)
gen port_sim = 0.5*SP_sim + 0.5*FTSE_sim
_pctile port_sim, p(5)
scalar VaR_MC = r(r1)

// =============================================================================
// END OF CODE
// =============================================================================