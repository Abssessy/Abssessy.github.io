// =============================================================================
// This code sample demonstrates proficiency in Stata for econometric analysis,
// covering data management, time-series operations, and panel data techniques.
// All identifying information has been removed using generic variables and
// sample operations to maintain anonymity. The final part includes some case 
// analysis.
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

// ******************************
// Case Analysis
// ******************************

// ******************************
// COVID-19 Early Warning System Social Value Analysis
// ******************************

clear all
set more off

// Define paths
global data_dir "Data"
global output_dir "Output" 
global figures_dir "Output/Figures"

// Create output directories
capture mkdir "$output_dir"
capture mkdir "$figures_dir"

// Load and prepare wave deaths data (only 1-15 weeks)
import delimited "$data_dir/wave_deaths.csv", clear
rename deaths baseline_deaths
rename week relative_week
save "$output_dir/wave_deaths_clean.dta", replace

// Load and prepare vaccination data
import delimited "$data_dir/COVID-19_Vaccinations_in_the_United_States_Jurisdiction.csv", clear

// Keep only national data
keep if location == "US"
keep date additional_doses_vax_pct

// Convert date and calculate weeks since booster authorization
gen date_numeric = date(date, "MDY")
format date_numeric %td
sort date_numeric

gen ref_date = date("22sep2021", "DMY")
gen week_num = floor((date_numeric - ref_date)/7) + 1

// Keep ONLY the last observation per week (latest date in each week)
sort week_num date_numeric
by week_num: keep if _n == _N

// Use the cumulative vaccination percentage
rename additional_doses_vax_pct cum_vax_pct
replace cum_vax_pct = min(cum_vax_pct, 100)

// Keep only first 17 weeks for template (max needed for scenario B)
keep if week_num <= 17

// Create vaccination template
preserve
    keep week_num cum_vax_pct
    rename week_num vaccination_week
    save "$output_dir/vaccination_template.dta", replace
restore

// Define model parameters
local num_waves 4
local wave_starts 40 82 130 205
local vaccine_dev_time 10
local vaccine_efficacy 0.9
local vsl 13500000
local annual_discount_rate 0.04
local weekly_discount_rate = (1 + `annual_discount_rate')^(1/52) - 1

// Initialize results dataset
clear
set obs `num_waves'
gen wave = _n
gen start_week = .
gen deaths_averted = .
gen social_value = .

local i = 1
foreach start in `wave_starts' {
    replace start_week = `start' in `i'
    local i = `i' + 1
}

save "$output_dir/results.dta", replace

// Load vaccination template into memory
use "$output_dir/vaccination_template.dta", clear
tempfile vax_template
save `vax_template'

// Process each variant wave
forvalues w = 1/`num_waves' {
    
    use "$output_dir/results.dta", clear
    local start = start_week[`w']
    
    // Define vaccine timelines
    local vax_available_a = `start' - 4 + `vaccine_dev_time'  // Status quo
    local vax_available_b = `start' - 12 + `vaccine_dev_time' // New technology
    
    display "Processing Wave `w'"
    display "Start week: `start'"
    display "Vaccine available - Scenario A: week `vax_available_a'"
    display "Vaccine available - Scenario B: week `vax_available_b'"
    
    // Create wave timeline (15 weeks as per data)
    clear
    set obs 15
    gen relative_week = _n
    gen absolute_week = `start' + relative_week - 1
    
    // Merge with death pattern
    merge 1:1 relative_week using "$output_dir/wave_deaths_clean.dta", nogen keep(match)
    
    // Initialize vaccination coverage
    gen cumulative_vax_a = 0
    gen cumulative_vax_b = 0
    
    // Calculate vaccination coverage for each week in the wave
    forvalues wk = 1/15 {
        local abs_week = absolute_week[`wk']
        
        // Scenario A: Status quo
        if `abs_week' >= `vax_available_a' {
            local weeks_since_vax = `abs_week' - `vax_available_a' + 1
            
            // Get vaccination rate from template
            preserve
                use `vax_template' if vaccination_week == `weeks_since_vax', clear
                if _N > 0 {
                    local vax_rate_a = cum_vax_pct[1]
                }
                else {
                    // If beyond template, use maximum observed
                    use `vax_template', clear
                    sum cum_vax_pct
                    local vax_rate_a = r(max)
                }
            restore
            
            replace cumulative_vax_a = `vax_rate_a' in `wk'
        }
        
        // Scenario B: New technology  
        if `abs_week' >= `vax_available_b' {
            local weeks_since_vax = `abs_week' - `vax_available_b' + 1
            
            // Get vaccination rate from template
            preserve
                use `vax_template' if vaccination_week == `weeks_since_vax', clear
                if _N > 0 {
                    local vax_rate_b = cum_vax_pct[1]
                }
                else {
                    // If beyond template, use maximum observed
                    use `vax_template', clear
                    sum cum_vax_pct
                    local vax_rate_b = r(max)
                }
            restore
            
            replace cumulative_vax_b = `vax_rate_b' in `wk'
        }
    }
    
    // Calculate deaths under each scenario
    gen deaths_a = baseline_deaths * (1 - `vaccine_efficacy' * (cumulative_vax_a/100))
    gen deaths_b = baseline_deaths * (1 - `vaccine_efficacy' * (cumulative_vax_b/100))
    
    // Calculate discounted deaths averted
    gen deaths_averted_wk = deaths_a - deaths_b
    gen discount_factor = 1/(1 + `weekly_discount_rate')^(absolute_week - 1)
    gen discounted_deaths_averted = deaths_averted_wk * discount_factor
    
    // Summarize for this wave
    sum discounted_deaths_averted
    local wave_deaths_averted = r(sum)
    local wave_social_value = `wave_deaths_averted' * `vsl'
    
    // Store results
    use "$output_dir/results.dta", clear
    replace deaths_averted = `wave_deaths_averted' in `w'
    replace social_value = `wave_social_value' in `w'
    save "$output_dir/results.dta", replace
    
    display "Wave `w': " `wave_deaths_averted' " deaths averted, $" `wave_social_value' " social value"
}

// Calculate and display final results
use "$output_dir/results.dta", clear
egen total_deaths_averted = total(deaths_averted)
egen total_social_value = total(social_value)

gen social_value_millions = social_value / 1000000
gen total_social_value_millions = total_social_value / 1000000
format social_value_millions total_social_value_millions %12.2f

list wave start_week deaths_averted social_value_millions

display "=== FINAL RESULTS ==="
display "Total deaths averted: " total_deaths_averted[1]
display "Total social value: $" total_social_value[1]
display "Total social value (millions): $" total_social_value_millions[1] " million"

// Save detailed results
save "$output_dir/results.dta", replace
export delimited "$output_dir/final_results.csv", replace

// Create charts
use "$output_dir/results.dta", clear
graph bar social_value_millions, over(wave) ///
    title("Social Value by Variant Wave (Millions USD)") ///
    ytitle("Social Value (Millions USD)") ///
    blabel(bar, format(%9.1f))
graph export "$figures_dir/social_value_by_wave.png", replace

use "$output_dir/results.dta", clear
graph bar deaths_averted, over(wave) ///
    title("Deaths Averted by Variant Wave") ///
    ytitle("Number of Deaths Averted") ///
    blabel(bar, format(%9.1f))
graph export "$figures_dir/deaths_averted_by_wave.png", replace

// Create summary table
use "$output_dir/results.dta", clear
local total_deaths = total_deaths_averted[1]
local total_social = total_social_value[1]
local total_social_millions = total_social_value_millions[1]

clear
set obs 1
gen total_deaths_averted = `total_deaths'
gen total_social_value = `total_social'
gen total_social_value_millions = `total_social_millions'

export delimited "$output_dir/summary_results.csv", replace

display "Analysis completed successfully"

* Clear memory and set up environment
clear all
set more off

* ==================================================
* PART 1: SOCIAL CAPITAL
* ==================================================

* Import social capital data
import delimited "q1_social_capital_county.csv", clear

* Create FIPS code with 5-digit format
generate str5 FIPS_code = string(county, "%05.0f")

* Save the social capital dataset
save social_capital, replace

* Import Opportunity Atlas data
import delimited "q1_atlas_outcomes.csv", clear 

* Create FIPS code with 5-digit format
generate str2 state_code = string(state, "%02.0f")
generate str3 county_code = string(county, "%03.0f")
generate str5 FIPS_code = state_code + county_code

* Merge with social capital data
merge 1:1 FIPS_code using social_capital

* Keep only successfully matched observations
keep if _merge == 3

* Clean up variables
drop state county state_code county_code _merge

* Reorder variables for better readability
order FIPS_code county_name ec_county ec_high_county kfr_pooled_pooled_p25 kfr_pooled_pooled_p75

* Remove observations with missing values in key variables
drop if strpos(ec_county, "NA") | strpos(ec_high_county, "NA") | strpos(kfr_pooled_pooled_p25, "NA") | strpos(kfr_pooled_pooled_p75, "NA")

* Convert string variables to numeric
destring ec_county ec_high_county kfr_pooled_pooled_p25 kfr_pooled_pooled_p75, replace

* Save the merged dataset
save data_q1, replace

* Generate summary statistics without National and State average data (where county code is "000")
estpost summarize ec_county ec_high_county kfr_pooled_pooled_p25 kfr_pooled_pooled_p75 if substr(FIPS_code, 3, 3) != "000", detail
esttab . using mysum.tex, ///
    cells("mean(fmt(%9.3f)) min(fmt(%9.3f)) p25(fmt(%9.3f)) p50(fmt(%9.3f)) p75(fmt(%9.3f)) max(fmt(%9.3f)) count(fmt(%9.0f))") ///
    noobs nomtitle nonumber replace ///
    booktabs ///
    varlabels(, blist(ec_county "\multicolumn{8}{l}{\textbf{Panel A: Social Capital}} \\" ///
                      kfr_pooled_pooled_p25 "\multicolumn{8}{l}{\textbf{Panel B: Opportunity Atlas}} \\")) ///
    collabels("Mean" "Min" "p25" "Median" "p75" "Max" "N") ///
    prehead("\begin{table}[htbp]\centering\caption{Summary Statistics}\label{tab:summary}\begin{tabular}{l*{7}{r}}\toprule") ///
    posthead("\midrule") ///
    postfoot("\bottomrule \end{tabular} \\ \footnotesize{\textit{Notes: This table presents summary statistics for key variables.}} \end{table}") ///
    mlabels(none)
			
* Regression analysis for low SES
reg kfr_pooled_pooled_p25 ec_county
local beta_low = round(_b[ec_county], 0.001)
local se_low = round(_se[ec_county], 0.001)
local r2_low = round(e(r2), 0.001)

* Regression analysis for high SES
reg kfr_pooled_pooled_p75 ec_county  
local beta_high = round(_b[ec_county], 0.001)
local se_high = round(_se[ec_county], 0.001)
local r2_high = round(e(r2), 0.001)

* Create a binned scatterplot with two series
binscatter kfr_pooled_pooled_p25 kfr_pooled_pooled_p75 ec_county, ///
    n(20) ///
    line(lfit) ///
    xtitle("Economic Connectedness") ///
    ytitle("Intergenerational Mobility") ///
    title("Economic Connectedness and Intergenerational Mobility by SES") ///
    legend(order(1 "Low SES" 2 "High SES")) ///
    note("Low SES: Slope = `beta_low' (SE = `se_low'), R² = `r2_low'" ///
         "High SES: Slope = `beta_high' (SE = `se_high'), R² = `r2_high'", ///
         size(small))

* Save the binned scatterplot
graph export "binned scatterplot_q1.pdf", replace

* ==================================================
* PART 2: Measurement Error and Attenuation Bias
* ==================================================

* Clear memory and set random seed for reproducibility
clear all
set seed 123
set more off

* Define simulation parameters
local alpha = 0
local beta = 0.002
local income_mean = 50000
local income_sd = 10000
local epsilon_sd = 100

* Create empty dataset to store results
clear
set obs 0
gen n = .
gen sigma_u = .
gen expected_beta_star = .
gen mean_beta_hat = .
gen mean_se = .

* Save empty template
tempfile results
save `results'

* Define combinations of sample size and measurement error levels
foreach n in 250 1000 5000 {
    foreach sigma_u in 0 5000 15000 {
        
        * Create temporary storage for regression coefficients and standard errors
        tempname beta_store se_store
        postfile `beta_store' beta_hat using temp_beta, replace
        postfile `se_store' se_hat using temp_se, replace
        
        * Run 1000 simulations for each parameter combination
        forvalues sim = 1/1000 {
            clear
            set obs `n'
            
            generate Income = rnormal(`income_mean', `income_sd')
            generate u = rnormal(0, `sigma_u')
            generate Income_star = Income + u
            generate epsilon = rnormal(0, `epsilon_sd')
            generate TestScore = `alpha' + `beta' * Income + epsilon
            
            regress TestScore Income_star
            local b_hat = _b[Income_star]
            local se_hat = _se[Income_star]
            
            post `beta_store' (`b_hat')
            post `se_store' (`se_hat')
        }
        
        postclose `beta_store'
        postclose `se_store'
        
        use temp_beta, clear
        sum beta_hat
        local mean_beta_hat = r(mean)
        
        use temp_se, clear
        sum se_hat
        local mean_se_val = r(mean)
        
        local var_income = `income_sd'^2
        local var_u = `sigma_u'^2
        local expected_beta_star_val = `beta' * (`var_income' / (`var_income' + `var_u'))
        
        use `results', clear
        set obs `=_N+1'
        replace n = `n' in `=_N'
        replace sigma_u = `sigma_u' in `=_N'
        replace expected_beta_star = `expected_beta_star_val' in `=_N'
        replace mean_beta_hat = `mean_beta_hat' in `=_N'
        replace mean_se = `mean_se_val' in `=_N'
        
        save `results', replace
    }
}

* Prepare data for esttab
use `results', clear

* Create string variables for display
gen n_str = string(n)
gen sigma_u_str = string(sigma_u)
gen expected_str = string(expected_beta_star, "%9.6f")
gen mean_beta_str = string(mean_beta_hat, "%9.6f")
gen mean_se_str = string(mean_se, "%9.6f")

* Create a single string variable for each row
gen row = n_str + " & " + sigma_u_str + " & " + expected_str + " & " + mean_beta_str + " & " + mean_se_str + " \\"

* Create LaTeX table using file write
file open latex_table using "simulation_results.tex", write replace

file write latex_table "\begin{table}[htbp]" _n
file write latex_table "\centering" _n
file write latex_table "\caption{Simulation Results}" _n
file write latex_table "\label{tab:simulation_results}" _n
file write latex_table "\vspace{0.1cm}" _n  
file write latex_table "\begin{tabular}{ccccc}" _n
file write latex_table "\toprule" _n
file write latex_table "Sample Size & \$\sigma_u\$ & Expected \$\beta^*\$ & Mean \$\hat{\beta}^*\$ & Mean SE \\" _n
file write latex_table "\midrule" _n

* Add midrule between different sample sizes
local current_n = n[1]
forvalues i = 1/`=_N' {
    if `i' > 1 & n[`i'] != `current_n' {
        file write latex_table "\midrule" _n
        local current_n = n[`i']
    }
    file write latex_table (row[`i']) _n
}

file write latex_table "\bottomrule" _n
file write latex_table "\end{tabular}" _n 
file write latex_table "\\" _n 
file write latex_table "\footnotesize{\textit{Notes: This table presents simulation results for different sample sizes and measurement error levels. The expected $\beta^*$ is calculated using the theoretical attenuation formula.}}" _n
file write latex_table "\end{table}" _n

file close latex_table

* Display the table in Stata for verification
display "LaTeX table saved to simulation_results.tex"
list n sigma_u expected_beta_star mean_beta_hat mean_se, clean noobs

* Clean up temporary files
capture erase temp_beta.dta
capture erase temp_se.dta

* ==================================================
* PART 3: Intergenerational Mobility
* ==================================================

* Clear memory and set up environment
clear all
set more off

* Import transition matrices data
import delimited "q3_transition_matrices.csv", clear

* Filter for White and Black combined gender data
keep if (kid_race == "White" | kid_race == "Black") & gender == "P"

* Store row indices for White and Black
gen row_num = _n
sum row_num if kid_race == "White" & gender == "P"
local white_row = r(mean)
sum row_num if kid_race == "Black" & gender == "P"
local black_row = r(mean)

* ==================================================
* PART 3.1: Children's Income Distribution
* ==================================================

* Create a dataset to store children's income distribution results
clear
set obs 10
gen race = ""
gen quintile = .
gen proportion = .

* Fill in race and quintile information
replace race = "White" in 1/5
replace race = "Black" in 6/10

forvalues i = 1/5 {
    replace quintile = `i' in `i'
    replace quintile = `i' in `=`i'+5'
}

* Load data for matrix operations
preserve
import delimited "q3_transition_matrices.csv", clear
keep if (kid_race == "White" | kid_race == "Black") & gender == "P"

* CORRECTED: Construct transition matrices with proper orientation
* For White children - each ROW represents parent quintile, each COLUMN represents child quintile
matrix T_white = [kfr_q1_cond_par_q1[`white_row'], kfr_q2_cond_par_q1[`white_row'], kfr_q3_cond_par_q1[`white_row'], kfr_q4_cond_par_q1[`white_row'], kfr_q5_cond_par_q1[`white_row'] \ ///
                  kfr_q1_cond_par_q2[`white_row'], kfr_q2_cond_par_q2[`white_row'], kfr_q3_cond_par_q2[`white_row'], kfr_q4_cond_par_q2[`white_row'], kfr_q5_cond_par_q2[`white_row'] \ ///
                  kfr_q1_cond_par_q3[`white_row'], kfr_q2_cond_par_q3[`white_row'], kfr_q3_cond_par_q3[`white_row'], kfr_q4_cond_par_q3[`white_row'], kfr_q5_cond_par_q3[`white_row'] \ ///
                  kfr_q1_cond_par_q4[`white_row'], kfr_q2_cond_par_q4[`white_row'], kfr_q3_cond_par_q4[`white_row'], kfr_q4_cond_par_q4[`white_row'], kfr_q5_cond_par_q4[`white_row'] \ ///
                  kfr_q1_cond_par_q5[`white_row'], kfr_q2_cond_par_q5[`white_row'], kfr_q3_cond_par_q5[`white_row'], kfr_q4_cond_par_q5[`white_row'], kfr_q5_cond_par_q5[`white_row']]
                   
matrix P0_white = [par_q1[`white_row'], par_q2[`white_row'], par_q3[`white_row'], par_q4[`white_row'], par_q5[`white_row']]

* For Black children  
matrix T_black = [kfr_q1_cond_par_q1[`black_row'], kfr_q2_cond_par_q1[`black_row'], kfr_q3_cond_par_q1[`black_row'], kfr_q4_cond_par_q1[`black_row'], kfr_q5_cond_par_q1[`black_row'] \ ///
                  kfr_q1_cond_par_q2[`black_row'], kfr_q2_cond_par_q2[`black_row'], kfr_q3_cond_par_q2[`black_row'], kfr_q4_cond_par_q2[`black_row'], kfr_q5_cond_par_q2[`black_row'] \ ///
                  kfr_q1_cond_par_q3[`black_row'], kfr_q2_cond_par_q3[`black_row'], kfr_q3_cond_par_q3[`black_row'], kfr_q4_cond_par_q3[`black_row'], kfr_q5_cond_par_q3[`black_row'] \ ///
                  kfr_q1_cond_par_q4[`black_row'], kfr_q2_cond_par_q4[`black_row'], kfr_q3_cond_par_q4[`black_row'], kfr_q4_cond_par_q4[`black_row'], kfr_q5_cond_par_q4[`black_row'] \ ///
                  kfr_q1_cond_par_q5[`black_row'], kfr_q2_cond_par_q5[`black_row'], kfr_q3_cond_par_q5[`black_row'], kfr_q4_cond_par_q5[`black_row'], kfr_q5_cond_par_q5[`black_row']]
                  
matrix P0_black = [par_q1[`black_row'], par_q2[`black_row'], par_q3[`black_row'], par_q4[`black_row'], par_q5[`black_row']]

* Calculate children's distribution: P_child = P_parent × T
matrix White_child = P0_white * T_white
matrix Black_child = P0_black * T_black

* Display for verification
display "White Children Distribution (Generation 1):"
matrix list White_child
display "Black Children Distribution (Generation 1):"
matrix list Black_child

restore

* Fill in the proportions from matrix calculations
forvalues q = 1/5 {
    replace proportion = White_child[1, `q'] if race == "White" & quintile == `q'
    replace proportion = Black_child[1, `q'] if race == "Black" & quintile == `q'
}

* Format the proportion variable
format proportion %9.4f

* Display the calculated distributions
display "Children's Income Distribution by Race (Generation 1)"
display "=================================================="
list race quintile proportion, sepby(race) noobs

* Create bar chart for children's income distribution
graph bar proportion, over(quintile) over(race) ///
    title("Children's Income Distribution by Race", size(medium)) ///
    subtitle("Generation 1: Based on Intergenerational Transition Matrices", size(small)) ///
    ytitle("Proportion", size(small)) ///
    blabel(bar, format(%9.3f) size(small)) ///
    bar(1, color(navy)) bar(2, color(maroon)) ///
    legend(label(1 "White Children") label(2 "Black Children") size(small)) ///
    ylabel(0(0.1)0.5, angle(0) labsize(small)) ///
    plotregion(margin(medium)) ///
    graphregion(color(white)) ///
    note("Q1 = Bottom quintile (poorest), Q5 = Top quintile (richest)" ///
         "Data source: Chetty et al. (2020)", size(vsmall))

* Export the graph as PDF
graph export "children_income_distribution.pdf", replace

* Verify that proportions sum to 1 for each race
bysort race: egen total_prop = total(proportion)
display "Verification: Proportions sum to 1 for each race"
list race total_prop, noobs

* ==================================================
* PART 3.5: Top Quintile Share Evolution Over Generations
* ==================================================

* Create dataset to store generational evolution results
clear
set obs 22
gen race = ""
gen generation = .
gen top_quintile_share = .

* Fill in race and generation information
forvalues i = 0/10 {
    local pos1 = `i' + 1
    local pos2 = `i' + 12
    replace race = "White" in `pos1'
    replace race = "Black" in `pos2'
    replace generation = `i' in `pos1'
    replace generation = `i' in `pos2'
}

* Use the same matrices from Part 3.1 for consistency
preserve
import delimited "q3_transition_matrices.csv", clear
keep if (kid_race == "White" | kid_race == "Black") & gender == "P"

* White population (same as above)
matrix T_white = [kfr_q1_cond_par_q1[`white_row'], kfr_q2_cond_par_q1[`white_row'], kfr_q3_cond_par_q1[`white_row'], kfr_q4_cond_par_q1[`white_row'], kfr_q5_cond_par_q1[`white_row'] \ ///
                  kfr_q1_cond_par_q2[`white_row'], kfr_q2_cond_par_q2[`white_row'], kfr_q3_cond_par_q2[`white_row'], kfr_q4_cond_par_q2[`white_row'], kfr_q5_cond_par_q2[`white_row'] \ ///
                  kfr_q1_cond_par_q3[`white_row'], kfr_q2_cond_par_q3[`white_row'], kfr_q3_cond_par_q3[`white_row'], kfr_q4_cond_par_q3[`white_row'], kfr_q5_cond_par_q3[`white_row'] \ ///
                  kfr_q1_cond_par_q4[`white_row'], kfr_q2_cond_par_q4[`white_row'], kfr_q3_cond_par_q4[`white_row'], kfr_q4_cond_par_q4[`white_row'], kfr_q5_cond_par_q4[`white_row'] \ ///
                  kfr_q1_cond_par_q5[`white_row'], kfr_q2_cond_par_q5[`white_row'], kfr_q3_cond_par_q5[`white_row'], kfr_q4_cond_par_q5[`white_row'], kfr_q5_cond_par_q5[`white_row']]
                   
matrix P0_white = [par_q1[`white_row'], par_q2[`white_row'], par_q3[`white_row'], par_q4[`white_row'], par_q5[`white_row']]

* Black population (same as above)
matrix T_black = [kfr_q1_cond_par_q1[`black_row'], kfr_q2_cond_par_q1[`black_row'], kfr_q3_cond_par_q1[`black_row'], kfr_q4_cond_par_q1[`black_row'], kfr_q5_cond_par_q1[`black_row'] \ ///
                  kfr_q1_cond_par_q2[`black_row'], kfr_q2_cond_par_q2[`black_row'], kfr_q3_cond_par_q2[`black_row'], kfr_q4_cond_par_q2[`black_row'], kfr_q5_cond_par_q2[`black_row'] \ ///
                  kfr_q1_cond_par_q3[`black_row'], kfr_q2_cond_par_q3[`black_row'], kfr_q3_cond_par_q3[`black_row'], kfr_q4_cond_par_q3[`black_row'], kfr_q5_cond_par_q3[`black_row'] \ ///
                  kfr_q1_cond_par_q4[`black_row'], kfr_q2_cond_par_q4[`black_row'], kfr_q3_cond_par_q4[`black_row'], kfr_q4_cond_par_q4[`black_row'], kfr_q5_cond_par_q4[`black_row'] \ ///
                  kfr_q1_cond_par_q5[`black_row'], kfr_q2_cond_par_q5[`black_row'], kfr_q3_cond_par_q5[`black_row'], kfr_q4_cond_par_q5[`black_row'], kfr_q5_cond_par_q5[`black_row']]
                  
matrix P0_black = [par_q1[`black_row'], par_q2[`black_row'], par_q3[`black_row'], par_q4[`black_row'], par_q5[`black_row']]

restore

* Calculate distribution evolution over generations
* For White population
matrix P_current_white = P0_white
replace top_quintile_share = P_current_white[1,5] if race == "White" & generation == 0

forvalues gen = 1/10 {
    matrix P_new_white = P_current_white * T_white
    replace top_quintile_share = P_new_white[1,5] if race == "White" & generation == `gen'
    matrix P_current_white = P_new_white
}

* For Black population
matrix P_current_black = P0_black
replace top_quintile_share = P_current_black[1,5] if race == "Black" & generation == 0

forvalues gen = 1/10 {
    matrix P_new_black = P_current_black * T_black
    replace top_quintile_share = P_new_black[1,5] if race == "Black" & generation == `gen'
    matrix P_current_black = P_new_black
}

* Format the share variable
format top_quintile_share %9.4f

* Verify consistency between Part 3.1 and Part 3.5
display "Consistency Check: Generation 1 in Part 3.5 should match Part 3.1 results"
list race generation top_quintile_share if generation == 0 | generation == 1, noobs

* Create line graph showing evolution over generations with data labels
twoway (line top_quintile_share generation if race == "White", ///
        lcolor(navy) lwidth(medium) lpattern(solid)) ///
       (scatter top_quintile_share generation if race == "White", ///
        mcolor(navy) msymbol(O) msize(small) ///
        mlabel(top_quintile_share) mlabcolor(navy) mlabposition(12) mlabsize(small) mlabformat(%9.3f)) ///
       (line top_quintile_share generation if race == "Black", ///
        lcolor(maroon) lwidth(medium) lpattern(dash)) ///
       (scatter top_quintile_share generation if race == "Black", ///
        mcolor(maroon) msymbol(D) msize(small) ///
        mlabel(top_quintile_share) mlabcolor(maroon) mlabposition(6) mlabsize(small) mlabformat(%9.3f)), ///
    title("Evolution of Top Income Quintile Share Over Generations", size(medium)) ///
    subtitle("By Race, Assuming Constant Intergenerational Transition Matrix", size(small)) ///
    xtitle("Generation", size(small)) ///
    ytitle("Proportion in Top Quintile", size(small)) ///
    xlabel(0(1)10, labsize(small)) ///
    ylabel(0(0.05)0.5, angle(0) format(%9.2f) labsize(small)) ///
    legend(order(1 "White" 3 "Black") position(6) rows(1) size(small) region(lstyle(none))) ///
    plotregion(margin(small)) ///
    graphregion(color(white) margin(small)) ///
    note("Data source: Chetty et al. (2020)", size(vsmall))

* Export the graph as PDF
graph export "top_quintile_evolution.pdf", replace

* Display numerical results
display "Top Income Quintile Share Evolution Over Generations"
display "=================================================================="
display "Generation    White        Black"
display "------------------------------------------------------------------"
forvalues i = 0/10 {
    local white_share = top_quintile_share[`i' + 1]
    local black_share = top_quintile_share[`i' + 12]
    display "   `i'         " %6.4f `white_share' "      " %6.4f `black_share'
}

* Final completion message
display ""
display "Analysis completed successfully!"
display "Graphs exported as:"
display "- children_income_distribution.pdf (Part 3.1)"
display "- top_quintile_evolution.pdf (Part 3.5)"

// =============================================================================
// END OF CODE
// =============================================================================