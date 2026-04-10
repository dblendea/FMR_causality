*===========================================================================
*  COMPREHENSIVE ANALYSIS — Beat-to-beat Granger Framework for FMR
*  Reproduces ALL variability results from the manuscript:
*    • Table 2  — Individual patient VAR(3) Granger causality (fixed lag 3,
*                 ~100 ms at 30 fps; sensitivity at lags 1 and 2)
*    • Table 3  — Pooled panel VAR at lags 2, 6, 7 (systolic, patient FE)
*    • Table 3  — Full-cycle ALPM analysis (forward/reverse Granger)
*    • Section B2 — Diastolic-only VAR (ALPM sensitivity analysis)
*    • Lag-specific analysis (lag 2 vs lag 6)
*    • Cross-sectional correlations
*    • Figure 4 — Individual-patient Granger-predictive phenotypes
*    • Figure 5 — Two-timescale MR predictability
*
*  Input : MRV_vertical_10mar2026.csv  (same directory)
*  Output: Tables → CSV/Excel; Figures → PDF/PNG
*
*  NOTE: Table 1 (baseline cross-sectional characteristics) is excluded
*       
*  VERSION: 1 Apr 2026 — Fixed lag 3 primary for Table 2
*           Section A: fixed lag 3 primary + sensitivity at lags 1, 2
*===========================================================================

clear all
set more off
capture log close
log using "FMR_Granger_analysis.log", replace text

timer clear
timer on 1

*───────────────────────────────────────────────────────────────
* 0.  LOAD DATA AND LOWERCASE ALL VARIABLE NAMES
*───────────────────────────────────────────────────────────────
import delimited "MRV_vertical_10mar2026.csv", clear varnames(1) encoding("utf-8")
rename *, lower

* Force numeric on key columns (CSV may contain Excel formula strings)
destring mr_area_cm2 la_volume_ml lv_volume_ml alpm_length_mm ///
         mv_annulus_mm phase_systole1_diastole0 patient_id ///
         time_index_patient cycle mr_type_a1v0 frame, ///
         replace force

* Drop any trailing empty rows
drop if missing(patient_id)

sort patient_id time_index_patient

di as text ""
di as result "============================================"
di as result "  DATA LOADED: `c(N)' observations"
di as result "============================================"

*───────────────────────────────────────────────────────────────
* 0a. INFER CARDIAC CYCLE WHERE MISSING
*     Cycle is NaN for many AFMR patients.
*     Infer from phase transitions: 0→1 = new cycle start.
*───────────────────────────────────────────────────────────────
gen cycle_orig = cycle           /* preserve original */

* Within each patient, detect 0→1 phase transition
sort patient_id time_index_patient
by patient_id: gen _phase_lag = phase_systole1_diastole0[_n-1] if _n > 1
by patient_id: gen _new_cycle = (_phase_lag == 0 & phase_systole1_diastole0 == 1) if _n > 1
replace _new_cycle = 0 if _new_cycle == .

* Cumulative sum of transitions = cycle number
by patient_id: gen _cycle_inf = 1 + sum(_new_cycle)

* Use inferred cycle only where original is missing
replace cycle = _cycle_inf if missing(cycle)
drop _phase_lag _new_cycle _cycle_inf

di as text "Cycle inference complete. Remaining missing:"
count if missing(cycle)

*───────────────────────────────────────────────────────────────
* 0b. COMPUTE WITHIN-PATIENT Z-SCORES FROM RAW DATA
*     Systolic-only z-scores (for individual + pooled systolic models)
*     Full-cycle z-scores (for full-cycle ALPM analysis)
*───────────────────────────────────────────────────────────────
capture drop z_*

* --- SYSTOLIC Z-SCORES ---

* MR area: systolic-only z-score
bysort patient_id: egen _mr_sys_mean = mean(mr_area_cm2) if phase_systole1_diastole0 == 1
bysort patient_id: egen _mr_sys_sd   = sd(mr_area_cm2)   if phase_systole1_diastole0 == 1
bysort patient_id: egen mr_sys_mean  = max(_mr_sys_mean)
bysort patient_id: egen mr_sys_sd    = max(_mr_sys_sd)
gen z_mr_sys = (mr_area_cm2 - mr_sys_mean) / mr_sys_sd if phase_systole1_diastole0 == 1 & mr_sys_sd > 0
drop _mr_sys_mean _mr_sys_sd mr_sys_mean mr_sys_sd

* LA volume: all-frames z-score
bysort patient_id: egen _la_mean = mean(la_volume_ml)
bysort patient_id: egen _la_sd   = sd(la_volume_ml)
gen z_lavol_sys = (la_volume_ml - _la_mean) / _la_sd if _la_sd > 0
gen z_lavol_full = z_lavol_sys
drop _la_mean _la_sd

* LV volume: all-frames z-score
bysort patient_id: egen _lv_mean = mean(lv_volume_ml)
bysort patient_id: egen _lv_sd   = sd(lv_volume_ml)
gen z_lvvol_sys = (lv_volume_ml - _lv_mean) / _lv_sd if _lv_sd > 0
gen z_lvvol_full = z_lvvol_sys
drop _lv_mean _lv_sd

* ALPM length: all-frames z-score
bysort patient_id: egen _alpm_mean = mean(alpm_length_mm)
bysort patient_id: egen _alpm_sd   = sd(alpm_length_mm)
gen z_alpm_sys = (alpm_length_mm - _alpm_mean) / _alpm_sd if _alpm_sd > 0
gen z_alpm_full = z_alpm_sys
drop _alpm_mean _alpm_sd

* Annulus: all-frames z-score
bysort patient_id: egen _ann_mean = mean(mv_annulus_mm)
bysort patient_id: egen _ann_sd   = sd(mv_annulus_mm)
gen z_annulus_sys = (mv_annulus_mm - _ann_mean) / _ann_sd if _ann_sd > 0
drop _ann_mean _ann_sd

* --- FULL-CYCLE Z-SCORES (for Section C) ---

* MR area - full cycle
bysort patient_id: egen _mr_fc_mean = mean(mr_area_cm2)
bysort patient_id: egen _mr_fc_sd   = sd(mr_area_cm2)
gen z_mr_fc = (mr_area_cm2 - _mr_fc_mean) / _mr_fc_sd if _mr_fc_sd > 0
drop _mr_fc_mean _mr_fc_sd

* LA volume - full cycle
bysort patient_id: egen _la_fc_mean = mean(la_volume_ml)
bysort patient_id: egen _la_fc_sd   = sd(la_volume_ml)
gen z_lavol_fc = (la_volume_ml - _la_fc_mean) / _la_fc_sd if _la_fc_sd > 0
drop _la_fc_mean _la_fc_sd

* LV volume - full cycle
bysort patient_id: egen _lv_fc_mean = mean(lv_volume_ml)
bysort patient_id: egen _lv_fc_sd   = sd(lv_volume_ml)
gen z_lvvol_fc = (lv_volume_ml - _lv_fc_mean) / _lv_fc_sd if _lv_fc_sd > 0
drop _lv_fc_mean _lv_fc_sd

* ALPM length - full cycle
bysort patient_id: egen _alpm_fc_mean = mean(alpm_length_mm)
bysort patient_id: egen _alpm_fc_sd   = sd(alpm_length_mm)
gen z_alpm_fc = (alpm_length_mm - _alpm_fc_mean) / _alpm_fc_sd if _alpm_fc_sd > 0
drop _alpm_fc_mean _alpm_fc_sd

di as text "Z-scores computed (systolic + full-cycle)."
summarize z_mr_sys z_lavol_sys z_lvvol_sys z_alpm_sys z_annulus_sys
summarize z_mr_fc z_lavol_fc z_lvvol_fc z_alpm_fc


*###########################################################################
*##                                                                       ##
*##    SECTION A:  INDIVIDUAL PATIENT VAR  (TABLE 2)                      ##
*##                                                                       ##
*##    Strategy:                                                           ##
*##      1. For each patient, run VAR at fixed lags 1, 2, and 3           ##
*##      2. Primary analysis: fixed lag 3 (~100 ms at 30 fps)             ##
*##      3. Sensitivity: lags 1 and 2 phenotypes reported alongside       ##
*##      4. IRF signs extracted at lag 3 for all patients                  ##
*##                                                                       ##
*##    Output: Table2_individual_VAR.xlsx (single consolidated table)      ##
*##                                                                       ##
*###########################################################################

di as text ""
di as result "============================================"
di as result "  TABLE 2: INDIVIDUAL PATIENT VAR"
di as result "  (Fixed lag 3 primary; sensitivity 1, 2)"
di as result "============================================"

* Save full dataset before filtering
tempfile _fulldata
save `_fulldata'

keep if phase_systole1_diastole0 == 1
sort patient_id time_index_patient

bysort patient_id: egen _alpm_count = count(z_alpm_sys)

*--- Postfile: one row per patient, all results ---
tempfile t2_results
postfile t2h ///
    patient_id mr_type n_obs model_type ///
    stable_l3 max_eig_l3 ///
    p_LA_l3 sign_LA_l3 fevd_LA_l3 ///
    p_LV_l3 sign_LV_l3 fevd_LV_l3 ///
    p_ALPM_l3 sign_ALPM_l3 fevd_ALPM_l3 ///
    p_ALL_l3 ///
    stable_l1 stable_l2 ///
    p_LA_l1 p_LV_l1 p_ALPM_l1 p_ALL_l1 ///
    p_LA_l2 p_LV_l2 p_ALPM_l2 p_ALL_l2 ///
    sign_LA_l1 sign_LV_l1 sign_ALPM_l1 ///
    sign_LA_l2 sign_LV_l2 sign_ALPM_l2 ///
    fevd_LA_l1 fevd_LV_l1 fevd_ALPM_l1 ///
    fevd_LA_l2 fevd_LV_l2 fevd_ALPM_l2 ///
    using `t2_results', replace

levelsof patient_id, local(patients)

foreach pt of local patients {

    di as text "  Table 2 — Patient `pt'..."

    quietly {

    preserve
    keep if patient_id == `pt'
    sort time_index_patient

    local mr_type  = mr_type_a1v0[1]
    local nobs     = _N
    local has_alpm = (_alpm_count[1] > 0)

    tsset time_index_patient

    * --- Skip if fewer than 20 systolic frames ---
    if `nobs' < 20 {
        post t2h (`pt') (`mr_type') (`nobs') (.) ///
            (.) (.) ///
            (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) ///
            (.) (.) ///
            (.) (.) (.) (.) ///
            (.) (.) (.) (.) ///
            (.) (.) (.) ///
            (.) (.) (.) ///
            (.) (.) (.) (.) (.) (.)
        restore
        continue
    }

    local model_type = cond(`has_alpm' == 1, 1, 2)

    * ──────────────────────────────────────────────────
    *  RUN VAR AT FIXED LAGS 1, 2, 3
    *  Collect p-values, signs, FEVD, stability for each
    * ──────────────────────────────────────────────────

    forvalues LAGNUM = 1/3 {

        * --- Initialize all locals for this lag ---
        local _pLA = .
        local _pLV = .
        local _pALPM = .
        local _pALL = .
        local _sLA = .
        local _sLV = .
        local _sALPM = .
        local _fLA = .
        local _fLV = .
        local _fALPM = .
        local _stab = .
        local _meig = .

        * --- Fit VAR ---
        if `has_alpm' == 1 {
            capture var z_mr_sys z_lavol_sys z_lvvol_sys z_alpm_sys, lags(1/`LAGNUM')
        }
        else {
            capture var z_mr_sys z_lavol_sys z_lvvol_sys, lags(1/`LAGNUM')
        }
        if _rc != 0 {
            * Store results for this lag
            local stable_`LAGNUM' = .
            local p_LA_`LAGNUM' = .
            local p_LV_`LAGNUM' = .
            local p_ALPM_`LAGNUM' = .
            local p_ALL_`LAGNUM' = .
            local sign_LA_`LAGNUM' = .
            local sign_LV_`LAGNUM' = .
            local sign_ALPM_`LAGNUM' = .
            local fevd_LA_`LAGNUM' = .
            local fevd_LV_`LAGNUM' = .
            local fevd_ALPM_`LAGNUM' = .
            local max_eig_`LAGNUM' = .
            continue
        }

        * --- Stability ---
        local _stab = 0
        local _meig = .
        capture varstable
        if _rc == 0 {
            capture matrix _eig = r(Modulus)
            if _rc == 0 {
                local _stab = 1
                local neig = rowsof(_eig)
                local _meig = 0
                forvalues e = 1/`neig' {
                    if _eig[`e',1] > `_meig' local _meig = _eig[`e',1]
                    if _eig[`e',1] >= 1 local _stab = 0
                }
            }
        }

        * --- Granger causality ---
        local la_terms ""
        local lv_terms ""
        local alpm_terms ""
        forvalues l = 1/`LAGNUM' {
            local la_terms "`la_terms' [z_mr_sys]L`l'.z_lavol_sys"
            local lv_terms "`lv_terms' [z_mr_sys]L`l'.z_lvvol_sys"
            if `has_alpm' == 1 {
                local alpm_terms "`alpm_terms' [z_mr_sys]L`l'.z_alpm_sys"
            }
        }

        capture test `la_terms'
        if _rc == 0  local _pLA = r(p)

        capture test `lv_terms'
        if _rc == 0  local _pLV = r(p)

        if `has_alpm' == 1 {
            capture test `alpm_terms'
            if _rc == 0  local _pALPM = r(p)

            capture test `la_terms' `lv_terms' `alpm_terms'
            if _rc == 0  local _pALL = r(p)
        }
        else {
            capture test `la_terms' `lv_terms'
            if _rc == 0  local _pALL = r(p)
        }

        * --- IRF / FEVD ---
        local irf_name = "pt`pt'_l`LAGNUM'"
        capture irf create `irf_name', step(10) set(`irf_name') replace nose
        if _rc == 0 {
            tempfile _main_irf
            save `_main_irf'
            quietly use "`irf_name'.irf", clear

            quietly sum oirf if step==1 & response=="z_mr_sys" & impulse=="z_lavol_sys"
            local _sLA = sign(r(mean))
            quietly sum oirf if step==1 & response=="z_mr_sys" & impulse=="z_lvvol_sys"
            local _sLV = sign(r(mean))

            quietly sum fevd if step==5 & response=="z_mr_sys" & impulse=="z_lavol_sys"
            local _fLA = r(mean) * 100
            quietly sum fevd if step==5 & response=="z_mr_sys" & impulse=="z_lvvol_sys"
            local _fLV = r(mean) * 100

            if `has_alpm' == 1 {
                quietly sum oirf if step==1 & response=="z_mr_sys" & impulse=="z_alpm_sys"
                local _sALPM = sign(r(mean))
                quietly sum fevd if step==5 & response=="z_mr_sys" & impulse=="z_alpm_sys"
                local _fALPM = r(mean) * 100
            }

            quietly use `_main_irf', clear
            capture erase "`irf_name'.irf"
        }

        * --- Store results for this lag ---
        local stable_`LAGNUM' = `_stab'
        local max_eig_`LAGNUM' = `_meig'
        local p_LA_`LAGNUM' = `_pLA'
        local p_LV_`LAGNUM' = `_pLV'
        local p_ALPM_`LAGNUM' = `_pALPM'
        local p_ALL_`LAGNUM' = `_pALL'
        local sign_LA_`LAGNUM' = `_sLA'
        local sign_LV_`LAGNUM' = `_sLV'
        local sign_ALPM_`LAGNUM' = `_sALPM'
        local fevd_LA_`LAGNUM' = `_fLA'
        local fevd_LV_`LAGNUM' = `_fLV'
        local fevd_ALPM_`LAGNUM' = `_fALPM'

    }   /* end LAGNUM loop */

    * ──────────────────────────────────────────────────
    *  POST RESULTS: lag 3 = primary; lags 1,2 = sensitivity
    * ──────────────────────────────────────────────────

    post t2h ///
        (`pt') (`mr_type') (`nobs') (`model_type') ///
        (`stable_3') (`max_eig_3') ///
        (`p_LA_3') (`sign_LA_3') (`fevd_LA_3') ///
        (`p_LV_3') (`sign_LV_3') (`fevd_LV_3') ///
        (`p_ALPM_3') (`sign_ALPM_3') (`fevd_ALPM_3') ///
        (`p_ALL_3') ///
        (`stable_1') (`stable_2') ///
        (`p_LA_1') (`p_LV_1') (`p_ALPM_1') (`p_ALL_1') ///
        (`p_LA_2') (`p_LV_2') (`p_ALPM_2') (`p_ALL_2') ///
        (`sign_LA_1') (`sign_LV_1') (`sign_ALPM_1') ///
        (`sign_LA_2') (`sign_LV_2') (`sign_ALPM_2') ///
        (`fevd_LA_1') (`fevd_LV_1') (`fevd_ALPM_1') ///
        (`fevd_LA_2') (`fevd_LV_2') (`fevd_ALPM_2')

    restore
    }   /* end quietly */
}   /* end patient loop */

postclose t2h


*───────────────────────────────────────────────────────────────
*  LABEL, CLASSIFY PHENOTYPES, PRINT, AND EXPORT
*───────────────────────────────────────────────────────────────

preserve
use `t2_results', clear

* --- Labels ---
label var patient_id       "Patient ID"
label var mr_type          "MR type (1=AFMR, 0=VFMR)"
label var n_obs            "N systolic frames"
label var model_type       "Model (1=Full, 2=Reduced)"
label var stable_l3        "VAR stable at lag 3"
label var max_eig_l3       "Max eigenvalue at lag 3"
label var p_LA_l3          "Granger p: LA vol -> MR (lag 3)"
label var sign_LA_l3       "IRF sign: LA -> MR (lag 3)"
label var fevd_LA_l3       "FEVD% h=5: LA vol (lag 3)"
label var p_LV_l3          "Granger p: LV vol -> MR (lag 3)"
label var sign_LV_l3       "IRF sign: LV -> MR (lag 3)"
label var fevd_LV_l3       "FEVD% h=5: LV vol (lag 3)"
label var p_ALPM_l3        "Granger p: ALPM -> MR (lag 3)"
label var sign_ALPM_l3     "IRF sign: ALPM -> MR (lag 3)"
label var fevd_ALPM_l3     "FEVD% h=5: ALPM (lag 3)"
label var p_ALL_l3         "Granger p: ALL -> MR (lag 3)"

* --- Phenotype classification: PRIMARY = lag 3 ---
gen str20 phenotype = ""
replace phenotype = "LV+PM"        if model_type==1 & p_LV_l3 < 0.05 & p_ALPM_l3 < 0.05
replace phenotype = "LA+LV"        if model_type==1 & p_LA_l3 < 0.05 & p_LV_l3 < 0.05 & (p_ALPM_l3 >= 0.05 | p_ALPM_l3==.)
replace phenotype = "LV-dominant"  if model_type==1 & p_LV_l3 < 0.05 & (p_ALPM_l3 >= 0.05 | p_ALPM_l3==.) & p_LA_l3 >= 0.05
replace phenotype = "PM-dominant"  if model_type==1 & p_ALPM_l3 < 0.05 & (p_LV_l3 >= 0.05 | p_LV_l3==.) & p_LA_l3 >= 0.05
replace phenotype = "LA-dominant"  if model_type==1 & p_LA_l3 < 0.05 & p_LV_l3 >= 0.05 & (p_ALPM_l3 >= 0.05 | p_ALPM_l3==.)
replace phenotype = "Neither"      if model_type==1 & p_LV_l3 >= 0.05 & (p_ALPM_l3 >= 0.05 | p_ALPM_l3==.) & p_LA_l3 >= 0.05
replace phenotype = "LV-dominant"  if model_type==2 & p_LV_l3 < 0.05 & p_LA_l3 >= 0.05
replace phenotype = "LA+LV"        if model_type==2 & p_LV_l3 < 0.05 & p_LA_l3 < 0.05
replace phenotype = "LA-dominant"  if model_type==2 & p_LA_l3 < 0.05 & p_LV_l3 >= 0.05
replace phenotype = "Neither"      if model_type==2 & p_LV_l3 >= 0.05 & p_LA_l3 >= 0.05
replace phenotype = "Excluded"     if n_obs < 20 | missing(stable_l3)
* Mark unstable lag-3 models
replace phenotype = phenotype + " {&dagger;}" if stable_l3 == 0 & phenotype != "Excluded" & phenotype != ""

* Sensitivity phenotypes at lags 1 and 2
foreach L in 1 2 {
    gen str20 phenotype_l`L' = ""
    replace phenotype_l`L' = "LV+PM"        if model_type==1 & p_LV_l`L' < 0.05 & p_ALPM_l`L' < 0.05
    replace phenotype_l`L' = "LA+LV"        if model_type==1 & p_LA_l`L' < 0.05 & p_LV_l`L' < 0.05 & (p_ALPM_l`L' >= 0.05 | p_ALPM_l`L'==.)
    replace phenotype_l`L' = "LV-dominant"  if model_type==1 & p_LV_l`L' < 0.05 & (p_ALPM_l`L' >= 0.05 | p_ALPM_l`L'==.) & p_LA_l`L' >= 0.05
    replace phenotype_l`L' = "PM-dominant"  if model_type==1 & p_ALPM_l`L' < 0.05 & (p_LV_l`L' >= 0.05 | p_LV_l`L'==.) & p_LA_l`L' >= 0.05
    replace phenotype_l`L' = "LA-dominant"  if model_type==1 & p_LA_l`L' < 0.05 & p_LV_l`L' >= 0.05 & (p_ALPM_l`L' >= 0.05 | p_ALPM_l`L'==.)
    replace phenotype_l`L' = "Neither"      if model_type==1 & p_LV_l`L' >= 0.05 & (p_ALPM_l`L' >= 0.05 | p_ALPM_l`L'==.) & p_LA_l`L' >= 0.05
    replace phenotype_l`L' = "LV-dominant"  if model_type==2 & p_LV_l`L' < 0.05 & p_LA_l`L' >= 0.05
    replace phenotype_l`L' = "LA+LV"        if model_type==2 & p_LV_l`L' < 0.05 & p_LA_l`L' < 0.05
    replace phenotype_l`L' = "LA-dominant"  if model_type==2 & p_LA_l`L' < 0.05 & p_LV_l`L' >= 0.05
    replace phenotype_l`L' = "Neither"      if model_type==2 & p_LV_l`L' >= 0.05 & p_LA_l`L' >= 0.05
    replace phenotype_l`L' = "Excluded"     if n_obs < 20 | missing(stable_l3)
    * Mark unstable models
    replace phenotype_l`L' = phenotype_l`L' + " {&dagger;}" if stable_l`L' == 0 & phenotype_l`L' != "Excluded" & phenotype_l`L' != ""
}

gen str4 irf_la_str   = cond(sign_LA_l3 == 1, "{&uarr;}", cond(sign_LA_l3 == -1, "{&darr;}", ""))
gen str4 irf_lv_str   = cond(sign_LV_l3 == 1, "{&uarr;}", cond(sign_LV_l3 == -1, "{&darr;}", ""))
gen str4 irf_alpm_str = cond(sign_ALPM_l3 == 1, "{&uarr;}", cond(sign_ALPM_l3 == -1, "{&darr;}", ""))
gen str5 group_str    = cond(mr_type == 1, "AFMR", "VFMR")

* --- Print Table 2 ---
di as text ""
di as result "═══════════════════════════════════════════════════════════"
di as result "  TABLE 2: Individual Patient VAR(3) Granger Causality"
di as result "  Primary: fixed lag 3 (~100 ms at 30 fps)"
di as result "═══════════════════════════════════════════════════════════"
list patient_id group_str n_obs stable_l3 phenotype ///
     p_LA_l3 irf_la_str fevd_LA_l3 ///
     p_LV_l3 irf_lv_str fevd_LV_l3 ///
     p_ALPM_l3 irf_alpm_str fevd_ALPM_l3 ///
     phenotype_l1 phenotype_l2, ///
     noobs abbreviate(14) separator(0)

di as text ""
di as result "--- Phenotype Distribution (lag 3 primary) ---"
tab phenotype group_str, col

di as text ""
di as result "--- Stability Summary (lag 3) ---"
tab stable_l3 group_str if phenotype != "Excluded"

* --- Export ---
export delimited using "Table2_individual_VAR.csv", replace
export excel using "Table2_individual_VAR.xlsx", firstrow(varlabels) replace
di as result "  Table 2 exported."

* ─────────────────────────────────────────────
*  FIGURE 5: Individual-patient phenotype plot
* ─────────────────────────────────────────────
encode phenotype, gen(pheno_num)
gen pt_order = _n

graph hbar (count), over(phenotype, sort(1) descending label(labsize(small))) ///
    over(group_str, label(labsize(small))) ///
    bar(1, color(navy)) ///
    ytitle("Number of patients") ///
    title("Individual-Patient Granger-Predictive Phenotypes", size(medium)) ///
    subtitle("Classification by dominant predictor of MR area (p<0.05)", size(small)) ///
    note("Based on Wald {&chi}{sup:2} tests from individual VAR(3) models (~100 ms at 30 fps)." ///
         "FEVD at horizon 5. Excluded = <20 systolic frames.", size(vsmall)) ///
    legend(off) ///
    scheme(s2color)
graph export "Figure5_phenotypes.pdf", replace as(pdf)
graph export "Figure5_phenotypes.png", replace as(png) width(1920)
di as result "  Figure 5 exported."

tempfile t2_data
save `t2_data'

restore

* Reload full dataset for subsequent sections
use `_fulldata', clear



*###########################################################################
*###########################################################################
*##                                                                       ##
*##    SECTION B:  POOLED PANEL VAR  (TABLE 3)                            ##
*##                                                                       ##
*##    - 5 endogenous: MR, LA vol, LV vol, ALPM sys, annulus              ##
*##    - Exogenous: patient dummies (fixed effects)                       ##
*##    - Lags 2, 6, 7 (Table 3 multi-lag analysis)                        ##
*##    - Three subgroups: AFMR, VFMR, Overall                             ##
*##                                                                       ##
*###########################################################################
*###########################################################################

di as text ""
di as result "============================================"
di as result "  TABLE 3: POOLED PANEL VAR (LAGS 2, 6, 7)"
di as result "============================================"

preserve

* --- Step 1: Keep systole only ---
keep if phase_systole1_diastole0 == 1
sort patient_id time_index_patient

* --- Step 2: Drop patients with no systolic ALPM ---
bysort patient_id: egen _alpm_n = count(z_alpm_sys)
drop if _alpm_n == 0

di as text "  Obs after systolic ALPM filter: " _N

gen group = mr_type_a1v0

local varlist z_mr_sys z_lavol_sys z_lvvol_sys z_alpm_sys z_annulus_sys

tempfile pooled_data
save `pooled_data'

* --- Loop over AFMR, VFMR, Overall × Lags 2, 6, 7 ---
foreach grp in afmr vfmr overall {
    foreach lag in 2 6 7 {

    use `pooled_data', clear

    if "`grp'" == "afmr"  keep if group == 1
    if "`grp'" == "vfmr"  keep if group == 0

    * Patient dummies (fixed effects)
    capture drop pt_dum*
    quietly tabulate patient_id, gen(pt_dum)
    drop pt_dum1
    unab dumlist : pt_dum*

    * Time-set using the global systolic counter
    capture tsset time_index_combined_s
    if _rc != 0 {
        gen _t = _n
        tsset _t
    }

    local n_obs = _N
    local n_pts = 0
    quietly tab patient_id
    local n_pts = r(r)

    di as text ""
    di as result "════════════════════════════════════════"
    di as result "  `grp' — lag `lag' — N obs = `n_obs', N patients = `n_pts'"
    di as result "════════════════════════════════════════"

    * --- Fit VAR ---
    capture noisily var `varlist', lags(1/`lag') exog(`dumlist')
    if _rc != 0 {
        di as error "  VAR failed for `grp' lag `lag'"
        continue
    }

    * --- Stability ---
    di as text "--- Stability ---"
    varstable

    * --- Granger causality ---
    di as text "--- Granger causality -> MR area ---"
    foreach pred in z_lavol_sys z_lvvol_sys z_alpm_sys z_annulus_sys {
        local terms ""
        forvalues l = 1/`lag' {
            local terms "`terms' [z_mr_sys]L`l'.`pred'"
        }
        capture test `terms'
        if _rc == 0 {
            di as text "    `pred' -> z_mr_sys: chi2 = " %8.3f r(chi2) "  p = " %6.4f r(p)
        }
    }

    * --- Joint Granger test (ALL predictors -> MR) ---
    local all_terms ""
    foreach pred in z_lavol_sys z_lvvol_sys z_alpm_sys z_annulus_sys {
        forvalues l = 1/`lag' {
            local all_terms "`all_terms' [z_mr_sys]L`l'.`pred'"
        }
    }
    capture test `all_terms'
    if _rc == 0 {
        di as text "    ALL -> z_mr_sys (joint): chi2 = " %8.3f r(chi2) "  p = " %6.4f r(p)
    }

    * --- IRF/FEVD ---
    local irf_name = "`grp'_lag`lag'"
    capture irf create `irf_name', step(10) set(`irf_name') replace nose
    if _rc == 0 {
        di as text "--- FEVD at horizon 5 ---"
        irf table fevd, ///
            impulse(z_lavol_sys z_lvvol_sys z_alpm_sys z_annulus_sys z_mr_sys) ///
            response(z_mr_sys) set(`irf_name') noci

        di as text "--- OIRF at step 1 ---"
        irf table oirf, ///
            impulse(z_lavol_sys z_lvvol_sys z_alpm_sys z_annulus_sys) ///
            response(z_mr_sys) set(`irf_name') noci
    }

    }  /* end lag loop */
}  /* end grp loop */

restore     /* back to full dataset */


*###########################################################################
*###########################################################################
*##                                                                       ##
*##    SECTION B2: DIASTOLIC-ONLY VAR (ALPM sensitivity analysis)         ##
*##                                                                       ##
*##    - 4 endogenous: MR, LA vol, LV vol, ALPM (z_fc full-cycle vars)   ##
*##    - No annulus (z-scores computed from systolic frames only)          ##
*##    - Exogenous: patient dummies (fixed effects)                       ##
*##    - Lag 1, diastolic frames only                                     ##
*##    - Tests whether ALPM predicts diastolic MR                         ##
*##    - Three subgroups: AFMR, VFMR, Overall                             ##
*##                                                                       ##
*###########################################################################
*###########################################################################

di as text ""
di as result "============================================"
di as result "  DIASTOLIC-ONLY VAR (ALPM sensitivity)"
di as result "============================================"

preserve

* --- Keep diastole only ---
keep if phase_systole1_diastole0 == 0
sort patient_id time_index_patient

* --- Drop patients with no diastolic ALPM ---
bysort patient_id: egen _alpm_dia_n = count(z_alpm_fc)
drop if _alpm_dia_n == 0

di as text "  Diastolic obs with ALPM: " _N

gen group = mr_type_a1v0

* 4-variable model: no annulus (annulus z-scores from systolic frames only)
local dia_varlist z_mr_fc z_lavol_fc z_lvvol_fc z_alpm_fc

tempfile dia_data
save `dia_data'

foreach grp in afmr vfmr overall {

    use `dia_data', clear

    if "`grp'" == "afmr"  keep if group == 1
    if "`grp'" == "vfmr"  keep if group == 0

    * Patient dummies (fixed effects)
    capture drop pt_dum*
    quietly tabulate patient_id, gen(pt_dum)
    drop pt_dum1
    unab dumlist : pt_dum*

    * Time-set (use _n; time_index_combined_s has NaN for some diastolic frames)
    sort patient_id time_index_patient
    capture drop _t
    gen _t = _n
    tsset _t

    local n_obs = _N
    local n_pts = 0
    quietly tab patient_id
    local n_pts = r(r)

    di as text ""
    di as result "════════════════════════════════════════"
    di as result "  `grp' DIASTOLIC — lag 1 — N obs = `n_obs', N patients = `n_pts'"
    di as result "════════════════════════════════════════"

    * Count diastolic MR > 0
    count if z_mr_fc != 0
    di as text "  Frames with non-zero MR: " r(N) " of `n_obs'"

    * --- Fit VAR lag 1 ---
    capture noisily var `dia_varlist', lags(1) exog(`dumlist')
    if _rc != 0 {
        di as error "  VAR failed for `grp' diastolic"
        continue
    }

    * --- Stability ---
    di as text "--- Stability ---"
    varstable

    * --- Granger causality ---
    di as text "--- Granger causality -> diastolic MR ---"
    foreach pred in z_lavol_fc z_lvvol_fc z_alpm_fc {
        capture test [z_mr_fc]L1.`pred'
        if _rc == 0 {
            di as text "    `pred' -> z_mr_fc: chi2 = " %8.3f r(chi2) "  p = " %6.4f r(p)
        }
    }

    * --- Granger causality: reverse (each -> ALPM) ---
    di as text "--- Granger causality -> diastolic ALPM (reverse) ---"
    foreach pred in z_mr_fc z_lavol_fc z_lvvol_fc {
        capture test [z_alpm_fc]L1.`pred'
        if _rc == 0 {
            di as text "    `pred' -> z_alpm_fc: chi2 = " %8.3f r(chi2) "  p = " %6.4f r(p)
        }
    }

    * --- IRF/FEVD ---
    local irf_name = "`grp'_dia"
    capture irf create `irf_name', step(10) set(`irf_name') replace nose
    if _rc == 0 {
        di as text "--- FEVD at horizon 5 (MR equation) ---"
        irf table fevd, ///
            impulse(z_mr_fc z_lavol_fc z_lvvol_fc z_alpm_fc) ///
            response(z_mr_fc) set(`irf_name') noci

        di as text "--- FEVD at horizon 5 (ALPM equation) ---"
        irf table fevd, ///
            impulse(z_mr_fc z_lavol_fc z_lvvol_fc z_alpm_fc) ///
            response(z_alpm_fc) set(`irf_name') noci

        di as text "--- OIRF at step 1 (MR equation) ---"
        irf table oirf, ///
            impulse(z_lavol_fc z_lvvol_fc z_alpm_fc) ///
            response(z_mr_fc) set(`irf_name') noci

        di as text "--- OIRF at step 1 (ALPM equation) ---"
        irf table oirf, ///
            impulse(z_mr_fc z_lavol_fc z_lvvol_fc) ///
            response(z_alpm_fc) set(`irf_name') noci
    }
}

restore     /* back to full dataset */
*##    SECTION F:  LAG-SPECIFIC ANALYSIS (LAG 2 vs LAG 6)                 ##
*##                                                                       ##
*##    4-variable systolic-only model with patient fixed effects           ##
*##    No diastolic ALPM exogenous (pure endogenous comparison)            ##
*##                                                                       ##
*###########################################################################
*###########################################################################

di as text ""
di as result "============================================"
di as result "  LAG-SPECIFIC ANALYSIS"
di as result "============================================"

preserve
keep if phase_systole1_diastole0 == 1
sort patient_id time_index_patient

* Drop patients without ALPM
bysort patient_id: egen _alpm_n2 = count(z_alpm_sys)
drop if _alpm_n2 == 0

gen group = mr_type_a1v0

local lag_varlist z_mr_sys z_lavol_sys z_lvvol_sys z_alpm_sys

tempfile lag_data
save `lag_data'

* Create results table for Figure 6
tempfile lag_results
postfile lag_handle ///
    str6 grp_name lag n_obs ///
    p_LA p_LV p_ALPM ///
    using `lag_results', replace

foreach grp in afmr vfmr {
    foreach lag in 2 6 {

        use `lag_data', clear

        if "`grp'" == "afmr" keep if group == 1
        if "`grp'" == "vfmr" keep if group == 0

        capture drop pt_dum*
        quietly tabulate patient_id, gen(pt_dum)
        drop pt_dum1
        unab dumlist : pt_dum*

        capture tsset time_index_combined_s
        if _rc != 0 {
            gen _t = _n
            tsset _t
        }

        local n_obs = _N

        di as text ""
        di as result "  `grp' — lag `lag' — N = `n_obs'"

        capture noisily var `lag_varlist', lags(1/`lag') exog(`dumlist')
        if _rc != 0 {
            di as error "  VAR failed"
            post lag_handle ("`grp'") (`lag') (`n_obs') (.) (.) (.)
            continue
        }

        * Extract p-values via explicit test commands
        local p_LA = .
        local p_LV = .
        local p_ALPM = .
        
        * Build lag coefficient list dynamically
        local la_terms ""
        local lv_terms ""
        local alpm_terms ""
        forvalues l = 1/`lag' {
            local la_terms "`la_terms' [z_mr_sys]L`l'.z_lavol_sys"
            local lv_terms "`lv_terms' [z_mr_sys]L`l'.z_lvvol_sys"
            local alpm_terms "`alpm_terms' [z_mr_sys]L`l'.z_alpm_sys"
        }
        
        capture test `la_terms'
        if _rc == 0  local p_LA = r(p)
        
        capture test `lv_terms'
        if _rc == 0  local p_LV = r(p)
        
        capture test `alpm_terms'
        if _rc == 0  local p_ALPM = r(p)

        di as text "    LA -> MR: p = " %7.3f `p_LA'
        di as text "    LV -> MR: p = " %7.3f `p_LV'
        di as text "  ALPM -> MR: p = " %7.3f `p_ALPM'

        post lag_handle ("`grp'") (`lag') (`n_obs') (`p_LA') (`p_LV') (`p_ALPM')
    }
}

postclose lag_handle

* --- Load lag results for Figure 6 ---
use `lag_results', clear

di as text ""
di as result "--- Lag-Specific Summary ---"
list, noobs separator(2)

export delimited using "Lag_specific_results.csv", replace

* ─────────────────────────────────────────────
*  FIGURE 6: Two-timescale MR predictability
*  Grouped bar chart: predictors side-by-side
* ─────────────────────────────────────────────
* Reshape to long form: one row per grp × lag × predictor
gen id = _n
reshape long p_, i(id) j(predictor) string
rename p_ pval

* Create plotting variables
gen nlogp = -log10(pval)
gen grp_num = cond(grp_name == "afmr", 0, 1)
gen pred_num = cond(predictor == "LA", 1, cond(predictor == "LV", 2, 3))
gen str12 pred_label = cond(predictor == "LA", "LA vol", ///
    cond(predictor == "LV", "LV vol", "ALPM"))

* X position: group (0/1) + predictor offset (-0.25, 0, +0.25)
gen xpos = grp_num + (pred_num - 2) * 0.25

local sig_line = -log10(0.05)

* Panel A: Within-beat (lag 2)
twoway ///
    (bar nlogp xpos if lag == 2 & predictor == "LV", ///
        barwidth(0.2) color(navy) ) ///
    (bar nlogp xpos if lag == 2 & predictor == "LA", ///
        barwidth(0.2) color(cranberry) ) ///
    (bar nlogp xpos if lag == 2 & predictor == "ALPM", ///
        barwidth(0.2) color(forest_green) ) ///
    , ///
    xlabel(0 "AFMR" 1 "VFMR", labsize(small)) ///
    ylabel(0(0.5)3, labsize(small)) ///
    ytitle("-log{subscript:10}(p)", size(small)) ///
    title("A. Within-beat (lag 2, ~67 ms)", size(medsmall)) ///
    yline(`sig_line', lcolor(red) lpattern(dash) lwidth(thin)) ///
    legend(order(1 "LV vol" 2 "LA vol" 3 "ALPM") rows(1) size(vsmall)) ///
    scheme(s2color) ///
    name(lag2_panel, replace)

* Panel B: Cross-beat (lag 6)
twoway ///
    (bar nlogp xpos if lag == 6 & predictor == "LV", ///
        barwidth(0.2) color(navy) ) ///
    (bar nlogp xpos if lag == 6 & predictor == "LA", ///
        barwidth(0.2) color(cranberry) ) ///
    (bar nlogp xpos if lag == 6 & predictor == "ALPM", ///
        barwidth(0.2) color(forest_green) ) ///
    , ///
    xlabel(0 "AFMR" 1 "VFMR", labsize(small)) ///
    ylabel(0(0.5)3, labsize(small)) ///
    ytitle("-log{subscript:10}(p)", size(small)) ///
    title("B. Cross-beat (lag 6, ~200 ms)", size(medsmall)) ///
    yline(`sig_line', lcolor(red) lpattern(dash) lwidth(thin)) ///
    legend(order(1 "LV vol" 2 "LA vol" 3 "ALPM") rows(1) size(vsmall)) ///
    scheme(s2color) ///
    name(lag6_panel, replace)

* Panel C: Summary significance table
* (Printed to log; visual summary in combined figure)
di as text ""
di as result "--- Figure 6C: Significance summary ---"
di as text "  Lag 2 (within-beat):"
list grp_name predictor pval if lag == 2, noobs separator(3)
di as text "  Lag 6 (cross-beat):"
list grp_name predictor pval if lag == 6, noobs separator(3)

graph combine lag2_panel lag6_panel, ///
    rows(1) imargin(small) ///
    title("Two-Timescale MR Predictability", size(medium)) ///
    subtitle("Granger causality significance by lag and FMR subtype", size(small)) ///
    note("Red dashed line = p=0.05. Patient fixed effects included." ///
         "AFMR n=17, VFMR n=20 (patients with ALPM visibility).", size(vsmall)) ///
    scheme(s2color)
graph export "Figure6_lag_specific.pdf", replace as(pdf)
graph export "Figure6_lag_specific.png", replace as(png) width(1920)
di as result "  Figure 6 exported."

restore     /* back to full dataset */


*###########################################################################
*###########################################################################
*##                                                                       ##
*##    SECTION D:  CROSS-SECTIONAL CORRELATIONS                           ##
*##                                                                       ##
*###########################################################################
*###########################################################################

di as text ""
di as result "============================================"
di as result "  CROSS-SECTIONAL CORRELATIONS"
di as result "============================================"

preserve
keep if phase_systole1_diastole0 == 1

* Patient-level means
collapse (mean) mr_area_cm2 la_volume_ml lv_volume_ml ///
    alpm_length_mm mv_annulus_mm, by(patient_id mr_type_a1v0)

di as text ""
di as result "--- Correlation: LA volume vs MR area (all patients) ---"
pwcorr la_volume_ml mr_area_cm2, sig star(0.05)

di as text ""
di as result "--- Correlation: LV volume vs MR area ---"
pwcorr lv_volume_ml mr_area_cm2, sig star(0.05)

di as text ""
di as result "--- Correlation: ALPM length vs MR area ---"
pwcorr alpm_length_mm mr_area_cm2, sig star(0.05)

di as text ""
di as result "--- Correlation: Annulus vs MR area ---"
pwcorr mv_annulus_mm mr_area_cm2, sig star(0.05)

di as text ""
di as result "--- All correlations (matrix) ---"
pwcorr mr_area_cm2 la_volume_ml lv_volume_ml alpm_length_mm mv_annulus_mm, ///
    sig star(0.05) obs

restore


*###########################################################################
*###########################################################################
*##                                                                       ##
*##    SECTION E:  SUMMARY TABLE (COMBINED)                               ##
*##                                                                       ##
*###########################################################################
*###########################################################################

di as text ""
di as result "============================================"
di as result "  ANALYSIS SUMMARY"
di as result "============================================"

* Reload Table 2 for summary statistics
preserve
use `t2_data', clear

* Analysed patients
count if phenotype != "Excluded"
local n_analysed = r(N)
count
local n_total = r(N)

di as text "  Analysed: `n_analysed' / `n_total' patients"

* Stable models
count if stable_l3 == 1 & phenotype != "Excluded"
di as text "  Stable VAR(3) models: " r(N) " / `n_analysed'"

* Unstable models
list patient_id group_str max_eig_l3 if stable_l3 == 0 & phenotype != "Excluded", noobs

* Median systolic frames
summarize n_obs if phenotype != "Excluded", detail
di as text "  Median systolic frames: " r(p50) " (range " r(min) " - " r(max) ")"

* Phenotype distribution by subtype
di as text ""
di as result "--- AFMR phenotype distribution ---"
tab phenotype if mr_type == 1 & phenotype != "Excluded"
di as text ""
di as result "--- VFMR phenotype distribution ---"
tab phenotype if mr_type == 0 & phenotype != "Excluded"

* Key individual results
di as text ""
di as result "--- Key individual results ---"
foreach pt in 3 23 24 25 26 28 29 32 34 36 {
    list patient_id phenotype p_LA_l3 fevd_LA_l3 p_LV_l3 fevd_LV_l3 p_ALPM_l3 fevd_ALPM_l3 ///
        if patient_id == `pt', noobs
}

restore



*###########################################################################
*###########################################################################
*##                                                                       ##
*##    SECTION C:  FULL-CYCLE ALPM VAR/GRANGER ANALYSIS                   ##
*##                                                                       ##
*##    Purpose: Test whether ALPM predicts MR (and vice versa) when       ##
*##    diastolic frames are included, using z-scores standardised         ##
*##    across the full cardiac cycle.                                     ##
*##                                                                       ##
*##    Models:                                                            ##
*##      Multivariate: 4 endogenous (MR, LA vol, LV vol, ALPM)           ##
*##      + patient fixed effects                                          ##
*##      Lags: 1, 2, 4, 7 (Table 3 extension rows)                       ##
*##      Subgroups: AFMR, VFMR                                           ##
*##                                                                       ##
*##    Tests:                                                             ##
*##      A. ALPM -> MR (forward Granger)                                  ##
*##      B. MR -> ALPM (reverse Granger)                                  ##
*##      C. LV vol -> ALPM, LA vol -> ALPM (upstream pathways)            ##
*##      D. FEVD at horizon 10 for MR and ALPM equations                  ##
*##                                                                       ##
*##    Note: No annulus in this model (annulus z-scores were computed      ##
*##    from systolic frames only and are not available as full-cycle       ##
*##    z-scores for diastolic frames).                                    ##
*##                                                                       ##
*###########################################################################
*###########################################################################

preserve

bysort patient_id: egen _alpm_fc_n = count(z_alpm_fc)
drop if _alpm_fc_n == 0

di as text ""
di as text "  Full-cycle obs with ALPM: " _N

* Count by subtype
count if mr_type_a1v0 == 1
local n_afmr_obs = r(N)
quietly tab patient_id if mr_type_a1v0 == 1
local n_afmr_pts = r(r)

count if mr_type_a1v0 == 0
local n_vfmr_obs = r(N)
quietly tab patient_id if mr_type_a1v0 == 0
local n_vfmr_pts = r(r)

di as text "  AFMR: `n_afmr_pts' patients, `n_afmr_obs' observations"
di as text "  VFMR: `n_vfmr_pts' patients, `n_vfmr_obs' observations"

gen group_fc = mr_type_a1v0

sort patient_id time_index_patient

local fc_varlist z_mr_fc z_lavol_fc z_lvvol_fc z_alpm_fc

tempfile fc_data
save `fc_data'

*##    SECTION C1: FORWARD GRANGER (ALPM -> MR) + FEVD                     ##
*##               Multivariate model, lags 1, 2, 4, 7                     ##
*##               For Table 3 extension rows                              ##
*##                                                                       ##
*###########################################################################

di as text ""
di as result "╔══════════════════════════════════════════════════════════════╗"
di as result "║  SECTION C1: FULL-CYCLE MULTIVARIATE VAR (TABLE 3 ROWS)     ║"
di as result "║  4 endogenous: MR, LA vol, LV vol, ALPM (z full-cycle)     ║"
di as result "║  Lags: 1, 2, 4, 7                                          ║"
di as result "╚══════════════════════════════════════════════════════════════╝"

foreach grp in afmr vfmr {
    foreach lag in 1 2 4 7 {

    use `fc_data', clear

    if "`grp'" == "afmr"  keep if group_fc == 1
    if "`grp'" == "vfmr"  keep if group_fc == 0

    * Patient dummies (fixed effects)
    capture drop pt_dum*
    quietly tabulate patient_id, gen(pt_dum)
    drop pt_dum1
    unab dumlist : pt_dum*

    * Time-set: use sequential index (time_index_combined_s has NaN for some
    * diastolic frames; using _n after sort ensures all observations are included)
    sort patient_id time_index_patient
    capture drop _t
    gen _t = _n
    tsset _t

    local n_obs = _N
    local n_pts = 0
    quietly tab patient_id
    local n_pts = r(r)

    di as text ""
    di as result "════════════════════════════════════════"
    di as result "  `grp' — FULL CYCLE — lag `lag' — N obs = `n_obs', N patients = `n_pts'"
    di as result "════════════════════════════════════════"

    * --- Fit VAR ---
    capture noisily var `fc_varlist', lags(1/`lag') exog(`dumlist')
    if _rc != 0 {
        di as error "  VAR failed for `grp' lag `lag'"
        continue
    }

    * --- Stability ---
    di as text "--- Stability ---"
    varstable

    * --- Granger causality: each predictor -> MR ---
    di as text ""
    di as text "--- Granger causality -> MR area (forward) ---"
    foreach pred in z_lavol_fc z_lvvol_fc z_alpm_fc {
        local terms ""
        forvalues l = 1/`lag' {
            local terms "`terms' [z_mr_fc]L`l'.`pred'"
        }
        capture test `terms'
        if _rc == 0 {
            di as text "    `pred' -> z_mr_fc: chi2 = " %8.3f r(chi2) "  p = " %6.4f r(p)
        }
    }

    * --- Joint Granger test (ALL -> MR) ---
    local all_terms ""
    foreach pred in z_lavol_fc z_lvvol_fc z_alpm_fc {
        forvalues l = 1/`lag' {
            local all_terms "`all_terms' [z_mr_fc]L`l'.`pred'"
        }
    }
    capture test `all_terms'
    if _rc == 0 {
        di as text "    ALL -> z_mr_fc (joint): chi2 = " %8.3f r(chi2) "  p = " %6.4f r(p)
    }

    * --- Granger causality: each predictor -> ALPM (reverse) ---
    di as text ""
    di as text "--- Granger causality -> ALPM (reverse direction) ---"
    foreach pred in z_mr_fc z_lavol_fc z_lvvol_fc {
        local terms ""
        forvalues l = 1/`lag' {
            local terms "`terms' [z_alpm_fc]L`l'.`pred'"
        }
        capture test `terms'
        if _rc == 0 {
            di as text "    `pred' -> z_alpm_fc: chi2 = " %8.3f r(chi2) "  p = " %6.4f r(p)
        }
    }

    * --- IRF/FEVD ---
    local irf_name = "fc_`grp'_lag`lag'"
    capture irf create `irf_name', step(10) set(`irf_name') replace nose
    if _rc == 0 {

        * FEVD for MR equation at horizon 10
        di as text ""
        di as text "--- FEVD for MR area at horizon 10 ---"
        irf table fevd, ///
            impulse(z_mr_fc z_lavol_fc z_lvvol_fc z_alpm_fc) ///
            response(z_mr_fc) set(`irf_name') noci

        * FEVD for ALPM equation at horizon 10
        di as text ""
        di as text "--- FEVD for ALPM at horizon 10 ---"
        irf table fevd, ///
            impulse(z_mr_fc z_lavol_fc z_lvvol_fc z_alpm_fc) ///
            response(z_alpm_fc) set(`irf_name') noci

        * OIRF at step 1 for MR equation
        di as text ""
        di as text "--- OIRF at step 1 (MR equation) ---"
        irf table oirf, ///
            impulse(z_lavol_fc z_lvvol_fc z_alpm_fc) ///
            response(z_mr_fc) set(`irf_name') noci

        * OIRF at step 1 for ALPM equation
        di as text ""
        di as text "--- OIRF at step 1 (ALPM equation) ---"
        irf table oirf, ///
            impulse(z_mr_fc z_lavol_fc z_lvvol_fc) ///
            response(z_alpm_fc) set(`irf_name') noci
    }

    }  /* end lag loop */
}  /* end grp loop */


*###########################################################################
*##                                                                       ##
*##    SECTION C2: EXTENDED LAG SWEEP (lags 1-20)                          ##
*##               ALPM -> MR and MR -> ALPM p-values only                 ##
*##               For results text (significance counts)                  ##
*##                                                                       ##
*###########################################################################

di as text ""
di as result "╔══════════════════════════════════════════════════════════════╗"
di as result "║  SECTION C2: LAG SWEEP 1-20 (p-values for text)             ║"
di as result "╚══════════════════════════════════════════════════════════════╝"

foreach grp in afmr vfmr {

    di as text ""
    di as result "────────────────────────────────────────"
    di as result "  `grp' — lag sweep 1-20"
    di as result "────────────────────────────────────────"
    di as text ""
    di as text "  Lag  N_obs  N_pts  ALPM->MR_p     MR->ALPM_p     LV->MR_p       LA->MR_p"
    di as text "  ───  ─────  ─────  ──────────     ──────────     ──────────     ──────────"

    forvalues lag = 1/20 {

        use `fc_data', clear

        if "`grp'" == "afmr"  keep if group_fc == 1
        if "`grp'" == "vfmr"  keep if group_fc == 0

        capture drop pt_dum*
        quietly tabulate patient_id, gen(pt_dum)
        drop pt_dum1
        unab dumlist : pt_dum*

        sort patient_id time_index_patient
        capture drop _t
        gen _t = _n
        tsset _t

        local n_obs = _N
        quietly tab patient_id
        local n_pts = r(r)

        capture noisily var `fc_varlist', lags(1/`lag') exog(`dumlist')
        if _rc != 0 {
            di as text "  " %3.0f `lag' "  " %5.0f `n_obs' "  " %5.0f `n_pts' "  VAR FAILED"
            continue
        }

        * ALPM -> MR
        local alpm_mr_terms ""
        forvalues l = 1/`lag' {
            local alpm_mr_terms "`alpm_mr_terms' [z_mr_fc]L`l'.z_alpm_fc"
        }
        local p_alpm_mr = .
        capture test `alpm_mr_terms'
        if _rc == 0  local p_alpm_mr = r(p)

        * MR -> ALPM
        local mr_alpm_terms ""
        forvalues l = 1/`lag' {
            local mr_alpm_terms "`mr_alpm_terms' [z_alpm_fc]L`l'.z_mr_fc"
        }
        local p_mr_alpm = .
        capture test `mr_alpm_terms'
        if _rc == 0  local p_mr_alpm = r(p)

        * LV -> MR
        local lv_mr_terms ""
        forvalues l = 1/`lag' {
            local lv_mr_terms "`lv_mr_terms' [z_mr_fc]L`l'.z_lvvol_fc"
        }
        local p_lv_mr = .
        capture test `lv_mr_terms'
        if _rc == 0  local p_lv_mr = r(p)

        * LA -> MR
        local la_mr_terms ""
        forvalues l = 1/`lag' {
            local la_mr_terms "`la_mr_terms' [z_mr_fc]L`l'.z_lavol_fc"
        }
        local p_la_mr = .
        capture test `la_mr_terms'
        if _rc == 0  local p_la_mr = r(p)

        di as text "  " %3.0f `lag' "  " %5.0f `n_obs' "  " %5.0f `n_pts' ///
            "  " %10.4f `p_alpm_mr' ///
            "     " %10.4f `p_mr_alpm' ///
            "     " %10.4f `p_lv_mr' ///
            "     " %10.4f `p_la_mr'
    }
}


*###########################################################################
*##                                                                       ##
*##    SECTION C3: ALPM VARIANCE DECOMPOSITION SUMMARY                     ##
*##               FEVD of ALPM equation at horizon 10                     ##
*##               Key lags only (1, 4, 10, 20)                            ##
*##                                                                       ##
*###########################################################################

di as text ""
di as result "╔══════════════════════════════════════════════════════════════╗"
di as result "║  SECTION C3: ALPM FEVD SUMMARY (for Discussion text)        ║"
di as result "╚══════════════════════════════════════════════════════════════╝"

foreach grp in afmr vfmr {
    foreach lag in 1 4 10 20 {

    use `fc_data', clear

    if "`grp'" == "afmr"  keep if group_fc == 1
    if "`grp'" == "vfmr"  keep if group_fc == 0

    capture drop pt_dum*
    quietly tabulate patient_id, gen(pt_dum)
    drop pt_dum1
    unab dumlist : pt_dum*

    sort patient_id time_index_patient
    capture drop _t
    gen _t = _n
    tsset _t

    local n_obs = _N
    quietly tab patient_id
    local n_pts = r(r)

    di as text ""
    di as result "  `grp' — lag `lag' — N = `n_obs' (`n_pts' pts)"

    capture noisily var `fc_varlist', lags(1/`lag') exog(`dumlist')
    if _rc != 0 {
        di as error "  VAR failed"
        continue
    }

    local irf_name = "fevd_`grp'_lag`lag'"
    capture irf create `irf_name', step(10) set(`irf_name') replace nose
    if _rc == 0 {
        di as text "  FEVD for ALPM equation at horizon 10:"
        irf table fevd, ///
            impulse(z_mr_fc z_lavol_fc z_lvvol_fc z_alpm_fc) ///
            response(z_alpm_fc) set(`irf_name') noci
    }

    }
}


restore     /* back to full dataset from full-cycle preserve */

*───────────────────────────────────────────────────────────────
*  CLEANUP AND FINISH
*───────────────────────────────────────────────────────────────
timer off 1
timer list

di as text ""
di as result "╔══════════════════════════════════════════════════╗"
di as result "║  ANALYSIS COMPLETE                              ║"
di as result "║                                                 ║"
di as result "║  Output files:                                  ║"
di as result "║    Table2_individual_VAR.csv / .xlsx             ║"
di as result "║    Lag_specific_results.csv                      ║"
di as result "║    Figure5_phenotypes.pdf / .png                 ║"
di as result "║    Figure6_lag_specific.pdf / .png               ║"
di as result "║    FMR_Granger_analysis.log                      ║"
di as result "║                                                 ║"
di as result "║  Table 3 results are printed in the log.        ║"
di as result "║  Full-cycle ALPM results are printed in the log.║"
di as result "║                                                 ║"
di as result "║  Fixed lag 3 primary for Table 2 (Section A)              ║"
di as result "╚══════════════════════════════════════════════════╝"

log close
