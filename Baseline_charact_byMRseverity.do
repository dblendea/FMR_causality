*###########################################################################
*  Table 1. Baseline Characteristics by MR Severity — Mild vs Moderate–Severe
*
*  Database: MRV_horizontal_10mar26.csv
*  Design:   3-panel table (AFMR, VFMR, Overall) × (Mild vs Mod–Severe)
*  Tests:    t-test (continuous), Fisher's exact (categorical)
*###########################################################################

clear all
set more off

log using "Table1_baseline.log", replace text

*───────────────────────────────────────────────────────────────
* 0.  LOAD DATA
*───────────────────────────────────────────────────────────────
import delimited "MRV_horizontal_10mar26.csv", clear varnames(1) encoding("utf-8")

* Drop trailing empty rows
drop if missing(id)

* Force numeric on key columns
destring age genderm0f1 mratrial1ventricular0 cad01 htn01 dm01 ///
    obesity01 sleepapnea01 smoking01 dyslipidemia01 nyhaclass ///
    afibhistory pacemaker01 crt01 hrbpm qrsms lbbb ///
    mr mr0mild1modsevere tr phtn ///
    lveddmm lvesdmm lvedvml ef ivs pwlv ///
    lvapex_base_lenght_earlysys_cm lvapex_base_length_endsys_cm ///
    la lavolume ///
    mv_diam_end_sys_cm_4ch mv_diam_end_diast ///
    mv_tentheight_earlysys_cm_psal mv_tenting_height_end_syst_psal ///
    mv_tentarea_earlysys_cm2_psal mv_tenting_area_end_syst_psal ///
    alpmlength_4ch_cm pmpmlength_pslongaxis_cm, ///
    replace force

* Rename long variable names to stay within Stata's 32-char macro limit
rename lvapex_base_lenght_earlysys_cm lv_len_earlysys
rename lvapex_base_length_endsys_cm   lv_len_endsys
rename mv_diam_end_sys_cm_4ch         mv_ann_sys
rename mv_diam_end_diast              mv_ann_diast
rename mv_tentheight_earlysys_cm_psal tent_ht_early
rename mv_tenting_height_end_syst_psal tent_ht_endsys
rename mv_tentarea_earlysys_cm2_psal  tent_area_early
rename mv_tenting_area_end_syst_psal  tent_area_endsys
rename alpmlength_4ch_cm              alpm_4ch
rename pmpmlength_pslongaxis_cm       pmpm_ps

* Verify
di as result "  N = " _N " patients"
tab mratrial1ventricular0

* Create readable group variables
gen subtype = cond(mratrial1ventricular0 == 1, "AFMR", "VFMR")
gen severity = mr0mild1modsevere
label define sev_lbl 0 "Mild" 1 "Mod-Severe"
label values severity sev_lbl

di as text ""
di as result "============================================"
di as result "  TABLE 1: BASELINE CHARACTERISTICS"
di as result "============================================"


*───────────────────────────────────────────────────────────────
* 1.  DEFINE VARIABLE LISTS
*───────────────────────────────────────────────────────────────

* Continuous variables: varname | label
* Categorical variables (binary 0/1): varname | label

* Demographics — continuous
local cont_demo age nyhaclass hrbpm qrsms
global lab_age         "Age (years)"
global lab_nyhaclass   "NYHA class"
global lab_hrbpm       "Heart rate (bpm)"
global lab_qrsms       "QRS duration (ms)"

* Demographics — categorical
local cat_demo genderm0f1 cad01 htn01 dm01 obesity01 dyslipidemia01 ///
    sleepapnea01 smoking01 afibhistory pacemaker01 crt01 lbbb
global lab_genderm0f1      "Female sex"
global lab_cad01           "CAD"
global lab_htn01           "Hypertension"
global lab_dm01            "Diabetes"
global lab_obesity01       "Obesity"
global lab_dyslipidemia01  "Dyslipidemia"
global lab_sleepapnea01    "Sleep apnea"
global lab_smoking01       "Smoking"
global lab_afibhistory     "AF history"
global lab_pacemaker01     "Pacemaker"
global lab_crt01           "CRT"
global lab_lbbb            "LBBB"

* Note: genderm0f1 is male=0, female=1, so the count of "1" = female count

* MR Characteristics — continuous
local cont_mr mr tr
global lab_mr  "MR grade (1-4)"
global lab_tr  "TR grade"

* MR Characteristics — categorical
local cat_mr phtn
global lab_phtn "Pulmonary hypertension"

* LV Parameters — continuous
local cont_lv lveddmm lvesdmm lvedvml ef ///
    lv_len_earlysys lv_len_endsys ///
    ivs pwlv
global lab_lveddmm        "LVEDD (mm)"
global lab_lvesdmm        "LVESD (mm)"
global lab_lvedvml        "LVEDV (ml)"
global lab_ef             "Ejection fraction (%)"
global lab_lv_len_earlysys "LV length early sys (cm)"
global lab_lv_len_endsys  "LV length end sys (cm)"
global lab_ivs            "IVS thickness (mm)"
global lab_pwlv           "PW thickness (mm)"

* LA Parameters — continuous
local cont_la la lavolume
global lab_la       "LA diameter (mm)"
global lab_lavolume "LA volume (ml)"

* MV Geometry — continuous
local cont_mv mv_ann_sys mv_ann_diast ///
    tent_ht_early tent_ht_endsys ///
    tent_area_early tent_area_endsys
global lab_mv_ann_sys      "MV annulus end-sys (cm)"
global lab_mv_ann_diast    "MV annulus end-diast (cm)"
global lab_tent_ht_early   "Tenting ht early sys (cm)"
global lab_tent_ht_endsys  "Tenting ht end-sys (cm)"
global lab_tent_area_early "Tenting area early (cm2)"
global lab_tent_area_endsys "Tenting area end-sys (cm2)"

* PM Length — continuous
local cont_pm alpm_4ch pmpm_ps
global lab_alpm_4ch  "ALPM length 4-ch (cm)"
global lab_pmpm_ps   "PMPM length PS (cm)"


*───────────────────────────────────────────────────────────────
* 2.  HELPER PROGRAMS
*───────────────────────────────────────────────────────────────

capture program drop table1_continuous
program define table1_continuous
    * Usage: table1_continuous varname
    args varname

    local lbl "${lab_`varname'}"
    
    quietly summarize `varname' if severity == 0 & $cond
    local m0 : di %4.1f r(mean)
    local s0 : di %4.1f r(sd)
    
    quietly summarize `varname' if severity == 1 & $cond
    local m1 : di %4.1f r(mean)
    local s1 : di %4.1f r(sd)
    
    * t-test
    capture ttest `varname' if $cond, by(severity)
    if _rc == 0 {
        local pval : di %5.3f r(p)
        if r(p) < 0.001  local pval "<0.001"
    }
    else {
        local pval "."
    }
    
    di as text "  `lbl'" _col(35) "`m0' (`s0')" _col(55) "`m1' (`s1')" _col(75) "`pval'"
end

capture program drop table1_categorical
program define table1_categorical
    * Usage: table1_categorical varname
    args varname
    
    local lbl "${lab_`varname'}"
    
    * Count value==1 in each severity group
    quietly count if `varname' == 1 & severity == 0 & $cond
    local n1_mild = r(N)
    quietly count if severity == 0 & $cond & !missing(`varname')
    local ntot_mild = r(N)
    local pct_mild : di %3.0f (`n1_mild'/`ntot_mild'*100)
    
    quietly count if `varname' == 1 & severity == 1 & $cond
    local n1_ms = r(N)
    quietly count if severity == 1 & $cond & !missing(`varname')
    local ntot_ms = r(N)
    local pct_ms : di %3.0f (`n1_ms'/`ntot_ms'*100)
    
    * Fisher's exact test
    * Handle case where all values are identical (e.g., 0 vs 0) — tab fails
    capture tab `varname' severity if $cond, exact
    if _rc == 0 {
        local pval : di %5.3f r(p_exact)
        if r(p_exact) < 0.001  local pval "<0.001"
    }
    else {
        * If tab fails (no variation), report p=1.000
        local pval "1.000"
    }
    
    di as text "  `lbl'" _col(35) "`n1_mild' (`pct_mild'%)" _col(55) "`n1_ms' (`pct_ms'%)" _col(75) "`pval'"
end


*───────────────────────────────────────────────────────────────
* 3.  PRINT TABLE — THREE PANELS
*───────────────────────────────────────────────────────────────

foreach panel in "AFMR" "VFMR" "Overall" {

    if "`panel'" == "AFMR"    global cond "mratrial1ventricular0 == 1"
    if "`panel'" == "VFMR"    global cond "mratrial1ventricular0 == 0"
    if "`panel'" == "Overall" global cond "1 == 1"
    
    * Count per severity
    quietly count if severity == 0 & $cond
    local n_mild = r(N)
    quietly count if severity == 1 & $cond
    local n_ms = r(N)
    quietly count if $cond
    local n_tot = r(N)
    
    di as text ""
    di as result "═══════════════════════════════════════════════════════════════════════════"
    di as result "  `panel' (n=`n_tot')    Mild (n=`n_mild')    Mod-Sev (n=`n_ms')"
    di as result "═══════════════════════════════════════════════════════════════════════════"
    di as text _col(35) "Mild" _col(55) "Mod-Sev" _col(75) "p"
    di as text "  ─────────────────────────────────────────────────────────────────────────"
    
    * --- Demographics ---
    di as result "  Demographics"
    foreach v of local cont_demo {
        table1_continuous `v'
    }
    foreach v of local cat_demo {
        table1_categorical `v'
    }
    
    * --- MR Characteristics ---
    di as text ""
    di as result "  MR Characteristics"
    foreach v of local cont_mr {
        table1_continuous `v'
    }
    foreach v of local cat_mr {
        table1_categorical `v'
    }
    
    * --- LV Parameters ---
    di as text ""
    di as result "  Left Ventricular Parameters"
    foreach v of local cont_lv {
        table1_continuous `v'
    }
    
    * --- LA Parameters ---
    di as text ""
    di as result "  Left Atrial Parameters"
    foreach v of local cont_la {
        table1_continuous `v'
    }
    
    * --- MV Geometry ---
    di as text ""
    di as result "  Mitral Valve Geometry"
    foreach v of local cont_mv {
        table1_continuous `v'
    }
    
    * --- PM Length ---
    di as text ""
    di as result "  Papillary Muscle Length"
    foreach v of local cont_pm {
        table1_continuous `v'
    }
}


*───────────────────────────────────────────────────────────────
* 4.  EXPORT TO CSV FOR EASY COMPARISON
*───────────────────────────────────────────────────────────────
di as text ""
di as result "============================================"
di as result "  EXPORTING TABLE 1 TO CSV"
di as result "============================================"

* All continuous variables in one list
local all_cont age nyhaclass hrbpm qrsms ///
    mr tr ///
    lveddmm lvesdmm lvedvml ef ///
    lv_len_earlysys lv_len_endsys ///
    ivs pwlv ///
    la lavolume ///
    mv_ann_sys mv_ann_diast ///
    tent_ht_early tent_ht_endsys ///
    tent_area_early tent_area_endsys ///
    alpm_4ch pmpm_ps

* All categorical variables
local all_cat genderm0f1 cad01 htn01 dm01 obesity01 dyslipidemia01 ///
    sleepapnea01 smoking01 afibhistory pacemaker01 crt01 lbbb phtn

tempfile t1_results
postfile t1_handle ///
    str40 variable str8 type ///
    str8 panel ///
    str20 mild str20 modsev ///
    double pval ///
    using `t1_results', replace

foreach panel in "AFMR" "VFMR" "Overall" {

    if "`panel'" == "AFMR"    local cond "mratrial1ventricular0 == 1"
    if "`panel'" == "VFMR"    local cond "mratrial1ventricular0 == 0"
    if "`panel'" == "Overall" local cond "1 == 1"
    
    * Continuous
    foreach v of local all_cont {
        quietly summarize `v' if severity == 0 & `cond'
        local m0 : di %4.1f r(mean)
        local s0 : di %4.1f r(sd)
        local mild_str = "`m0' (`s0')"
        
        quietly summarize `v' if severity == 1 & `cond'
        local m1 : di %4.1f r(mean)
        local s1 : di %4.1f r(sd)
        local ms_str = "`m1' (`s1')"
        
        capture ttest `v' if `cond', by(severity)
        local p = .
        if _rc == 0  local p = r(p)
        
        post t1_handle ("${lab_`v'}") ("cont") ("`panel'") ("`mild_str'") ("`ms_str'") (`p')
    }
    
    * Categorical
    foreach v of local all_cat {
        quietly count if `v' == 1 & severity == 0 & `cond'
        local n1 = r(N)
        quietly count if severity == 0 & `cond' & !missing(`v')
        local nt = r(N)
        local pct : di %3.0f (`n1'/`nt'*100)
        local mild_str = "`n1' (`pct'%)"
        
        quietly count if `v' == 1 & severity == 1 & `cond'
        local n1 = r(N)
        quietly count if severity == 1 & `cond' & !missing(`v')
        local nt = r(N)
        local pct : di %3.0f (`n1'/`nt'*100)
        local ms_str = "`n1' (`pct'%)"
        
        capture tab `v' severity if `cond', exact
        local p = 1
        if _rc == 0  local p = r(p_exact)
        
        post t1_handle ("${lab_`v'}") ("cat") ("`panel'") ("`mild_str'") ("`ms_str'") (`p')
    }
}

postclose t1_handle

use `t1_results', clear
export delimited using "Table1_baseline.csv", replace
list, noobs separator(0) abbreviate(40)

di as result "  Table 1 exported to Table1_baseline.csv"

log close
