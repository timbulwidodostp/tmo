*! version 0.9.0b3 2025-04-22

capture program drop tmo
program define tmo, rclass
    version 13

    syntax, ///
        cmd(str) x(varname) ylist(varlist) Idvar(varname) ///
        [Timevar(varname)] ///
        [LATitude(varname)] [LONgitude(varname)] [DISTTHREShold(real 0)] [miles] ///
        [THREShold(real -9)] [thresholdoff] ///
        [MISSlimit(real 0.1)] ///
        [FILEsuffix(str)] ///
        [savedyad] [load(str)] ///
        [plotq] ///
        [plothist] [plothistnbins(int 10000)] ///
        [plotse] [saveplotseest] ///
        [saveest] ///
        [scpc_cmd(str)] [scpc_uncond]
   
    ********************************
    *** RUN AND SAVE CMD OPTIONS ***
    ********************************

    qui `cmd'
    local spec `e(cmd)'
    
    * Assert that cmd uses supported command
    if !inlist("`spec'","regress","reghdfe","ivreghdfe","ivreg2") {
            di as error "`spec' not supported"
            exit
    }

    * Store original results and options
    local y `e(depvar)'
    scalar beta = _b[`x']
    scalar se = _se[`x']
    scalar N_obs = e(N)
    scalar N_clust = e(N_clust)
    scalar df_r = e(df_r)
    cap drop __tmo_sample
    gen byte __tmo_sample = e(sample)

    * Extract clusters
    local cluster `e(clustvar)'
    
    * Extract absorb vars for (iv)reghdfe
    local absorb_vars `e(absvars)'

    * Extract weights
    if "`e(wexp)'"!="" {
        local weightvar = subinstr("`e(wexp)'","=","",1)
        local weightexp [`e(wtype)'`e(wexp)']
    }
    else {
        local weightvar
        local weightexp
    }

    *** END RUN AND SAVE CMD OPTIONS ***



    *************************
    *** RUN SCPC IF GIVEN ***
    *************************
    
    if "`scpc_cmd'"!="" {
        * Install edited version of scpc that stores critical values
        qui net install scpc, from("https://raw.githubusercontent.com/wjnkim/tmo/master/scpc_tmo") replace

        preserve
            qui keep if !missing(`longitude') & !missing(`latitude') & __tmo_sample

            rename `latitude' s_1
            rename `longitude' s_2 
            
            qui hashsort `idvar'
            qui `scpc_cmd'
            qui keep if e(sample)
            if "`scpc_uncond'"=="" {
                qui scpc, k(1) latlong
            }
            else {
                qui scpc, k(1) latlong uncond
            }
            scalar scpc_se=e(scpcstats)[1,2]
        
            mata: id_scpc = st_data(.,"`idvar'")
            clear

            mata: id_scpc_uniq = uniqrows(id_scpc)
            mata: Wfin_sum_vec = vec(Wfin[.,2::cols(Wfin)]*Wfin[.,2::cols(Wfin)]'):/(cols(Wfin)-1)
            mata: id_scpc_rowvec = vec(J(1,rows(Wfin),id_scpc_uniq))
            mata: id_scpc_colvec = vec(J(rows(Wfin),1,id_scpc_uniq'))
                
            mata: st_local("scpc_obsN",strofreal(rows(Wfin_sum_vec),"%50.0f"))

            gen id1=.
            gen id2=.
            gen Wfin=.
                
            qui set obs `scpc_obsN'
            
            mata: st_store(.,.,(id_scpc_rowvec,id_scpc_colvec,Wfin_sum_vec))
            
            qui keep if id1>=id2 // keep only lower triangular

            qui compress
            tempfile Wfin
            qui save `Wfin'
        restore
    }
    else {
        global scpc_cv = .
    }

    *** END RUN SCPC IF GIVEN ***



    *********************
    *** OPTION CHECKS ***
    *********************

    * Assert gtools installed
    cap which gtools
    if _rc {
        di as error "tmo requires gtools package -- please run: ssc install gtools"
        exit
    }

    * Assert that y and x are in cmd and y is the dependent variable
    local ycheck: word 2 of `cmd'
	if strpos("`cmd'","`y'")==0 | strpos("`cmd'","`x'")==0 | "`y'"!="`ycheck'" {
		di as error "cmd must contain `y' and `x' and `y' must be independent var"
		exit
	}

    * Assert that y appears only once in cmd
	if (strlen("`cmd'")-strlen(subinstr("`cmd'"," `y' ","",.)) != strlen(" `y' ")) {
		di as error "cmd contains multiple instances of `y'"
		exit
	}

    * Assert that x appears only once in cmd
	if (strlen("`cmd'")-strlen(subinstr("`cmd'","`x'","",.)) != strlen("`x'")) {
		di as error "cmd contains multiple instances of `x', please rename variables"
		exit
	}

    * Assert no duplicates in y ylist
    local done ""
    local dups ""
    foreach var in `y' `ylist' {
        confirm var `var', exact
        if strpos("`done'"," `var' ")>0 {
            local dups "`dups' `var'"
        }
        local done "`done' `var'"
    }
    if "`dups'"!="" {
        di as error "Duplicated variables in depvar/ylist: `dups'"
        exit
    }

    * Assert misslimit is between 0 and 1
    if `misslimit'<0 | `misslimit'>1 {
        di as error "misslimit() value must be between 0 and 1"
        exit
    }

    * Assert plothistnbins is positive if given
    if `plothistnbins'<=0 {
        di as error "plothistnbins() value must be positive"
        exit
    }

    * saveplotse requires plotse option
    if "`plotse'"=="" & "`saveplotse'"!="" {
        di as error "saveplotse option requires plotse option"
        exit
    }

    * Require file path for saving figures/data
    if "`filesuffix'"=="" & ("`plotse'"!="" | "`plothist'"!="" | "`savedyad'"!="" | "`saveplotseest'"!="" | "`saveest'"!="") {
        di as error "filesuffix() required for `plotse' `plothist' `savedyad' `saveplotseest' `saveest'"
    }

    * Require clustering (at least at location level) if panel
    if "`timevar'"!="" & "`cluster'"=="" {
        di as error "Clustering required for panel case (at least at `id' level)"
        exit
    }

    * Require absorb if (iv)reghdfe
    if inlist("`spec'","reghdfe","ivreghdfe") & "`absorb_vars'"=="" {
        di as error "absorb() required for `spec'"
        exit
    } 

    * Assert idvar is unique (within timevar if panel)
	if "`idvar'"!="" & "`timevar'"=="" {
		gisid `idvar' if __tmo_sample
	}
	if "`idvar'"!="" & "`timevar'"!="" {
		gisid `idvar' `timevar' if __tmo_sample
	}

    * Assert longitude and latitude provided if distthreshold!=0
    if ("`longitude'"=="" | "`latitude'"=="") & `distthreshold'!=0 {
        di as error "longitude() and latitude() required for distthreshold() option"
        exit
    }

    * Assert longitude and latitude provided if scpc
    if ("`longitude'"=="" | "`latitude'"=="") & "`scpc_cmd'"!="" {
        di as error "longitude() and latitude() required for scpc_cmd() option"
        exit
    }

    * Assert scpc_cmd() if scpc_uncond option
    if "`scpc_cmd'"=="" & "`scpc_uncond'"!="" {
        di as error "scpc_cmd() required for scpc_uncond() option"
        exit
    }

    * Assert no weights if scpc
    if "`scpc_cmd'"!="" & "`weightvar'"!="" {
        di as error "weights not allowed for scpc"
        exit
    }

    * Assert distthreshold>0 if provided
    if `distthreshold'<0 {
        di as error "distthreshold() must be greater than 0"
        exit
    }

    * Assert geodist package installed if distthreshold>0
    if `distthreshold'>0 {
        cap which geodist
        if _rc {
            di as error "distthreshold() requires geodist package -- please run: ssc install geodist"
            exit
        }
    }

    * Assert custom threshold is between 0 and 1 if provided
    if `threshold'!=-9 & (`threshold'<0 | `threshold'>1) {
        di as error "Custom threshold() must be between 0 and 1"
        exit
    }

    * Store number of locations and time periods
    qui gdistinct `idvar' if __tmo_sample
    local N = r(ndistinct)
    scalar N = `r(ndistinct)'
    mata: N = `N'
    if "`timevar'"!="" {
        qui gdistinct `timevar' if __tmo_sample
        local T = r(ndistinct)
        scalar T = r(ndistinct)
        mata: T = `T'
    }
    else {
        local T = 1
        scalar T = 1
        mata: T = 1
    }

    *** END OPTION CHECKS ***



    *****************
    *** CLEAN CMD ***
    *****************
    
    if inlist("`spec'","reghdfe","ivreghdfe") {
        * Remove any resid specified in cmd already
        local comma_start = strpos("`cmd'", ",")
        local before_comma = substr("`cmd'", 1, `comma_start'-1)
        local after_comma = substr("`cmd'", `comma_start'+1, .)
        local after_comma = regexr("`after_comma'", "(res|resi|resid|residu|residua|residual|residuals)\([^)]*\)", "")
        local after_comma = regexr("`after_comma'", "(res|resi|resid|residu|residua|residual|residuals)", "")

        * Remove saving FE options
        local after_comma = regexr("`after_comma'", "(a|ab|abs|abso|absor|absorb)\([^)]*\)", "")

        local cmd `before_comma', absorb(`absorb_vars') `after_comma'
    }

    *** END CLEAN CMD ***



    ****************************
    *** STORE ID-CLUSTER XW ****
    ****************************

    if "`cluster'"!="" & "`cluster'"!="`idvar'" {
        preserve
            qui keep if __tmo_sample
            keep `idvar' `cluster'

            qui gduplicates drop
            cap gisid `idvar'

            if _rc {
                di as error "Only clustering by groups of locations is supported (`cluster' must be constant within `idvar')"
                exit
            }

            local cln=1
            foreach var in `cluster' {
                rename `var' __tmo_cl1_`cln'
                local ++cln
            }
            rename `idvar' id1
            qui compress
            tempfile cl1
            qui save `cl1'

            local cln=1
            foreach var in `cluster' {
                rename __tmo_cl1_`cln' __tmo_cl2_`cln'
                local ++cln
            }
            rename id1 id2
            tempfile cl2
            qui save `cl2'
        restore
    }

    *** END STORE ID-CLUSTER XW ***



    *********************************
    *** STORE LOCATION DISTANCES ****
    *********************************

    if "`longitude'"!="" & "`latitude'"!="" {
        preserve
            qui keep if __tmo_sample
            keep `idvar' `longitude' `latitude'

            qui gduplicates drop
            cap gisid `idvar'

            if _rc {
                di as error "`longitude' and/or `latitude' not constant within some `idvar'"
                exit
            }

            gen n=_n
		
            rename `latitude' lat2
            rename `longitude' lon2
            rename `idvar' id2
            qui compress

            tempfile dist
            qui save `dist'
            
            rename id2 id1
            rename lat2 lat1
            rename lon2 lon1	
            
            qui sum n
            qui expand `r(max)'
            drop n
            
            qui hashsort id1
            by id1: gen n=_n
            
            qui merge m:1 n using `dist', assert(3) nogen
            drop n

            qui keep if id1>=id2 // keep only lower triangular
            
            qui geodist lat1 lon1 lat2 lon2, gen(dist) `miles' sphere
            qui replace dist=0 if id1==id2

            keep id1 id2 dist

            qui compress   
            tempfile dist
            qui save `dist'
        restore
    }

    *** END STORE LOCATION DISTANCES ***



    ****************************
    *** WRITE CMD FOR XTILDE ***
    ****************************

    if  inlist("`spec'","regress","reghdfe") {
        local xtildecmd = subinstr("`cmd'"," `y' "," `x' ",1)
        local xtildecmd = subinstr("`xtildecmd'"," `x' "," ",2)
        if "`spec'"=="reghdfe" {
            local xtildecmd `xtildecmd' resid
        }
    }
    if "`spec'"=="ivreghdfe" {        
        * Create xtildecmd using reghdfe
        local parens_start = strpos("`cmd'", "(")
        local parens_end = strpos("`cmd'", ")")
        local paren_content = substr("`cmd'", `parens_start' + 1, `parens_end' - `parens_start' - 1)
        local instr_start = strpos("`paren_content'", "=")
        local endog_part = substr("`paren_content'", 1, `instr_start'-1)
        local endog = word("`endog_part'", 1)
        if wordcount("`endog_part'") > 1 {
            di as error "Multiple endogenous regressors not supported"
            exit
        }
        if "`endog'"!="`x'" {
            di as error "Endogenous regressor `endog' is not `x'"
            exit
        }
        local instr = substr("`paren_content'", `instr_start' + 1, .)

        local remaining = subinstr("`cmd'", "ivreghdfe", "", 1)
        local remaining = subinstr("`remaining'", "`y'", "", 1)
        local comma_start = strpos("`remaining'", ",")
        local after_comma = substr("`remaining'", `comma_start'+1, .)
        local remaining = substr("`remaining'", 1, `comma_start'-1)
        local paren_start = strpos("`remaining'", "(")
        local paren_end = strpos("`remaining'", ")")
        local controls = substr("`remaining'", 1, `paren_start'-1) + " " + substr("`remaining'", `paren_end'+1, .)
        local xtildecmd1 reghdfe `endog' `instr' `controls', `after_comma'
        local xtildecmd2 reghdfe __tmo_xhat `controls', `after_comma' resid
    }
    if "`spec'"=="ivreg2" {
        * Create xtildecmd using reg
        local parens_start = strpos("`cmd'", "(")
        local parens_end = strpos("`cmd'", ")")
        local paren_content = substr("`cmd'", `parens_start' + 1, `parens_end' - `parens_start' - 1)
        local instr_start = strpos("`paren_content'", "=")
        local endog_part = substr("`paren_content'", 1, `instr_start'-1)
        local endog = word("`endog_part'", 1)
        if wordcount("`endog_part'") > 1 {
            di as error "Multiple endogenous regressors not supported"
            exit
        }
        if "`endog'"!="`x'" {
            di as error "Endogenous regressor `endog' is not `x'"
            exit
        }
        local instr = substr("`paren_content'", `instr_start' + 1, .)

        local remaining = subinstr("`cmd'", "ivreg2", "", 1)
        local remaining = subinstr("`remaining'", "`y'", "", 1)
        local comma_start = strpos("`remaining'", ",")
        local remaining = substr("`remaining'", 1, `comma_start'-1)
        local paren_start = strpos("`remaining'", "(")
        local paren_end = strpos("`remaining'", ")")
        local controls = substr("`remaining'", 1, `paren_start'-1) + " " + substr("`remaining'", `paren_end'+1, .)

        local comma_start = strpos("`cmd'", ",")
        local after_comma = substr("`cmd'", `comma_start'+1, .)
        local has_nocon = regexm("`after_comma'", "(,noc|,noco|,nocon |,nocons|,noconst|,noconsta|,noconstan|,noconstant| noc| noco| nocon| nocons| noconst| noconsta| noconstan| noconstant)")
        if `has_nocon' local nocons , nocons 
        else local nocons

        local xtildecmd1 reg `endog' `instr' `controls' `nocons'
        local xtildecmd2 reg __tmo_xhat `controls' `nocons'
    }
    
    *** END WRITE CMD FOR XTILDE ***

    

    ***************
    *** RUN TMO ***
    ***************
    preserve
        if "`load'"=="" { // if dyad data already exists, can load and skip this part (programmer option)
            qui keep if __tmo_sample

            * Estimate __tmo_xtilde
            if  inlist("`spec'","regress","reghdfe") {
                qui `xtildecmd'
                qui predict __tmo_xtilde, resid
            }
            if inlist("`spec'","ivreghdfe","ivreg2") {
                qui `xtildecmd1'
                qui predict __tmo_xhat
                qui `xtildecmd2'
                qui predict __tmo_xtilde, resid
                drop __tmo_xhat
            }

            * Check __tmo_xtilde is correct
            qui reg `y' __tmo_xtilde `weightexp'
            if abs(beta-_b[__tmo_xtilde])>1e-5 {
                di as error "__tmo_xtilde is incorrect"
                exit
            }

            * Loop through auxiliary outcomes and save residuals
            cap drop __tmo_resid*

            if  inlist("`spec'","regress","ivreg2") {
                * Only keep cmd before comma (faster runtime)
                local comma_start = strpos("`cmd'",",")
                if `comma_start'>0 {
                    local before_comma = substr("`cmd'",1,`comma_start'-1)
                    local after_comma = substr("`cmd'", `comma_start'+1, .)

                    * Check whether there is nocons option and include if so
                    local has_nocon = regexm("`after_comma'", "(,noc|,noco|,nocon |,nocons|,noconst|,noconsta|,noconstan|,noconstant| noc| noco| nocon| nocons| noconst| noconsta| noconstan| noconstant)")
                    if `has_nocon' {
                        local nocons , nocons 
                    }
                    else local nocons
                }
                else {
                    local before_comma `cmd'
                    local nocons
                }
                
                local cmd_toloop `before_comma' `nocons'

                local ynum=1
                foreach aux_y in `y' `ylist' {
                    local cmd_inloop = subinstr("`cmd_toloop'","`y'","`aux_y'",1)
                    qui `cmd_inloop'
                    qui predict __tmo_resid`ynum', resid
                    local ++ynum
                }
            }
            
            if "`spec'"=="reghdfe" {
                * Remove any clustering (faster runtime)
                local cmd_toloop `cmd'
                while regexm("`cmd_toloop'", "(cl|clu|clus|clust|cluste|cluster|vce)\([^)]*\)") {
                    local cmd_toloop = regexr("`cmd_toloop'", "(cl|clu|clus|clust|cluste|cluster|vce)\([^)]*\)", "")
                }

                local ynum=1
                foreach aux_y in `y' `ylist' {
                    local cmd_inloop = subinstr("`cmd_toloop'","`y'","`aux_y'",1)
                    qui `cmd_inloop' resid
                    qui predict __tmo_resid`ynum', resid
                    local ++ynum
                }
            }

            if "`spec'"=="ivreghdfe" {
                * Remove any clustering (faster runtime)
                local cmd_toloop `cmd'
                while regexm("`cmd_toloop'", "(cl|clu|clus|clust|cluste|cluster|vce)\([^)]*\)") {
                    local cmd_toloop = regexr("`cmd_toloop'", "(cl|clu|clus|clust|cluste|cluster|vce)\([^)]*\)", "")
                }

                * Specify FEs to save
                local absorb_start = regexm("`cmd'", "(a|ab|abs|abso|absor|absorb)\(([^)]+)\)")
                local absorb_vars "`=regexs(2)'"
                local absorb_vars_savefe
                local absorb_vars_fesum
                local fe=1
                foreach fe_var in `absorb_vars' {
                    local absorb_vars_savefe `absorb_vars_savefe' __tmo_fe`fe'=`fe_var'
                    local absorb_vars_fesum `absorb_vars_fesum' + __tmo_fe`fe'
                    local ++fe
                }
                local cmd_toloop = regexr("`cmd_toloop'", "(a|ab|abs|abso|absor|absorb)\([^)]*\)", "")
                local cmd_toloop `cmd_toloop' absorb(`absorb_vars_savefe')

                local ynum=1
                foreach aux_y in `y' `ylist' {
                    local cmd_inloop = subinstr("`cmd_toloop'","`y'","`aux_y'",1)
                    cap drop __tmo_fe*
                    qui `cmd_inloop'
                    cap drop __tmo_xb 
                    qui predict __tmo_xb
                    qui gen __tmo_resid`ynum' = `aux_y' - (__tmo_xb`absorb_vars_fesum')
                    local ++ynum
                }
            }
            scalar D = `ynum'-1

            * Calculate correlation in residuals and contribution to variance
            keep `idvar' `timevar' `weightvar' __tmo_xtilde __tmo_resid*
            qui hashsort `idvar' `timevar'

            if "`timevar'"=="" { // Cross-sectional case
                * Store data in Mata
                mata: id = st_data(.,"`idvar'")
                mata: xtilde = st_data(.,"__tmo_xtilde")

                if "`weightvar'"!="" {
                    mata: wgt = st_data(.,"`weightvar'")
                }
                else {
                    mata: wgt = J(rows(xtilde),1,1)
                }
                mata: xtilde_wgt = xtilde:*wgt
                
                mata: Res1 = st_data(.,"__tmo_resid1")
                mata: Res = st_data(.,"__tmo_resid*")

                * Compute correlation in outcomes
                clear
                mata: CovEpsVec = DenomVec = LowTriInd = id_widerowvec = id_widecolvec = J(0,0,.)
                mata: corr_resid(id, Res, `misslimit', CovEpsVec, DenomVec, LowTriInd, id_widerowvec, id_widecolvec)

                * Compute contribution to SE for each pair of locations
                mata: ResXtildeVec = J(0,0,.)
                mata: sandwich_crosssec(id_widerowvec, id_widecolvec, LowTriInd, xtilde_wgt, Res1, ResXtildeVec)
            }
            else { // Panel case
                * Store data in Mata
                mata: id = st_data(.,"`idvar'")
                mata: t = st_data(.,"`timevar'")
                mata: xtilde = st_data(.,"__tmo_xtilde")
            
                if "`weightvar'"!="" {
                    mata: wgt = st_data(.,"`weightvar'")
                }
                else {
                    mata: wgt = J(rows(xtilde),1,1)
                }
                mata: xtilde_wgt = xtilde:*wgt

                mata: Res1 = st_data(.,"__tmo_resid1")
                
                ** Reshape Res
                * To make sure each location is shown for all time periods
                qui  tsset `idvar' `timevar'
                tsfill, full

                * Drop time periods that are all missing for resids of main outcome 
                qui gegen __tmo_resid1_missing = sum(!missing(__tmo_resid1)), by(`timevar')
                qui drop if __tmo_resid1_missing==0
                drop __tmo_resid1_missing
                
                * Input NT x D matrix of residuals to Mata
                mata: Res = st_data(.,"__tmo_resid*")
                mata: Dn = cols(Res)
                mata: idfull = st_data(.,"`idvar'")
                
                * Reshape matrix of residuals to N x TD
                mata: Res = rowshape(Res,N)
                mata: id_wide = rowshape(idfull,N)
                mata: id_wide = id_wide[,1]
                mata: assert(rows(uniqrows(id_wide))==N)
                mata: assert(cols(Res)==Dn*T)

                * Compute correlation in outcomes
                clear
                mata: CovEpsVec = DenomVec = LowTriInd = id_widerowvec = id_widecolvec = J(0,0,.)
                mata: corr_resid(id_wide, Res, `misslimit', CovEpsVec, DenomVec, LowTriInd, id_widerowvec, id_widecolvec)

                * Compute contribution to SE for each pair of locations
                sandwich_panel, rows(20000000) nloc(`N') ntime(`T') noi
            }

            * Normalize contribution to variance
            mata: resxtildenorm(ResXtildeVec, xtilde, xtilde_wgt)

            * Calculate part of contribution to scpc SE if specified
            if "`scpc_cmd'"!="" {
                mata: y_scpc = (sqrt(N)/(xtilde'*xtilde)):*Res1:*xtilde	
                mata: st_local("scpc_obsN",strofreal(rows(id),"%50.0f"))
                
                clear
                gen id1=.
                gen y_scpc1=.
                
                qui set obs `scpc_obsN'
                
                mata: st_store(.,.,(id,y_scpc))
                
                cap gisid id1
                if _rc gcollapse (sum) y_scpc1, by(id1)

                qui compress
                tempfile scpc1
                qui save `scpc1'
                
                rename id1 id2
                rename y_scpc1 y_scpc2
                
                tempfile scpc2
                qui save `scpc2'
            }

            * Bring Mata data into Stata
            if "`cluster'"!="" & "`cluster'"!="`idvar'" {
                local cl1p cl1path(`cl1')
                local cl2p cl2path(`cl2') 
            }
            else {
                local cl1p
                local cl2p
            }
            if "`filesuffix'"!="" {
                local savepath filesuffix(`filesuffix') 
            }
            else {
                local savepath
            }
            if "`longitude'"!="" & "`latitude'"!="" {
                local distp distpath(`dist')
            }
            else {
                local distp
            }
            if "`scpc_cmd'"!="" {
                local scpcw scpcwfin(`Wfin')
                local scpcy1 scpcy1path(`scpc1')
                local scpcy2 scpcy2path(`scpc2')
            }
            else {
                local scpcw
                local scpcy1
                local scpcy2
            }
            load_data, nloc(`N') `cl1p' `cl2p' `distp' `savedyad' `savepath' `scpcw' `scpcy1' `scpcy2'
        }
        else {
            if "`filesuffix'"!="" {
                local savepath filesuffix(`filesuffix') 
            }
            else {
                local savepath
            }

            use "`load'_est.dta", clear
            scalar D = N_outcomes[1]

            use "`load'_dyad.dta", clear
        }

        * Compute degrees of freedom and threshold
        dfqt, `plotq' `plothist' nbins(`plothistnbins') nloc(`N') `savepath'

        * Compute TMO SE
        est_tmo_se, dist_cutoff(`distthreshold') custom_thres(`threshold') `thresholdoff'

        * Plot TMO SE over threshold 
        if "`plotse'"!="" {
            tmo_over_thres, `savepath' `saveplotseest' noi dist_cutoff(`distthreshold')
        }
        


        **********************
        *** OUTPUT RESULTS ***
        **********************
        
        clear
        return clear
        ereturn clear
        
        mat tmo_results = J(1,6,.)
        mat rownames tmo_results = "`x'"
        mat colnames tmo_results = "Coef" "TMO SE" "t" "P>|t|" "95% Conf" "Interval"
        mat tmo_results[1,1] = beta
        mat tmo_results[1,2] = tmo_se
        mat tmo_results[1,3] = beta/tmo_se
        mat tmo_results[1,4] = 2*ttail(df_r, abs(beta/tmo_se))
        scalar lb = beta - invttail(df_r,0.025)*tmo_se
        scalar ub = beta + invttail(df_r,0.075)*tmo_se
        mat tmo_results[1,5] = lb
        mat tmo_results[1,6] = ub
        matlist tmo_results, border(all) cspec(o2& %20s | %9.3f o2 & %9.3f o2 & %6.2f o2 & %4.3f o2 & %9.3f o2 & %9.3f o2 &) rspec(&-&)

        mat tmo_details = J(5,1,.)
        mat rowname tmo_details = "Optimal threshold" "% of off-diag in SE est." "% >= threshold (excl. clusters/Conley)" "# outcomes" "Degrees of freedom" 
        mat tmo_details[1,1] = thres
        mat tmo_details[2,1] = offdP*100
        mat tmo_details[3,1] = offdPnocl*100
        mat tmo_details[4,1] = D
        mat tmo_details[5,1] = df
        
        matlist tmo_details, cspec(& %38s | %9.3f &) rspec(& & & & - & &) coleqonly

        return scalar beta = beta
        return scalar orig_se = se
        return scalar tmo_se = tmo_se
        return scalar lb = lb
        return scalar ub = ub
        return scalar threshold = thres
        return scalar pct_ge_thres = offdP*100
        return scalar pct_ge_thres_nocl = offdPnocl*100
        return scalar T = T
        return scalar N_loc = N
        return scalar N_clust = N_clust
        return scalar N_outcomes = D
        return scalar N = N_obs
        return scalar dof = df
        return scalar finite_sample_dof = dof_adj
        return scalar df_r = df_r
        return scalar scpc_cv = ${scpc_cv}
    restore

    if "`saveest'"!="" {
        tmo_save, `savepath'
    }
end



***********************
*** HELPER COMMANDS ***
***********************

* Function for computing correlation of residuals
cap mata: mata drop corr_resid()
mata 
	void corr_resid(real matrix id, real matrix Res, real scalar misslimit, real matrix CovEpsVec, real matrix DenomVec, real matrix LowTriInd, real matrix id_widerowvec, real matrix id_widecolvec)
	{
		real matrix Res_ms, Res_ms_lethres, Res_ms_ind, Res_no_ms, Denom, ResSum_DinBoth, ResMean_DinBoth, ResMeanProd_DinBoth, CovEps, Sum_ResSq_DinBoth, ResSD_DinBoth, DenomCorr
		real scalar demean

        // Drop outcomes that are missing for more than `misslimit'
        Res_ms = colsum(Res:==.)
        Res_ms_lethres = (Res_ms:/rows(Res)):<=misslimit
        Res = select(Res,Res_ms_lethres)

        // standardize residuals
        Res = Res:-J(rows(Res),1,colsum(Res):/colsum(Res:!=.))
        Res = Res:/J(rows(Res),1,(colsum(Res:^2):/colsum(Res:!=.)):^0.5) // Studentize residuals
        assert (sum(abs(colsum(Res)):<1e-5)==cols(Res))
        assert (sum(abs((colsum(Res:^2):/colsum(Res:!=.)):^0.5 :- 1):<1e-5)==cols(Res))

        // covariance of residuals -> correlations
        Res_ms_ind = Res:!=. // keep track of missing (0==missing)
        Res_no_ms = editmissing(Res,0) // Res with missing replaced with 0
        Denom = Res_ms_ind*Res_ms_ind' // number of both nonmissing
        
        ResSum_DinBoth = Res_no_ms * Res_ms_ind'
        ResMean_DinBoth = ResSum_DinBoth :/ Denom
        ResMeanProd_DinBoth = ResMean_DinBoth :* ResMean_DinBoth'
        
        demean=1 // set to 1 to make residuals mean 0 within each location
        if (demean==1) {
            CovEps = ((Res_no_ms*Res_no_ms'):/Denom) - (ResMean_DinBoth:*(ResSum_DinBoth':/Denom)) - (ResMean_DinBoth':*(ResSum_DinBoth:/Denom)) + ResMeanProd_DinBoth
        }
        else {
            CovEps = ((Res_no_ms*Res_no_ms'):/Denom)
        }
        
        Sum_ResSq_DinBoth = (Res_no_ms:^2) * Res_ms_ind'
        ResSD_DinBoth = ((Sum_ResSq_DinBoth:/Denom) + (-2:*ResMean_DinBoth:*(ResSum_DinBoth:/Denom)) + ResMean_DinBoth:^2) :^ 0.5
        
        DenomCorr = ResSD_DinBoth :* ResSD_DinBoth'
        
        if (demean==1) {
            CovEps = CovEps:/DenomCorr
        }
        else {
            CovEps = CovEps :/ ((diagonal(CovEps):^0.5) * (diagonal(CovEps):^0.5)')
        }

        // vectorize off-diag covariance of residuals
        CovEpsVec = vec(CovEps)
        LowTriInd = J(rows(CovEps),cols(CovEps),1)
        _lowertriangle(LowTriInd,1)
        LowTriInd = vec(LowTriInd)
        CovEpsVec = select(CovEpsVec,LowTriInd)
        assert (rows(CovEpsVec)==rows(Res)*(rows(Res)+1)/2)
        
        // vectorize Denom
        DenomVec = vec(Denom)
        DenomVec = select(DenomVec,LowTriInd)

        // store id vectors
        id_widerowvec = select(vec(J(1,rows(id),id)),LowTriInd)
	    id_widecolvec = select(vec(J(rows(id),1,id')),LowTriInd)
	}
end

* Function for computing contribution to SE for each pair of locations in cross-sectional case
cap mata: mata drop sandwich_crosssec()
mata
    void sandwich_crosssec(real matrix id_widerowvec, real matrix id_widecolvec, real matrix LowTriInd, real matrix xtilde, real matrix Res1, real matrix ResXtildeVec)
	{
		real matrix Xtilde, XtildeVec, ResVec
		
		// xx'	
		Xtilde = xtilde*xtilde'
		XtildeVec = vec(Xtilde)
			
		// off-diag product of residuals for main outcome
		ResVec = vec(Res1*Res1')
		
		// contribution to SE
		ResXtildeVec = ResVec:*XtildeVec
		
		// create lower diag incl diag
		ResXtildeVec = select(ResXtildeVec,LowTriInd)

        // with id numbers
        ResXtildeVec = id_widerowvec, id_widecolvec, ResXtildeVec
    }
end

* Functions for computing contribution to SE for each pair of locations in panel case
cap program drop sandwich_panel
program define sandwich_panel,
    syntax, rows(int) nloc(int) ntime(int) [NOIsily]

    clear
    mata: ResXtildeVec = J(0,3,.)

    * Loop through location-pairs for RAM purposes
    local itersize = ceil(`rows'/(`nloc'*`ntime'*`ntime'))
    forv i = 1 (`itersize') `nloc' {
        clear

        local start_i = `i'
		local end_i = min(`i'+`itersize'-1, `nloc')
        
        if "`noisily'"!="" {
            local donepct = string(`start_i'*100/`nloc', "%5.2f")
            di "Computed `donepct'% of sandwich"
        }
        
        mata: Ivec = Jvec = T1vec = T2vec = resx_i = J(0,0,.)
        mata: sandwich_panel_loop(`start_i',`end_i',id,t,id_wide,xtilde_wgt,Res1,Ivec,Jvec,T1vec,T2vec,resx_i)

        gen id1 = .
        gen id2 = .
        gen resxx = .
    
        qui set obs `obsN'
        
        mata: st_store(.,.,(Ivec,Jvec,resx_i))
        qui drop if id2>id1

        * Collapse contribution to SE across time
        gcollapse (sum) resxx, by(id1 id2)

        * Append to master
        mata: ResXtildeVec = ResXtildeVec \ st_data(.,"id1 id2 resxx")
    }
end

cap mata: mata drop sandwich_panel_loop()
mata
     void sandwich_panel_loop(real scalar start_i, real scalar end_i, real matrix id, real matrix t, real matrix id_wide, real matrix xtilde, real matrix Res1, real matrix Ivec, real matrix Jvec, real matrix T1vec, real matrix T2vec, real matrix resx_i)
	{
		real matrix ind_row, ind_col, xtilde_i, res_i

        // Index of I in NT rows
		ind_row = ((id:>=id_wide[start_i]) + (id:<=id_wide[end_i])) :== 2
				
		// Index of J<=I in NT cols
		ind_col = id:<=id_wide[end_i]
		
		// Select i's rows in xtilde and res for y of interest
		xtilde_i = select(xtilde,ind_row)
		res_i = select(Res1,ind_row)
			
		// Compute product with j cols where j<i
		xtilde_i = xtilde_i*select(xtilde,ind_col)'
        res_i = res_i*select(Res1,ind_col)'
		resx_i = xtilde_i:*res_i
		
		// Get id numbers
		Ivec = vec(J(1,cols(resx_i),select(id,ind_row)))
		Jvec = vec(J(rows(resx_i),1,select(id,ind_col)'))
		
		// Get time periods 
		T1vec = vec(J(1,cols(resx_i),select(t,ind_row)))
		T2vec = vec(J(rows(resx_i),1,select(t,ind_col)'))
		
		// Vectorize resx
		resx_i = vec(resx_i)
		st_local("obsN",strofreal(rows(resx_i),"%50.0f"))
    }
end

* Function to normalize contribution to SE by 2*[X'X]^-1 and multiply off-diag by 2 (since only lower triangular)
cap mata: mata drop resxtildenorm()
mata
    void resxtildenorm(real matrix ResXtildeVec, real matrix xtilde, real matrix xtilde_wgt)
	{
		real matrix denom, offdiag

        denom = 1/(xtilde_wgt'*xtilde)^2
		offdiag = ResXtildeVec[,1]:!=ResXtildeVec[,2]

        ResXtildeVec[,3] = ResXtildeVec[,3]:*denom
        ResXtildeVec[,3] = ResXtildeVec[,3] + (ResXtildeVec[,3]:*offdiag)
    }
end

* Function to load Mata data into Stata
cap program drop load_data
program define load_data,
    syntax, nloc(int) [cl1path(str)] [cl2path(str)] [distpath(str)] [savedyad] [filesuffix(str)] [scpcwfin(str)] [scpcy1path(str)] [scpcy2path(str)]

    clear

    mata: st_local("obsN",strofreal(rows(CovEpsVec),"%50.0f"))
		
    gen id1=.
    gen id2=.
    gen corr=.

    qui set obs `obsN'

    mata: st_store(.,.,(id_widerowvec,id_widecolvec,CovEpsVec))

    qui compress
    tempfile corr
    qui save `corr'

    clear

    mata: st_local("obsN",strofreal(rows(ResXtildeVec),"%50.0f"))

    gen id1=.
    gen id2=.
    gen xxresxx=.

    qui set obs `obsN'

    mata: st_store(.,.,(ResXtildeVec))

    qui compress
    qui merge 1:1 id1 id2 using `corr', assert(2 3) nogen

    if _N!=(`nloc'^2+`nloc')/2 {
        di as error "Number of location pairs is incorrect"
        exit
    }

    if "`cl1path'"!="" {
        qui merge m:1 id1 using `cl1path', assert(3) nogen
        qui merge m:1 id2 using `cl2path', assert(3) nogen
        gen byte same_cl = 0
        qui ds __tmo_cl1_*
        foreach var1 in `r(varlist)' {
            local var2 = subinstr("`var1'","_cl1_","_cl2_",1)
            qui replace same_cl = 1 if `var1'==`var2'
        }
    }

    if "`distpath'"!="" {
        qui merge 1:1 id1 id2 using `distpath', assert(3) nogen
    }

    if "`scpcwfin'"!="" {
        qui merge 1:1 id1 id2 using `scpcwfin', assert(1 3) nogen
        qui merge m:1 id1 using `scpcy1path', assert(1 3) nogen
        qui merge m:1 id2 using `scpcy2path', assert(1 3) nogen
        qui gen ryyr = (Wfin*y_scpc1*y_scpc2)*(1 + (id1!=id2))
    }
    
    * Fisher-transform correlation
    qui count if abs(corr)>1 & !missing(corr)
    if `r(N)'>0 {
        di as error "Correlations <-1 or >1 exist"
        exit
    }

    qui gen corr_fisher = 0.5 * ln((1+corr)/(1-corr))
    qui compress

    if "`savedyad'"!="" {
        save "`filesuffix'_dyad.dta", replace
    }
end

* Function to compute degrees of freedom and optimal threshold
cap program drop dfqt
program define dfqt
    syntax, [plotq] [plothist] nbins(int) nloc(int) [filesuffix(str)]

    cap drop offdiag
    qui gen byte offdiag = !missing(corr) & id1!=id2 & abs(corr)!=1

    qui gstats sum corr_fisher if offdiag
    scalar sd = (r(p75)-r(p25))/(invnormal(0.75)-invnormal(0.25))
	scalar df = 1/(sd^2)
	
    cap drop corr_fisher_abs 
    qui gen corr_fisher_abs = abs(corr_fisher)
    qui replace corr_fisher_abs = . if !offdiag // will be at end
    qui hashsort corr_fisher_abs

    cap drop cdf_emp_abs_1min cdf_iqr_abs_1min q_iqr_abs
    qui count if offdiag
    qui gen cdf_emp_abs_1min = 1 - (_n/`r(N)') if offdiag
	qui gen cdf_iqr_abs_1min = 1 - 2*(normal(corr_fisher_abs/sd)-0.5) if offdiag
    qui gen q_iqr_abs = cdf_emp_abs_1min - 2*cdf_iqr_abs_1min if offdiag

    qui sum q_iqr_abs if offdiag
    qui sum corr_fisher_abs if abs(q_iqr_abs-`r(max)')<=1e-10 & offdiag
    scalar fthres = `r(min)'
	scalar thres = tanh(fthres)

    cap drop pdf_iqr
    qui gen pdf_iqr = normalden(corr_fisher,sd) if offdiag
    
    if "`plotq'"!="" {        
        qui sum corr_fisher_abs if q_iqr_abs>=-0.002 & offdiag
        local xstart=floor(`r(min)'*10)/10
        local fthres=fthres

        twoway ///
                (line q_iqr_abs corr_fisher_abs if q_iqr_abs>=-0.002 & corr_fisher_abs<=1, lcolor(blue) lwidth(medthick) sort(corr_fisher_abs)), ///
                graphregion(color(white)) ///
                yline(0, lcolor(gray)) ///
                xtitle("Threshold {it:{&delta}}") ytitle("") ///
                ylab(, angle(horizontal) format("%04.3f")) ///
                xlab(`xstart'(0.1)1, grid gmin gmax format("%02.1f")) ///
                xline(`fthres', lcolor(red)) ///
                xsize(16) ysize(9)
        graph export "`filesuffix'_qt.pdf", as(pdf) replace
    }

    if "`plothist'"!="" {
        qui sum corr_fisher if offdiag
        local binwidth = (`r(max)'-`r(min)')/`nbins'
        qui gen bin_corr = floor((corr_fisher-`r(min)')/`binwidth') if offdiag
        qui gen bin_cent = bin_corr*`binwidth' + `binwidth'/2 + `r(min)'
        qui sum bin_corr if offdiag
        assert abs(`r(max)'-`nbins')<=1
        qui hashsort bin_corr corr_fisher
        qui by bin_corr: gen binN=_N if offdiag
        qui gen bin_dens = binN/(_N*`binwidth') if offdiag
        qui gegen byte bin_tag = tag(bin_corr) if offdiag
        
        local thres_str = string(thres,"%03.2f")
        local fthres_str = string(fthres,"%03.2f")
        local df_str = string(df,"%05.2f")
        local fthres = fthres
        
        twoway 	(bar bin_dens bin_cent if bin_tag==1 & corr_fisher>=-1 & corr_fisher<=1, base(0) barwidth(`binwidth') color(midgreen%30)) ///
                    (line pdf_iqr corr_fisher if corr_fisher>=-1 & corr_fisher<=1, sort(corr_fisher) lcolor(blue%90)) ///
                    , graphregion(color(white)) ///
                    xtitle("Fisher transformed correlation") ///
                    ytitle("Density") ///
                    ylab(, format(%02.1f) angle(horizontal)) ///
                    xlab(-1(0.2)1, format(%02.1f)) ///
                    xline(0, lcolor(gray) lpattern(longdash)) ///
                    xline(`fthres', lcolor(red)) ///
                    xline(-`fthres', lcolor(red)) ///
                    legend(order(2 "IQR df=`df_str', {it:{&delta}}{sup:*}=`thres_str' (`fthres_str' Fisher transformed)") ///
                            pos(6) nobox region(color(none))) ///
                    xsize(16) ysize(9)
        graph export "`filesuffix'_hist.png", as(png) replace
    }
end

* Function to estimate TMO SE
cap program drop est_tmo_se
program define est_tmo_se,
    syntax, [dist_cutoff(real 0)] [custom_thres(real -9)] [thresholdoff] [`scpc_uncond']

    if `custom_thres'!=-9 {
        scalar thres = `custom_thres'
        scalar fthres = 0.5 * ln((1+thres)/(1-thres))
        scalar df = .
        scalar sd = .
    }

    if "`thresholdoff'"!="" {
        scalar thres = .
        scalar fthres = .
        scalar df = .
        scalar sd = .
    }

    cap confirm var same_cl
    if !_rc {
        local orig_cond (same_cl)
    }
    else {
        local orig_cond (id1==id2)
    }

    cap confirm var same_cl
    if !_rc {
        if `dist_cutoff'>0 {
            local keep_cond (same_cl | dist<=`dist_cutoff')
        }
        else {
            local keep_cond (same_cl)
        }
    }
    else {
        if `dist_cutoff'>0 {
            local keep_cond (id1==id2 | dist<=`dist_cutoff')
        }
        else {
            local keep_cond (id1==id2)
        }
    }

    * For SCPC option
    cap confirm var Wfin
    if !_rc {
        local scpc = 1
        * Check SCPC SE (might not equal if uncond option or if missing some locations due to missing coordinates)
        qui count if missing(ryyr)
        if r(N)==0 & "`scpc_uncond'"=="" {
            qui sum ryyr
            if abs(scpc_se - sqrt(r(sum)))>1e-5 {
                di as error "SCPC SE does not match"
                exit
            }
        }
    }
    else {
        local scpc = 0
    }

    * Back out finite sample degrees of freedom adjustment from original SE
    qui sum xxresxx if `orig_cond'
    scalar dof_adj = se^2/r(sum)

    * TMO SE
    qui sum xxresxx if ((abs(corr)>=thres) & !missing(corr)) | `keep_cond'
    local xxresxx_sum = r(sum)
    if `scpc'==1 {
        qui sum ryyr if (abs(corr)<thres | missing(corr)) & !`keep_cond'
        local ryyr_sum = r(sum)
    }
    else {
        local ryyr_sum = 0
    }
	scalar tmo_se = sqrt((`xxresxx_sum'+`ryyr_sum')*dof_adj)

    * Store no. off-diag
    qui count if id1!=id2
    scalar offdN = r(N)
    qui count if !`orig_cond'
    scalar offdNnocl = r(N)

    * Proportion off-diag included in SE calculation (including both clustering, Conley, and TMO)
	qui count if (((abs(corr)>=thres) & !missing(corr)) | `keep_cond') & id1!=id2
	scalar offdP = r(N)/offdN

    * Proportion off-diag over threshold outside clusters or Conley
	qui count if (abs(corr)>=thres) & !missing(corr) & !`keep_cond'
	scalar offdPnocl = r(N)/offdNnocl
end

* Function to plot TMO SE over thresholds
cap program drop tmo_over_thres
program define tmo_over_thres,
    syntax, filesuffix(str) [saveplotseest] [NOIsily] [dist_cutoff(real 0)]

    cap confirm var same_cl
    if !_rc {
        if `dist_cutoff'>0 {
            local keep_cond (same_cl | dist<=`dist_cutoff')
            local subtitle2 `" "(excluding within-cluster and -`dist_cutoff' correlations)" "'
        }
        else {
            local keep_cond (same_cl)
            local subtitle2 `" "(excluding within-cluster correlations)" "'
        }
    }
    else {
        if `dist_cutoff'>0 {
            local keep_cond (id1==id2 | dist<=`dist_cutoff')
            local subtitle2 `" "(excluding within-`dist_cutoff' correlations)" "'
        }
        else {
            local keep_cond (id1==id2)
            local subtitle2
        }
    }

    cap confirm var Wfin
    if !_rc {
        local scpc = 1
    }
    else {
        local scpc = 0 
    }

    * Initialize matrix to store results
    mat tmo_over_thres = J(101,4,.)
	mat colnames tmo_over_thres = delta tmo_se offdP offdPnocl 

    local row=1
    forv thr=0(0.01)1.01 {
        if "`noisily'"!="" {
            if mod(`row',20)==1 {
                local thr_str = string(`thr',"%03.2f")
                di "Calculating TMO SE over threshold at `thr_str'"
            }
        }

        mat tmo_over_thres[`row',1]=`thr'
		qui sum xxresxx if ((abs(corr)>=`thr') & !missing(corr)) | `keep_cond'
        local xxresxx_sum = r(sum)
        if `scpc'==1 {
            qui sum ryyr if (abs(corr)<thres | missing(corr)) & !`keep_cond'
            local ryyr_sum = r(sum)
        }
        else {
            local ryyr_sum = 0
        }    
		mat tmo_over_thres[`row',2]=sqrt((`xxresxx_sum'+`ryyr_sum')*dof_adj)
		qui count if (abs(corr)>=`thr') & !missing(corr) & !(id1==id2)
		mat tmo_over_thres[`row',3]=r(N)/offdN
		qui count if (abs(corr)>=`thr') & !missing(corr) & !(`keep_cond')
		mat tmo_over_thres[`row',4]=r(N)/offdNnocl

        local ++row
    }

    * Plot TMO SE over threshold
    clear
    qui svmat2 tmo_over_thres, names(col)
    qui compress

    if "`saveplotseest'"!="" {
        save "`filesuffix'_tmo_se_over_thres.dta", replace
    }

    qui gen tmo_orig_se_ratio = tmo_se/se
    qui replace tmo_orig_se_ratio = 0 if missing(tmo_orig_se_ratio)

    local thres = thres
    twoway 	(line tmo_orig_se_ratio delta, lcolor(black) lwidth(medthick)) ///
			, ///
			graphregion(color(white)) ///
			xlabel(0(0.1)1.0, format("%02.1f") grid gmin gmax) ///
			ylab(, angle(horizontal) gmin gmax) ///
			xtitle("Threshold {it:{&delta}}") ///
			subtitle("Ratio of TMO standard error to original", pos(11)) ytitle("") ///
			legend(off) ///
			yline(1, lcolor(gray%50) lwidth(thick)) ///
			xline(`thres', lcolor(red) lwidth(medthick)) ///
			xsize(16) ysize(9)
	graph export "`filesuffix'_se_ratio_over_thres.pdf", as(pdf) replace

    twoway 	(line offdPnocl delta if offdPnocl<=0.1, lcolor(black) lwidth(medthick)) ///
			, ///
			graphregion(color(white)) ///
			ylab(0(0.01)0.1, format("%03.2f") angle(horizontal) gmin gmax) ///
			xlabel(0(0.1)1.0, format("%02.1f") grid gmin gmax) ///
			xtitle("Threshold {it:{&delta}}") ///
			subtitle("Proportion of correlations {&ge} threshold" `subtitle2', pos(11)) ytitle("") ///
			legend(off) ///
			xline(`thres', lcolor(red) lwidth(medthick)) ///
			xsize(16) ysize(9)
	graph export "`filesuffix'_prop_above_thres.pdf", as(pdf) replace
end

* Function to save TMO results to dta file
cap program drop tmo_save
program define tmo_save,
    syntax, FILEsuffix(str)

    preserve
        clear

        qui set obs 1

        qui gen beta = beta
        qui gen orig_se = se
        qui gen tmo_se = tmo_se
        qui gen lb = lb
        qui gen ub = ub
        qui gen threshold = thres
        qui gen pct_ge_thres = offdP*100
        qui gen pct_ge_thres_nocl = offdPnocl*100
        qui gen T = T
        qui gen N_loc = N
        qui gen N = N_obs
        qui gen N_clust = N_clust
        qui gen N_outcomes = D
        qui gen dof = df
        qui gen finite_sample_dof = dof_adj
        qui gen df_r = df_r
        qui gen scpc_cv = ${scpc_cv}

        qui compress
        save "`filesuffix'_est.dta", replace
    restore

end
