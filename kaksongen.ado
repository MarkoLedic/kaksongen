* version 1.0.0, July 2022
* Marko Ledić, Ivica Rubil, Ivica Urban

cap prog drop kaksongen
prog kaksongen, rclass sortpreserve


	version 15.1

	 
	syntax varlist(min=2 max=2 numeric) [aweight] [if] [in],       ///
		   ALPHA(numlist min=1 max=101 ascending >=0 <=1)      ///
		   [                                                   ///
		   RHO(real 2)                                         /// 
		   ASVARiables                                         ///
		   NODIVbytax                                          ///
		   GRaph                                               ///
		   *                                                   ///
		   ]
	

	marksample touse

	
	// count how many alphas are specified
	loc nofalphas : list sizeof local(alpha)
	

	// create: (a) locals containing: 1, 2, ..., nofalphas
	//         (b) locals containing the values of the specified alphas
	tokenize `alpha'
	forval a = 1/`nofalphas' {
		loc a`a' `a'               // (a)
		loc alpha`a' ``a''         // (b)
	}

	
	// parse varlist
	tokenize `varlist'
	tempvar x t
	gen `x' = `1'
	gen `t' = `2'

	
	// option alpha() is required
	if "`alpha'" == "" {
		di as err "ERROR: option alpha() is required"
		di
		error 198
	}	
	
	
	// rho must be greater than 1
	if `rho' <= 1 {
		di as err "ERROR: rho must be greater than 1"
		di
		error 198
	}
	
	
	// tax must not be larger than pre_tax income
	tempname mux mut
	
	qui sum `x' [`weight' `exp'] if `touse'
	sca `mux' = r(mean)
	
	qui sum `t' [`weight' `exp'] if `touse'
	sca `mut' = r(mean)
	
	if `mux' < `mut' {
		di as err "ERROR: mean pre-tax income is smaller than mean tax"
		exit
	}
			
	
	// post-tax income
	tempvar y
	gen double `y' = `x' - `t'
	qui sum `y' [`weight' `exp'] if `touse'
	tempname muy
	sca `muy' = r(mean)
	

	// compute decomposition
	tempname Gx Gy Dy Dt	
	qui {
	sgini `x' [`weight' `exp'] if `touse', param(`rho') 
	sca `Gx' = r(coeff)
	sgini `y' [`weight' `exp'] if `touse', param(`rho') 
	sca `Gy' = r(coeff)
	sgini `y' [`weight' `exp'] if `touse', param(`rho') sort(`x')
	sca `Dy' = r(coeff)
	sgini `t' [`weight' `exp'] if `touse', param(`rho') sort(`x')
	sca `Dt' = r(coeff)
	}
	
	tempname Wx Wy
	sca `Wx' = `mux'*(1-`Gx')
	sca `Wy' = `muy'*(1-`Gy')
	
	tempname DeltaW H
	if "`nodivbytax'" == "" {
		sca `DeltaW'  = (1/`mut')*(`Wy'-`Wx')
		sca `H' = (`muy'/`mut')*(`Dy'-`Gy')
		forval a = 1/`nofalphas' {	
			tempname N_alpha`a' P_alpha`a' 
			sca `N_alpha`a'' = -(1 - `alpha`a'' * `Gx')  
			sca `P_alpha`a'' = `Dt' - `alpha`a'' * `Gx'   
		}
	}
	else {
		sca `DeltaW'  = (`Wy'-`Wx')
		sca `H' = `muy'*(`Dy'-`Gy')
		forval a = 1/`nofalphas' {	
			tempname N_alpha`a' P_alpha`a' 
			sca `N_alpha`a'' = -`mut'*(1 - `alpha`a'' * `Gx')  
			sca `P_alpha`a'' = `mut'*(`Dt' - `alpha`a'' * `Gx')   
		}	
	}

	forval a = 1/`nofalphas' {
		tempname delta_alpha`a' pi_alpha`a' eta_alpha`a'
		sca `delta_alpha`a''= `DeltaW'     / abs(`N_alpha`a'')
		sca `pi_alpha`a''   = `P_alpha`a'' / abs(`N_alpha`a'')
		sca `eta_alpha`a''  = `H'          / abs(`N_alpha`a'')
	}
	
	tempname alpha_zero
	sca `alpha_zero' = `Dt'/`Gx'
	
	tempname tau 
	sca `tau' = `mut'/`mux'
	
	forval a = 1/`nofalphas' {	
		tempname tau_star_alpha`a' tau_over_tau_star_alpha`a'
		sca `tau_star_alpha`a'' = -`delta_alpha`a''*`tau'
		sca `tau_over_tau_star_alpha`a'' = `tau'/`tau_star_alpha`a''
	}
	
	forval a = 1/`nofalphas' {
		mat results`a'   = `alpha`a''                                   \ ///
		                   `rho'                                        \ ///
						   `DeltaW'                     \ ///
						   `N_alpha`a''                 \ ///
						   `P_alpha`a''                 \ ///
						   `H'                          \ ///
						   `delta_alpha`a''             \ ///
						   `pi_alpha`a''                \ ///
						   `eta_alpha`a''               \ ///
						   `alpha_zero'                 \ ///
						   `tau'                        \ ///
						   `tau_star_alpha`a''          \ ///
						   `tau_over_tau_star_alpha`a'' \ ///
						   `Gx'                         \ ///
						   `Gy'                         \ ///
						   `Dy'                         \ ///
						   `Dt'                         \ ///
						   `mux'                        \ ///
						   `muy'                        \ ///
						   `mut'                      
						   
		mat rownames results`a' = alpha rho DeltaW N_alpha P_alpha H ///
		                          delta_alpha pi_alpha eta_alpha alpha_zero ///
								  tau tau_star_alpha tau_over_tau_star_alpha ///
								  Gx Gy Dy Dt mux muy mut 
		mat colnames results`a' = value
	}

	
	if `nofalphas' == 1 {
		mat results = results1
	}
	else {
		mat results = results1
		forval a = 2/`nofalphas' {
			mat results = results , results`a'
		}	
	}
	
	
	// required for option -variables- 
	mat resultsT  = results'
	
	
	// display results
	mat li results
	
	
	// return results as matrix
	ret mat results = results
		
	
	// return some macros
	ret loc rho = "`rho'"
	ret loc alpha = "`alpha'"
	ret loc cmd = "`0'"       // full command
	ret loc xvar = "`1'"
	ret loc tvar = "`2'"
	
	
	// optionally: results as variables
	if "`asvariables'" == "asvariables" {
		qui svmat resultsT, names(col)
	}
		
	// graph 
	if ("`graph'" == "graph" & "`asvariables'" != "asvariables") {
		di as err "ERROR: Option asvariables is required"
		di
		error 198
	}
	
	forval a = 1/`nofalphas' {
		loc oneover`a' = `a'/`nofalphas'
	}
	
	if ("`graph'" == "graph" & "`asvariables'" == "asvariables") {
			graph twoway con DeltaW N_alpha P_alpha H alpha,     ///
            xti("{&alpha}")                                                  ///
            leg(order(1 "{&Delta}W" 2 "N({&alpha})" 3 "P({&alpha})" 4 "H"))  ///
			`options'			
	}

end



