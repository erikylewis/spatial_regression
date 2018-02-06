/*
Eric Lewis
February 2018
Modified from code from Solomon Hsiang 

syntax: spatial_reg LHS Exog_RHS_vars (Endog_RHS_vars = Exog_Inst) [if] [in]
   lat(latitude variable) -- for calculate distance option
   		--> specify latitude in decimal degrees 
   lon(longitude variable) -- for calculate distance option
   		--> specify longitude in decimal degrees 
   indexvar -- for  where distance is already calculated
   		--> relies on a matrix that gives distances (called distmat)
   distcutoff -- at what distance are positive weights defined 
   incrementsize -- used for cases where indexvar option is used
   		--> affects how quickly distance extraction is done
   		--> too small --> number of loops is large, takes time
   		--> too large --> extracts too large of submatrices from distmat, 
   				which can be computationally taxing
   small -- adjusts variance-covariance matrix for small sample size
   		--> multiplies it by (n/(n-k))
   		--> similar to the "small" option for ivregress
   		--> should be used for OLS
   bartlett -- uses Bartlett weighting. Default is to use uniform weighting.
 
Allows for both OLS and 2SLS
	--> 2SLS includes the syntax (Endog_RHS_vars = Exog_Inst) 
	--> very similar to the difference in syntax for regress vs. ivregress




Code allows for two ways to calculate distances:
	use lat() and lon() options
	Or can use a pre-determined matrix 
		The predetermined matrix MUST be saved in Mata
		The predetermined matrix MUST be named distmat (in Mata)
		The variable indexvar indexes the position of an observation
			within the distmat matrix
		For example, if for a given observation, if indexvar=i, then submatrix:
			distmat[i,.] gives the vector of distances between that observation
			and all other locations observed
		This code allows for distance matrices that exceed Stata's maximum matrix size.
			This is because it relies on the MATA command fopen(), which can read and
			write large matrices without running them through Stata.


This regression specification does not include an automatic constant! If needed, a contant
     should be added in as one of the covariates.
			

*/ 

program spatial_reg, eclass byable(recall)
version 14

#delimit ;
syntax anything [if] [in] [, 
				lat(varname numeric)
				lon(varname numeric)
				INDEXvar(varname numeric)
				DISTcutoff(real 1) DIPROGress
				DISPlay star bartlett dropvar small 
				INCrementsize(int 300) ] 
				 ;
#delimit cr

/*--------PARSING COMMANDS AND SETUP-------*/

	// Finds the lhs var, exogenous variables, endogenous variables, and instruments
	gettoken y rest : anything
	// di "`y'"
	// di "`rest'"
	di as text "Left hand side variable: " as input "`y'"
	local rest_cmdline `rest'
	
	// get exogenous vars
	gettoken exog rest : rest, parse("(")
	// di "`exog'"
	// di "`rest'"
	gettoken open rest : rest, parse("(")
	// di "`rest'"
	
	// get endogenous vars
	gettoken endog rest : rest, parse("=")
	gettoken equal rest : rest, parse("=")
	// di "`endog'"
	// di "`rest'"
	
	// get instruments
	gettoken inst close : rest, parse(")")
	// di "`inst'"
	// di "`close'"
	
	

/*----------DETERMINES THE SAMPLE-----------------*/
	capture drop touse
	marksample touse				// indicator for inclusion in the sample
	markout `touse' `y' `exog' `endog' `inst'   
		// if variables in `y', `exog', `endog', or `inst' are missing, 
		// excludes those observations from the sample `touse'
	gen touse = `touse'



/*---Checks whether information is included to do distance calculations----*/
// If a distance matrix is specified, checks that it exists	
	#delimit ;
	if ("`lat'"=="" & "`indexvar'"=="") |
			(("`lat'"=="" & "`lon'"!="") | 
			 ("`lat'"!="" & "`lon'"=="")) { ;
			di as error "Must specify either (latitude AND longitude) OR indexvar" ;
		exit ;
	} ;
	#delimit cr
	if "`lat'"!="" & "`lon'"!="" local disttype = "lat_lon"
	if "`indexvar'"!="" local disttype = "distance_matrix"
	// Checks that the matrix exists if specified
	if "`disttype'" == "distance_matrix" {
		di as text "distance index variable: " as input "`indexvar'"
		mata: mata describe distmat
		capture mata: mata describe distmat
		if _rc!= 0 {
			di as error "There is no matrix loaded into Mata called distmat"
			exit
		}
	}
	if "`disttype'" == "lat_lon" {
		di as text "Uses latitude variable " as input "`lat'" as text "; longitude variable " as input "`lon'"
	}	

	// Determines whether OLS or 2SLS will be used:
	// Determines whether `endog' and `inst' are defined and whether parentheses were used
		// First determines whether "(", ")", and "=" are included:
		local has_leftp  = strpos("`anything'","(")>0 & strpos("`anything'","(")<.
		local has_rightp = strpos("`anything'",")")>0 & strpos("`anything'",")")<.
		local has_equal = strpos("`anything'","=")>0 & strpos("`anything'","=")<.
	
	if "`endog'"!="" & "`inst'"!="" & `has_leftp' & `has_rightp' & `has_equal'  {
		local reg_type = "2sls"
	}
	if "`endog'"=="" & "`inst'"=="" & ~`has_leftp' & ~`has_rightp' & ~`has_equal'  {
		local reg_type = "ols"
	}
	if ~inlist("`reg_type'","ols","2sls") {
		di as error "Invalid syntax: does not appear to be 2SLS or OLS syntax"
		exit
	}

/*----------DETERMINES IF indexvar VARIABLE IS WRONG-----------------*/
	if "`indexvar'"!="" {
		qui count if missing(`indexvar') & `touse'
		if r(N)>0 {
			di as error "indexvar variable has missing values"
			exit
		}
		mata: rowsof_distmat = rows(distmat)
		mata: st_matrix("rowsof", rowsof_distmat) 
		local rowsof_distmat = el(rowsof,1,1)
		
		qui count if (`indexvar'>`rowsof_distmat') & `touse'
		if r(N)>0 {
			di as error "indexvar variable exceeds the # of rows of distmat"
			exit
		}
		qui count if `indexvar'<1 & `touse'
		if r(N)>0 {
			di as error "indexvar variable has values less than 1"
			exit
		}
		qui count if `indexvar'!=ceil(`indexvar') & `touse'
		if r(N)>0 {
			di as error "indexvar variable has non-integer values"
			exit
		}
	}


/*--------CHECKS WHETHER exog variable has a constant as the last term-----*/
	local numexogvars :  list sizeof exog
	local lastexogvar : word `numexogvars' of `exog'
	qui sum `lastexogvar' if `touse'
	if r(sd) != 0 {
		di as error "Warning: " as text "Last variable in exogenous variable list is not a constant"
	}


//generating a function of the included obs
	quietly count if `touse'		
	scalar n = r(N)					// # obs
	scalar n_obs = r(N)
	local n_obs = r(N)


/*-------GETS BETA COEFFICIENTS AND RESIDUALS--------------*/
/*--------DOES OLS, STORE RESULTS-------*/
if "`reg_type'"=="ols" {
	quietly: reg `y' `exog' if `touse', hascons

		local r2_old = e(r2)
		local df_m_old = e(df_m)
		local rmse_old = e(rmse)
		local mss_old = e(mss)
		local rss_old = e(rss)
		local r2_a_old = e(r2_a)
		local n_old = e(n)
		local df_r_old = e(df_r)
		local rank_old = e(rank)

		tempvar e_resid
		qui predict `e_resid', xb
		qui replace `e_resid' = `y' - `e_resid'
		local e_resid_name "`e_resid'"
		
		// Gets rid of colinear variables--those that have o. at the beginning of the name
		matrix beta_spatial = e(b)
		local depvars_pre_decolin : colfullnames e(b)
		local depvars_post_decolin
		local i_colinearity = 1
		local depvars_notcolin_columns
		foreach word in `depvars_pre_decolin' {
			if substr("`word'",1,2)!="o." {
				local depvars_post_decolin `depvars_post_decolin' `word'
				local depvars_notcolin_columns `depvars_notcolin_columns' `i_colinearity'
			}
			local i_colinearity = `i_colinearity' + 1
		}		
		matselrc beta_spatial beta_spatial_decolin, col(`depvars_notcolin_columns')
		matrix colnames beta_spatial_decolin = `depvars_post_decolin'
		matrix beta_spatial = beta_spatial_decolin
		matrix drop beta_spatial_decolin
		
		/*--------New Locals with names of variables to plug into Mata----*/
		local X "`depvars_post_decolin'"
		local X_var_names  "`depvars_post_decolin'"
		local Y "`y'"  // This one isn't really needed
}	

/*--------DOES 2SLS, STORE RESULTS-------*/
if "`reg_type'"=="2sls" {
// di "C"
	quietly: ivregress 2sls `y' `exog' (`endog' = `inst') if `touse', hascons
	
		local r2_old = e(r2)
		local df_m_old = e(df_m)
		local rmse_old = e(rmse)
		local mss_old = e(mss)
		local rss_old = e(rss)
		local r2_a_old = e(r2_a)
		local n_old = e(n)
		local rank_old = e(rank)

		if "`small'"=="small" {
			local df_r_old = `n_old' - `rank_old'
		}

		tempvar e_resid
		qui predict `e_resid', xb
		qui replace `e_resid' = `y' - `e_resid'
		local e_resid_name "`e_resid'"

		matrix beta_spatial = e(b)
		// for 2SLS, Stata does not automatically put in the o._cons term
		// when hascons is specified. Unlike above.

		*####* [ Fixes this to get rid of colinear estimates...
		local depvars_pre_decolin : colfullnames beta_spatial
		local depvars_post_decolin
		local i_colinearity = 1
		local depvars_notcolin_columns
		foreach word in `depvars_pre_decolin' {
			if substr("`word'",1,2)!="o." {
				local depvars_post_decolin `depvars_post_decolin' `word'
				local depvars_notcolin_columns `depvars_notcolin_columns' `i_colinearity'
			}
			local i_colinearity = `i_colinearity' + 1
		}		
		matselrc beta_spatial beta_spatial_decolin, col(`depvars_notcolin_columns')
		di "`depvars_post_decolin'"
		di "`dpvars_notcolin_columns'"
		matrix colnames beta_spatial_decolin = `depvars_post_decolin'
		matrix beta_spatial = beta_spatial_decolin
		matrix drop beta_spatial_decolin

		* Creates a new list of endogenous variables:
		local endog_notcolinear : list depvars_post_decolin & endog 
		local exog_notcolinear : list depvars_post_decolin - endog_notcolinear

		

	/*--------Creates Predicted Endogenous variables----*/
	local X
	local endog_hat
	local X_var_names
	foreach i in `endog_notcolinear' {
		qui reg `i' `exog' `inst' if `touse', nocons
		predict `i'___HAT, xb
		local X "`X' `i'___HAT"
		local endog_hat "`endog_hat' `i'___HAT"
		local X_var_names "`X_var_names' `i'"
	}
	/*--------New Locals with names of variables to plug into Mata----*/
	local X "`X' `exog_notcolinear'"
	local X_var_names  "`X_var_names' `exog_notcolinear'"
	local Y "`y'" // This one isn't really needed

		*####* ]

	
}

// di "E"


	mata{
	// "F"

	/*--------IMPORT ALL VALUES INTO MATA-------*/
	disttype = st_local("disttype")
	// disttype
		
	// Don't really need Y:
	// Y_var = st_local("Y") //importing variable assignments to mata (e.g., the name)
	// st_view(Y=.,.,tokens(Y_var),"touse")
	// Y = st_data(.,tokens(Y_var),"touse")
	
	X_var = st_local("X") // Imports the names of all the predicted endogenous HAT variables
	st_view(X=.,.,tokens(X_var),"touse")
	// X = st_data(.,tokens(X_var),"touse")
	
	// e_resid_name:
	e_resid_var = st_local("e_resid_name")
	st_view(e=.,.,tokens(e_resid_var),"touse")
	
	//#st_view(e=.,.,("___e"),"touse")
	// e = st_data(.,("___e"),"touse")
	
	if (disttype=="lat_lon") {
		lat_var = st_local("lat")
		st_view(lat=.,.,tokens(lat_var),"touse")
		// lat = st_data(.,tokens(lat_var),"touse")
	
		lon_var = st_local("lon")
		st_view(lon=.,.,tokens(lon_var),"touse")
		// lon = st_data(.,tokens(lon_var),"touse")
	}
	if (disttype=="distance_matrix") {
		index_var = st_local("indexvar")
		st_view(dmat_index=.,.,tokens(index_var),"touse")
		// dmat_index = st_data(.,tokens(index_var),"touse")
	
	}
	
	// "G"
	
	k_variables = cols(X)				//importing other parameters
	
	n = st_numscalar("n")
	
	// lag_var = st_local("lagcutoff") // This is from S Hsiang's code
	// lag_cutoff = strtoreal(lag_var) // This is from S Hsiang's code
	dist_var = st_local("distcutoff")
	dist_cutoff = strtoreal(dist_var)
	
	
	
	/*-----------COMPUTE INNER TERM (X-prime e e-prime X) OF VARIANCE MATRIX--------------*/
	
	XeeX = J(k_variables,k_variables,0)
	// "H"
	// disttype 
	// should display what the disttype is--lat_lon or the other one
	
	// CASE WHERE CALCULATES DISTANCE MANUALLY
	if (disttype=="lat_lon") {
	//"Eye"
		// Loops over the observations
		for (i = 1; i <= n; i++) {
			if (i/1000==ceil(i/1000) & "`diprogress'"=="diprogress") {
				"Progress: on observation "+strofreal(i)
			}
	
			latdif_i = (lat[i,1] :- lat)*pi()/180
			londif_i = (lon[i,1] :- lon)*pi()/180
			// This uses the Haversine formula: http://en.wikipedia.org/wiki/Haversine_formula
			distance_i = 2*6371*asin((sin(latdif_i/2):^2+cos(lat[i,1]*pi()/180)*cos(lat*pi()/180):*(sin(londif_i/2):^2)):^.5)
			window_i = distance_i :<= dist_cutoff
			
			// Bartlett adjustment:
			if ("`bartlett'"=="bartlett") {
				weight_i = 1 :- distance_i:/dist_cutoff
				window_i = window_i:*weight_i
			} 
			XeeXh = ((X[i,.]'J(1,n,1)*e[i,1]):*(J(k_variables,1,1)*e':*window_i'))*X
			XeeX = XeeX + XeeXh
		}
	}
	
	
	if (disttype=="distance_matrix" ) {
	// "J"
		// loops over all observations
		// assumes that the matrix containing distance measures is named
		//      distmat
		incrementsize = strtoreal(st_local("incrementsize"))
		// Not sure if this works when increment size is less than n...
		// But this is unlikely to happen in normal life
	
		loop_ceiling = ceil(n/incrementsize)
	
	
	
			if ("`diprogress'"=="diprogress") {
				"Running through "+strofreal(loop_ceiling)+"*"+strofreal(loop_ceiling+1)+"/2 = "+strofreal(loop_ceiling*(loop_ceiling+1)/2)+" iterations"
			}
	
		for (i = 1; i <= loop_ceiling; i++){
			if ("`diprogress'"=="diprogress") {
				if (i==1) {
					numcompleted=0
					totcompleted=0
				}
				if (((i-1)/100)==ceil((i-1)/100)) {
					stata(`"di _newline as text c(current_time) _continue"')
				}
				stata(`"di as text "." _continue"')
				if (i>1) {
					numcompleted = 1+loop_ceiling-(i-1)
					totcompleted = totcompleted+numcompleted
					
				}
				if (i==loop_ceiling) {
					stata(`"di "" "')
				}
			}
	
			lower_i = incrementsize*i-incrementsize+1
			upper_i = min((incrementsize*i, n))
	
			for (j = i; j <= loop_ceiling; j++){
	
				lower_j = incrementsize*j-incrementsize+1
				upper_j = min((incrementsize*j, n))
				distance_ij = distmat[dmat_index[lower_i..upper_i],dmat_index[lower_j..upper_j]] 
				window_ij = distance_ij :<= dist_cutoff
	
				if (sum(window_ij)>0){
					// Bartlett adjustment:
					if (("`bartlett'"=="bartlett") & dist_cutoff>0){
						window_ij = window_ij:*(1 :- distance_ij:/dist_cutoff)
					} 
					if (j==i){
						e_inner_ij = (e[lower_i..upper_i,.]*(e[lower_j..upper_j,.])'):*window_ij
						XeeXh = (X[lower_i..upper_i,.])'*e_inner_ij*(X[lower_j..upper_j,.])
	
						if (j==i) XeeX = XeeX + XeeXh
						if (j>i) XeeX = XeeX + XeeXh + XeeXh'
					}
					if (j>i) {
						e_i = e[lower_i..upper_i,.]
						e_j = e[lower_j..upper_j,.]
						X_i = X[lower_i..upper_i,.]
						X_j = X[lower_j..upper_j,.]
						
						where_i = select(1..rows(window_ij), (rowsum(window_ij):>0)')'
						where_j = select(1..cols(window_ij), (colsum(window_ij):>0))
						
						e_i_sub = e_i[where_i]
						e_j_sub = e_j[where_j]
						X_i_sub = X_i[where_i,.]
						X_j_sub = X_j[where_j,.]
						
						window_ij_sub = window_ij[where_i, where_j]
						
						e_inner_ij = (e_i_sub*e_j_sub'):*window_ij_sub
						XeeXh = X_i_sub'*e_inner_ij*X_j_sub

						if (j==i) XeeX = XeeX + XeeXh
						if (j>i) XeeX = XeeX + XeeXh + XeeXh'

					}
				}
			} 
			// end of j loop
		} 
		// end of i loop
	} 
	// end of disttype=="distance_matrix" 
	
	/*-----------------Now calculates the rest of V---------------------*/
	
	// "K"
	
	XeeX = XeeX/n
	// should be symmetric already, but creates a matrix that is definitly symmetric
	XeeX = (XeeX + XeeX')/2
	
	invXX = luinv(X'*X/n)
	V = invXX * XeeX * invXX / n
	V = (V + V')/(2)
	
	st_matrix("V_spatial", V) 
	// Writes the matrix V_spatial to be available in Stata (not just Mata, where it has been known as V)
	
	} // ends mata
* ENDS MATA



if ("`small'"=="small" |"`reg_type'"=="ols" ) {
	matrix V_spatial = V_spatial*`n_obs'/(`n_obs'-rowsof(V_spatial))
}


// the row and column names of the new VCE must match the vector b

// matrix list V_spatial

scalar k = rowsof(V_spatial)
//local n_minus_k = n_old-k
//"got to N-K"
//local n_obs = n
//"got to N obs"
local n_minus_k = `n_obs'-rowsof(V_spatial)

// Temporary debugging code: displays row and column names of beta_spatial and V_spatial
local b_names : colfullnames beta_spatial
matrix colnames V_spatial = `b_names'
matrix rownames V_spatial = `b_names'


if "`small'"=="" {
	ereturn post beta_spatial V_spatial, esample(`touse') obs(`n_obs')
}
if "`small'"=="small" {
	ereturn post beta_spatial V_spatial, esample(`touse') dof(`n_minus_k') obs(`n_obs')
}

ereturn local cmd "spatial_reg"
ereturn local cmdline "spatial_reg `0'"


if "`reg_type'"=="2sls" {
	ereturn scalar rss = `rss_old'
	ereturn scalar N = `n_obs'
	ereturn scalar df_m = `df_m_old'
	ereturn scalar rmse = `rmse_old'
	ereturn scalar mss = `mss_old'
	ereturn scalar r2 = `r2_old'
	ereturn scalar r2_a = `r2_a_old'
	ereturn scalar r2_a = `r2_a_old'
	ereturn scalar rank = `rank_old'
	if "`small'"=="small" {
		ereturn scalar df_r = `df_r_old'
	}
}
if "`reg_type'"=="ols" {
	ereturn scalar rss = `rss_old'
	ereturn scalar N = `n_obs'
	ereturn scalar df_m = `df_m_old'
	ereturn scalar rmse = `rmse_old'
	ereturn scalar mss = `mss_old'
	ereturn scalar r2 = `r2_old'
	ereturn scalar r2_a = `r2_a_old'
	ereturn scalar rank = `rank_old'
}


if "`reg_type'"=="2sls" {
	disp as txt "2SLS regression: SE corrected for spatial dependence"
	disp as txt "     Dependent variable: `Y'"
	disp as txt "     Endogenous variables: `endog'"
	disp as txt "     Exogenous variables: `exog'"
	disp as txt "     Instrumental variables: `inst'"
	disp as txt "Spatial correlation kernal cutoff: `distcutoff' KM"
	if "`bartlett'" == "bartlett" {
		disp as txt "     (Note: Linear Bartlett window used for spatial kernal)"
	}
}
if "`reg_type'"=="ols" {
	disp as txt "OLS regression: SE corrected for spatial dependence"
	disp as txt "     Dependent variable: `Y'"
	disp as txt "     Independent variables: `exog'"
	disp as txt "Spatial correlation kernal cutoff: `distcutoff' KM"
	if "`bartlett'" == "bartlett" {
		disp as txt "     (Note: Linear Bartlett window used for spatial kernal)"
	}
}

ereturn display

/*-----------CLEANS THINGS UP---------------*/
#delimit ;
foreach item in V X X_var XeeX XeeXh dist_cutoff dist_var distance_i distance_ij 
	disttype dmat_index e e_resid_var e_inner_ij example_distmat i incrementsize 
	index_var invXX j k_variables lat lat_var latdif_i lon lon_var londif_i lower_i 
	lower_j n upper_i upper_j weight_i weight_ij window_i window_ij loop_ceiling
	numcompleted remaining totcompleted X_i X_i_sub X_j X_j_sub e_i e_i_sub
	e_j e_j_sub rowsof_distmat where_i where_j window_ij_sub { ;
	capture mata: mata drop `item' ;
} ; // crucially, keeps distmat so that it can be used in another regression ;
#delimit cr

//# cap drop ___e     // need to re-write code so that this is a temporary variable
cap drop `endog_hat' // need to re-write code so that these are temporary variables
cap drop touse
end // temporary end
exit // temporary exit




ereturn display // standard Stata regression table format

//------------------------------------------------------------------
// cleaning up Mata environment

capture mata mata drop V invXX  XeeX XeeXh XeeX_spatial_HAC window_t window_i weight t i ti pi X1 Y1 e1 time1 n1 lat lon lat1 lon1 lat_scale lon_scale rows_ti rows_pi timeUnique panelUnique Ntime Npanel X X_var XeeX_spatial Y Y_var b dist_cutoff dist_var distance_i k_variables lag_cutoff lag_var lat_var lon_var n panel panel_var time time_var weight_i


// drops created X_hat variables
qui drop `endog_hat'
qui drop touse
qui drop _est_TSLS



end


