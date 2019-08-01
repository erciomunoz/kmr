*! version 1.1 01August2019 - speeded up with some small changes (~8 times faster)
* version 1.0 19January2019

* This command implements Korinek, Mistiaen and Ravallion. Journal of Econometrics Vol. 135, Issue 1, Jan 2007

cap program drop kmr
program define kmr, eclass sortpreserve
	version 14.2
	
	syntax [varlist(default=none)] [in] [if], ///
			groups(string) ///
			Interviews(string) ///
			Nonresponse(string) [ ///
			NOCONstant		 ///
			SWeights(string) ///
			GENerate(string) ///
			graph(string)	///
			technique(string) ///
			delta(real 0.1) ///
			start(string) ///
			difficult	  ///
			maxiter(real 100) ///
			]
			
	marksample touse
			
	tempname betas V b n ngroups fvalue sigma aic schwarz new_group generate2
		
	if "`varlist'" == "" & "`noconstant'" != "" {
	di as error `"You need at least one explanatory variable with the option noconstant."'
	exit
	}
		
	if "`technique'" == "" {
	local technique = "nr"
	}
	
	if "`technique'" == "nm" & "`delta'" == "" {
	di as error `"You need to provide a value for delta to use with nm technique."'
	exit
	}
	
	if "`technique'" != "" & "`technique'" !="nr" & "`technique'" !="dfp" & "`technique'" !="bfgs" & "`technique'" !="bhhh" & "`technique'" !="nm" {
	di as error `"Invalid technique."'
	exit
	}
	
	qui: sum `varlist' if `touse'
	scalar `n' = r(N)
	
	// Create new group index from 1 to max(groups)
    qui: egen `new_group' = group(`groups') if `touse'
	sort `new_group'
	qui: ta `new_group'
	scalar `ngroups' = r(r)
	
	if "`generate'" == "" {
	local `generate2' = "P"
	} 
	else {
	local `generate2' = "`generate'"
	}
		
	mata: mywork("`varlist'","`new_group'","`interviews'","`nonresponse'","`b'","`touse'","`sweights'","`generate2'","`fvalue'","`V'","`sigma'","`aic'","`schwarz'","`noconstant'",`delta',"`start'","`difficult'",`maxiter')
	
	if "`generate'" != "" {
	gen `generate' = `generate2'
	gen `generate'_upper = `generate2'_upper
	gen `generate'_lower = `generate2'_lower
	}
	
	if "`noconstant'" == "" {
    local cnames "`varlist' _cons"
	matrix colnames `b' = `cnames'
	matrix colnames `V' = `cnames'
	matrix rownames `V' = `cnames'
    }
	else {
	local cnames "`varlist'"
    matrix colnames `b' = `cnames'
	matrix colnames `V' = `cnames'
	matrix rownames `V' = `cnames'	
	}
		
********************************************************************************
**** Show line graph ****
********************************************************************************

if "`graph'"!="" {

	sort `graph'
	line `generate2' `generate2'_upper `generate2'_lower `graph' if `touse', ytitle("Probability of response") xtitle(`graph') legend(off) lpattern("l" "." ".")

}
	
********************************************************************************
**** Return estimation and saved results ****
********************************************************************************

	di

	dis " "
	dis in green "Compliance function" ///
		_column (50) "Number of obs" _column(69) "=" _column(71) %8.0f in yellow `n'
	dis in green _column(50) "AIC" _column(69) "=" _column(71) %8.2f in yellow `aic'
	dis in green "Number of groups" _column(20) "=" _column(21) %5.0f in yellow `=`ngroups'' ///
      	in green _column(50) "Schwarz" _column(69) "=" _column(71)  %8.4f in yellow `schwarz'
		
	ereturn clear
	ereturn post `b' `V', esample(`touse') depname(`Probability of response')
	ereturn local technique `technique'
	ereturn local cmd     "kmr"
	ereturn local title       = "Compliance function estimate using group's response rates"
	ereturn local cmdline   = itrim("kmr `varlist', groups(`groups') i(`interviews') n(`nonresponse') gen(`generate') ")
	
	ereturn scalar value      = `fvalue'
	ereturn scalar sigmavalue = `sigma'
	ereturn scalar aic        = `aic'
	ereturn scalar schwarz    = `schwarz'
	ereturn scalar n          = `n'
	ereturn scalar ngroups	  = `ngroups'
	
    ereturn display
	
end

mata:

void mywork( string scalar indepvars,  string scalar new_group, 		///
			 string scalar interviews, string scalar nonresponse,       ///
             string scalar bname,      string scalar touse,             ///
			 string scalar w1,         string scalar w2,		    	///
			 string scalar fvname,	   string scalar variance,          ///
			 string scalar vsigma, 	   string scalar vaic,              ///
			 string scalar vschwarz,   string scalar nocons, 			///
			 numeric vector delta,     string scalar start,      		///
			 string scalar difficult,  numeric vector maxiter           ///
 )
{
 
	real scalar f, s
    real vector b
    real matrix X, y, Z1, Z2, groups, addobs, P, pop, g, G, sigma2, CI1, CI2
	string scalar Z2_name, ci_uname, ci_lname
 
	X = st_data(., (new_group,interviews,nonresponse), touse)
	
	if (indepvars == "" & nocons == "")  {
	y = J(rows(X),1,1)
	} 
	else if (indepvars != "" & nocons != "") {
    y = st_data(., indepvars, touse)
	}
	else if (indepvars != "" & nocons == "") {
    y = st_data(., indepvars, touse),J(rows(X),1,1)
	}
	
	groups = rows(mm_collapse(X[.,(2::3)], 1, X[.,1]))
	addobs = J(groups,rows(y),0)
	for (i=1; i<=groups; i++) {
	addobs[i,.] = X[.,1]' :== i
	}
	
	pop = mm_collapse(X[.,(2::3)], 1, X[.,1])
	pop = pop[.,2]+pop[.,3]
		
	S = optimize_init()

	optimize_init_argument(S, 1, y)
    optimize_init_argument(S, 2, X)
	optimize_init_argument(S, 3, pop)
	optimize_init_argument(S, 4, addobs)
	optimize_init_params(S,J(1,cols(y),0))
	if (start!="") {
	optimize_init_params(S,st_matrix(start))
	}
	optimize_init_evaluator(S, &plleval())
	optimize_init_which(S,"min") 
	optimize_init_technique(S, st_local("technique"))
	if (st_local("technique")=="bhhh") {
    optimize_init_evaluatortype(S,"gf0")
	}
	optimize_init_nmsimplexdeltas(S,delta)
	if (difficult!="") {
	optimize_init_singularHmethod(S,"hybrid")
	}
	optimize_init_conv_maxiter(S,maxiter)
	b = optimize(S)
	st_matrix(bname, b)
	
	f = optimize_result_value(S)
	st_numscalar(fvname, f)
	
	
	P = invlogit(y * b')
	g = pop-addobs*(1:/P) 
	s = sum(g :* g) / sum(pop)
	st_numscalar(vsigma, s)
	
	G = addobs * (y :/ (J(1,cols(y),1)#exp(y * b')) )
	sigma2 = s * invsym(G'*diag(1 :/ pop)*G)
	st_matrix(variance, sigma2)
	
	// Return compliance function
	if (w2!="") {
	st_view(Z1, ., st_addvar(("float"), (w2)),touse )
	Z1[., .] = P
	ci_uname = w2+"_upper"
	st_view(CI1, ., st_addvar(("float"), (ci_uname)),touse )
	CI1[., .] = invlogit(y * b' + 1.96*cross(((y*sigma2):*y)',J(cols(y),1,1)) ) 	
	ci_lname = w2+"_lower"
	st_view(CI2, ., st_addvar(("float"), (ci_lname)),touse )
	CI2[., .] = invlogit(y * b' - 1.96*cross(((y*sigma2):*y)',J(cols(y),1,1)) )
	}
	
	// Return corrected survey weights
	if (w1!="") {
	Z2_name = w1+"_c"
	st_view(Z2, ., st_addvar(("float"), (Z2_name)),touse )
	Z2[., .] = st_data(., w1, touse) :/ P
	}
		
	// Return information criteria	
	st_numscalar(vaic, rows(pop)*log(f/rows(pop))+2*cols(y))
	st_numscalar(vschwarz, rows(pop)*log(f/rows(pop))+cols(y)*log(cols(y)))
	
}

void plleval(real scalar todo,   real vector b, ///
             real matrix y,      real matrix X, ///
			 real matrix pop,    real matrix addobs, ///
		     crit, grad, hess)
{
    real matrix P, g
	
	P			= invlogit(y * b')
	g           = pop-addobs*(1:/P) 
	crit 		= quadcross(g, (1:/pop), g)
	
}

end


