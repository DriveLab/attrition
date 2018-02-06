/******************************************************************************
Overview:	This file creates a programme to implement Lee bounds for unbalanced
					attrition, allowing any number of treatment groups. The bounds and
					their associated standard errors are ereturned in the matrices
					b_upper, s_upper, b_lower, s_lower and the code for the reference
					group is returned in refgroup.
					See notes 3 and 4 for important limitations.
					The programme only accommodates the regress command.

Created:	2015/10/08
Updated:	2015/10/08 by Rob to finish programming
					2015/10/22 by Simon to avoid overly long variable names
					2016/03/14 by Rob to deal with bounded variable with mass points at
						one or both boundaries.
					2017/01/31 by Rob to calculate bounds on weekly in-person vs phone

Inputs:		NA
Outputs:	NA

Notes:		1) The programme works for any number of treatment groups and allows
					different levels of attrition in each group.
					2) All groups are trimmed to have the same fraction of missing data
					as the group with the highest level of attrition (the ``reference
					group''). If the reference group is not the control group, then the
					upper bound on the treatment effects for groups I and J may be from
					different regressions. The programme accounts for this and returns
					the correct bounds in the matrices b_upper and b_lower.
					3) The standard errors are not correct, because they do not take into
					account the variance component from the estimated trimming rate. They
					should only be treated as suggestive.
					4) The programme pools survey waves to attain balanced attrition for
					the full panel. It does not attain balanced attrition within each
					wave. This is probably not correct but it would require a much more
					complex code.
					5) The programme requires the following inputs to be specified:
							- scalar outcome variable (can be discrete or continuous)
							- categorical treatment variable (with 0 for the control group)
							- clustering variable
							- one or more control variables (optional)
					6) The programme saves the upper and lower bounds in the matrices
					e(upper) and e(lower) respectively, with the point estimates and
					standard errors in the first and second columns respectively.
					7) The programme sorts the data.
					8) The programme requires that the seed is set to produce consistent
					results, as the sort command breaks ties ``randomly.''

Syntax:		leebound <outcome> <treatment> <clustervar> <covlist>
******************************************************************************/



program define leebound, eclass
	args y treat clust xlist

	/** Calculate upper and lower bounds of outcome distribution. **/
	qui summ `y'
		local tempmax = r(max)
		local tempmin = r(min)

	/** Create placeholder matrices for reference group, point estimates and
			standard errors for coefficients in non-control group. **/
	qui summ `treat'
		local maxtreat = r(max)
		local treatcomp = `maxtreat' - 1
	forvalues i = 1(1)`maxtreat' {
		matrix temp`i' = J(2,2,.)
	}
	matrix upper = J(`maxtreat',2,.)
	matrix lower = J(`maxtreat',2,.)

	/** Create placeholder matrices for point estimates and standard errors for
			comparisons between treatment groups. **/
	forvalues i = 1(1)`treatcomp' {
		local j = min( `i' + 1, `maxtreat' )
		forvalues k = `j'(1)`maxtreat' {
			matrix comp`i'`k' = J(2,2,.)
			}
		}
		
	/** Calculate the non-missing rate in each group, identify the group with the
			lowest non-missing rate (``reference group''), and hence calculate the
			rate of excess non-missing data in the other groups. **/
	local tempnonmiss = 1
	forvalues i = 0(1)`maxtreat' {
		qui count if `treat'==`i'
			local tempden = r(N)
		qui count if `treat'==`i' & `y'!=.
			local tempnonmiss`i' = r(N)/`tempden'
			display "Rate of non-missing data in group " `i' " is " `tempnonmiss`i''
		if `tempnonmiss`i'' < `tempnonmiss' {
			local temprefgroup = `i'
				
			}
		}
	matrix temprefgroup = `temprefgroup'

	/** Calculate the fraction and number of observations in each group other
			than the reference group that should be trimmed. **/
	forvalues i = 0(1)`maxtreat' {
		local temptrimpct`i' = abs(`tempnonmiss'-`tempnonmiss`i'')/`tempnonmiss`i''
		qui count if `treat'==`i' & `y'!=.
		local tempnum`i' = r(N)
		local temptrimnum`i' = round( r(N)*`temptrimpct`i'', 1)
		}
	local tempmiss = 1-`tempnonmiss'

	/** Create trimmed outcome measures, by sorting the data by the outcome and
			trimming the required number of observations from the top/bottom. **/
	qui gen ttt`y' = `y'
	qui gen ttl`y' = `y'
	qui gen temprand = runiform()
	sort `treat' `y' temprand
	qui gen temporder = 1 if _n==1 | `treat'!=`treat'[_n-1]
	qui replace temporder = temporder[_n-1] + 1 if temporder==.
	forvalues i = 0(1)`maxtreat' {
		if `i'!=`temprefgroup' {
			qui replace ttl`y' = . if `treat'==`i' & temporder<`temptrimnum`i''
			qui replace ttt`y' = . if `treat'==`i' & temporder>`tempnum`i''-`temptrimnum`i''
			}
		}
	drop temprand temporder

	/** Generate treatment group dummies. **/
	forvalues i = 1(1)`maxtreat' {
		qui gen temptreat_`i' = (treatment==`i')
		}

	/** Calculate bounds from trimmed outcomes, for control vs treatment and then
			for treatment vs treatment. **/
	qui reg ttt`y' temptreat_* `xlist', cluster(`clust')
	forvalues i = 1(1)`maxtreat' {
		matrix temp`i'[1,1] = _b[temptreat_`i']
		matrix temp`i'[1,2] = _se[temptreat_`i']
		}
	forvalues i = 1(1)`treatcomp' {
		local j = min( `i' + 1, `maxtreat' )
		forvalues k = `j'(1)`maxtreat' {
			qui lincom temptreat_`i' - temptreat_`k'
			matrix comp`i'`k'[1,1] = r(estimate)
			matrix comp`i'`k'[1,2] = r(se)
			}
		}
	qui reg ttl`y' temptreat_* `xlist', cluster(`clust')
	forvalues i = 1(1)`maxtreat' {
		matrix temp`i'[2,1] = _b[temptreat_`i']
		matrix temp`i'[2,2] = _se[temptreat_`i']
		}
	forvalues i = 1(1)`treatcomp' {
		local j = min( `i' + 1, `maxtreat' )
		forvalues k = `j'(1)`maxtreat' {
			qui lincom temptreat_`i' - temptreat_`k'
			matrix comp`i'`k'[2,1] = r(estimate)
			matrix comp`i'`k'[2,2] = r(se)
			}
		}

	/** Save bounds on control vs treatment in permanent matrix. **/
	forvalues i = 1(1)`maxtreat' {
		mata : st_matrix("temp`i'", sort(st_matrix("temp`i'"), 1))
		matrix upper[`i',1] = temp`i'[2,1]
		matrix upper[`i',2] = temp`i'[2,2]
		matrix lower[`i',1] = temp`i'[1,1]
		matrix lower[`i',2] = temp`i'[1,2]
		}

	/** Save bounds on treatment vs treatment in permanent matrix. **/
	forvalues i = 1(1)`treatcomp' {
		local j = min( `i' + 1, `maxtreat' )
		forvalues k = `j'(1)`maxtreat' {
			mata : st_matrix("comp`i'`k'", sort(st_matrix("comp`i'`k'"), 1))
			}
		}

	/** Display analysis results for control vs treatment. **/
	matrix colnames upper = b se
	matrix colnames lower = b se
	local temprowname = ""
	forvalues i = 1(1)`maxtreat' {
		local temprowname = "`temprowname' treatment_`i' "
		}
	matrix rownames upper = `temprowname'
	matrix rownames lower = `temprowname'
	matrix list upper
	matrix list lower

	/** Display analysis results for control vs treatment. **/
	local temprowname = "lower_bound upper_bound"
	forvalues i = 1(1)`treatcomp' {
		local j = min( `i' + 1, `maxtreat' )
		forvalues k = `j'(1)`maxtreat' {
			matrix colnames comp`i'`k' = b se
			matrix rownames comp`i'`k' = `temprowname'
			matrix list comp`i'`k'
			}
		}

	/** Save point estimates and standard errors from trimmed regression. **/
	ereturn clear
	ereturn matrix refgroup temprefgroup
	ereturn matrix upper upper
	ereturn matrix lower lower
	forvalues i = 1(1)`treatcomp' {
		local j = min( `i' + 1, `maxtreat' )
		forvalues k = `j'(1)`maxtreat' {
			ereturn matrix comp`i'`k' comp`i'`k'
			}
		}
	capture drop ttt`y' ttl`y' temptreat_*

end



