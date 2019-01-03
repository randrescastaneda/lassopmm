*! version 0.2 December27-2018
*! R. Andres Castaneda - acastaneda@worldbank.org
*! Paul Corral - pcorralrodas@worldbank.org  -- Corresponding author
*! Leonardo Lucchetti - llucchetti@worldbank.org
*! World Bank Group - Poverty and Equity Global Practice 
*! Equity policy lab

cap program drop lassopmm
program define lassopmm, eclass byable(recall)
	version 11, missing
#delimit;
	syntax varlist(min=2 numeric fv) [if] [in] [aw],
		uniqid(varlist)
		add(integer)
	[
		SORTy
		knn(numlist max=1 int)
		lambda(real -1)
		numlambda(integer 100)
		numfolds(integer 10)
		psu(varlist max=1 numeric)
		seed(integer 12345)
		NOIsily
	];
#delimit cr
set more off

qui{
*===============================================================================
// 1. HOUSE KEEPING
*===============================================================================

	//Mlong?
	local mlong `_dta[_mi_style]'
	
	//confirm mi registered
	if ("`_dta[_mi_marker]'"==""){
		dis as error "You must mi set your data before using this command, see help file"	
		error 119
		exit
	}
	
	//mark sample for use
	marksample touse
	
	// Mark for the other variables in the estimation
	foreach j of varlist `psu' `uniqid'{
		replace `touse' = 0 if missing(`j')
	}
	
	
	tokenize `varlist'
	// Local for dependent variable
	local depvar `1'
	
	// obtain the independent variables
	macro shift 
	local _my_x `*'
	
	if ("`_dta[_mi_ivars]'"!="`depvar'"){
		dis as error "You must mi register `depvar' as imputed before using this command"	
		error 198
		exit
	}	
	
	if ("`sorty'"!="" & "`knn'"!=""){
		//Only display warning!
		dis as error "Note that KNN option only works when the default option of PMM is used."
	}
	
	//Weights for estimation
	tempvar wi _psu _group
	if missing(`"`weight'"') generate byte `wi' = 1
	else generate `wi' `exp'
	
	//Check unique id
	isid `uniqid' 
	
	//PSU for bootstrapping
	if ("`psu'"==""){
		sort `uniqid'
		gen `_group' = 1
	}
	else{
		sort `psu' `uniqid'
		clonevar `_group' = `psu'		
	}
	
*===============================================================================
// 2. Get bootstrap permutation vectors
*===============================================================================
	set seed `seed'	
	mata: st_view(mypsu=.,.,"`_group'", "`touse'")
	noi: mata: ind = _getmyboot(mypsu,`add')	
	
	tempfile mydata
	save `mydata'
	
	
	
*===============================================================================
// 3. Run lasso and imputation on different bootstrapped data
*===============================================================================
	//Get numsim vectors of Y	
	forval i=1/`add'{
		if (`i'!=1) use `mydata', clear
			mata: st_view(myvar=., ., "`depvar' `_my_x' `wi'","`touse'")
			mata: myvar[.,.] = myvar[ind[.,`i'],.]
		if ("`noisily'"!="") noi:lassoregress `depvar' `_my_x' if `touse'==1 [aw=`wi'], numlambda(`numlambda') numfolds(`numfolds') lambda(`lambda')
		else                 lassoregress `depvar' `_my_x' if `touse'==1 [aw=`wi'], numlambda(`numlambda') numfolds(`numfolds') lambda(`lambda')
		tempname _beta BB
		tempvar touse1 touse2
		
		mat `_beta' = e(b)
		
		gen `touse1' = e(sample)
		gen `touse2'= `touse1'==0 & missing(`depvar')
	
		local a=1
		local chosen
		foreach x of local _my_x{
			if (`_beta'[1,`a']!=0){
				local chosen `chosen' `x'
				mat `BB' = nullmat(`BB'),`_beta'[1,`a']
			}
			local ++a
		}
		mat `BB' = `BB',`_beta'[1,`a']
		
		foreach x of local chosen{
			replace `touse2' = 0 if missing(`x')
		}
		
			mata: b = st_matrix("`BB'")
			mata: st_view(x=., .,"`chosen'", "`touse1'")
			mata: st_view(w=., .,"`wi'", "`touse1'")
			mata: sd = quadmeanvariance(x,w)
			mata: mu = sd[1,.]
			mata: sd = sqrt(diagonal(sd[|2,1\.,.|]))'
			
			if (`i'==1)	mata: st_view(y=., .,"`depvar'", "`touse1'")

			
			mata: st_view(x1=.,.,"`chosen'", "`touse2'")
			mata: st_view(w1=., .,"`wi'", "`touse2'")
			
			
			mata: yhat1 = quadcross((((x:-mu):/sd) ,J(rows(x), 1,1))',b')
			mata: yhat2 = quadcross((((x1:-mu):/sd),J(rows(x1),1,1))',b')
	
		if (`i'==1){
			if ("`sorty'"=="") mata: y1 = y[_Mpmm(yhat1, yhat2, `knn')]
			else               mata: y1 = _randomleo(y,yhat2)
		}
		else{
			if ("`mlong'"=="mlong"){
				if ("`sorty'"=="") mata: y1 = y1 \	y[_Mpmm(yhat1, yhat2, `knn')]
				else               mata: y1 = y1 \	_randomleo(y,yhat2)
			}	
			else{	
				if ("`sorty'"=="") mata: y1 = y1 ,	y[_Mpmm(yhat1, yhat2, `knn')]
				else               mata: y1 = y1 ,	_randomleo(y,yhat2)			
			}
		}
	} //End of sim loop
*===============================================================================
// 3. Export simulations to mi set data
*===============================================================================
	if ("`mlong'"=="mlong"){
		use `mydata', clear
			tempvar touse2 nn
			gen `touse2' = `touse'==0 & missing(`depvar')
			qui:count 
			local hu0 = r(N)
				preserve
					qui:keep if `touse2'==1
					tempfile _nnio
					qui:save `_nnio'
					qui:count
					local hu=r(N)
				restore
			local a=1
			while (`a'<=`add'){
				qui:append using `_nnio', gen(`nn')
				qui:replace _mi_m=`a' if `nn'==1
				local ++a
				drop `nn'			
			}
			qui: char _dta[_mi_M] `add'
			qui: char _dta[_mi_n] `hu'
			qui: char _dta[_mi_N] `hu0'		
			replace `touse2'=`touse2'==1 & _mi_m>0
			mata: st_store(.,st_varindex(tokens("`depvar'")),"`touse2'",y1)				
	}
	else{
		use `mydata', clear
			tempvar touse2 nn
			gen `touse2' = `touse'==0 & missing(`depvar')
			replace _mi_miss = 1 if `touse2'==1
			local a=1
			local toimp
			while (`a'<=`add'){
				gen double _`a'_`depvar' = `depvar' if `touse'==1
				
				local toimp `toimp' _`a'_`depvar'
				local ++a				
			}
			qui: char _dta[_mi_M] `add'
			mata: st_store(.,st_varindex(tokens("`toimp'")),"`touse2'",y1)				

	}	
}
end

*===============================================================================
// Annex: MATA functions
*===============================================================================
mata
	//Function will return an index selection vector for PMM Y
	function _Mpmm(yh1, yh2, knn){
		//Search distance to yh2
		ry2 = rows(yh2)
		ry1 = rows(yh1)
		if (knn>1){
			for(i=1; i<=ry2; i++){
				myy = order(abs(yh1:-yh2[i]),1)[|1\knn|]				
				if (i==1) mynn = _f_sampleepsi(1,1,myy)
				else mynn = mynn \ _f_sampleepsi(1,1,myy)		
			}
		}
		else{
			for(i=1; i<=ry2; i++){
				if (i==1) mynn = order(abs(yh1:-yh2[i]),1)[1]	
				else mynn = mynn \ order(abs(yh1:-yh2[i]),1)[1]	
			}	
		}
		return(mynn)
	}	
	//Function, selects a random set of observed y, of length rows unobserved. 
	// sorts the set and assigns the value to the sorted xb from unobserved
	function _randomleo(yo, yh2){
		ry = rows(yh2)
		//random sample of y
		tosort      = sort((runningsum(J(ry,1,1)),yh2),2)
		tosort[.,2] = sort(_f_sampleepsi(1, ry, yo),1)
		_sort(tosort,1)
		
		return(tosort[.,2])
	}	
	//n is the number of simulations, dim is the number of rows of the output, epsi is the source
	function _f_sampleepsi(real scalar n, real scalar dim, real matrix eps){				  
		sige2 = J(dim,n,0)
		N = rows(eps)
		if (cols(eps)==1) for(i=1; i<=n; i++) sige2[.,i]=eps[ceil(N*runiform(dim,1)),1]
		else              for(i=1; i<=n; i++) sige2[.,i]=eps[ceil(N*runiform(dim,1)),i]
		//for(i=1; i<=n; i++) sige2[.,i]=eps[ceil(rows(eps)*runiform(dim,1)),i]
		return(sige2)	
	}
	//Function for bootstrap permutation vectors
	function _getmyboot(_psu,sim){
		info = panelsetup(_psu,1)
		_index = runningsum(J(rows(_psu),1,1))
		rr = rows(info)
		
		for(i=1; i<=rr; i++){
			m = info[i,1],1 \ info[i,2],1
			r1= rows(_index[|m|])
			if (i==1) _X = _f_sampleepsi(sim,r1,_index[|m|])
			else      _X = _X \ _f_sampleepsi(sim,r1,_index[|m|])
		}
		
		return(_X)
	}	
end



