*! version 0.1 December27-2018
*! Works only for one vector of pmm
*! Raul Andres Castaneda - acastaneda@worldbank.org
*! Paul Corral - pcorralrodas@worldbank.org 
*! Leonardo Lucchetti - llucchetti@worldbank.org
*! World Bank Group - Poverty and Equity Global Practice 
*! Equity policy lab

cap program drop lassopmm
program define lassopmm, eclass byable(recall)
	version 11, missing
#delimit;
	syntax varlist(min=2 numeric fv) [if] [in] [aw],
		uniqid(varlist)
	[
		SORTy
		knn(numlist max=1 int)
		lambda(real -1)
		numlambda(integer 100)
		numsim(integer 1)
		numfolds(integer 10)
		psu(varlist max=1 numeric)
		seed(integer 12345)
		mlong
	];
#delimit cr
set more off

qui{
	
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
	
	if ("`psu'"==""){
		sort `uniqid'
		gen `_group' = 1
	}
	else{
		sort `psu' `uniqid'
		clonevar `_group' = `psu'		
	}
	
	mata: st_view(mypsu=.,.,"`_group'", "`touse'")
	noi: mata: ind = _getmyboot(mypsu,`numsim')	
	
	tempfile mydata
	save `mydata'
	
	set seed `seed'	
	
	//Get numsim vectors of Y	
	forval i=1/`numsim'{
		if (`i'!=1) use `mydata', clear
			noi:mata: st_view(myvar=., ., "`depvar' `_my_x' `wi'","`touse'")
			mata: myvar[.,.] = myvar[ind[.,`i'],.]
		noi:lassoregress `depvar' `_my_x' if `touse'==1 [aw=`wi'], numlambda(`numlambda') numfolds(`numfolds') lambda(`lambda')
		tempname _beta BB
		tempvar touse1 touse2
		
		mat `_beta' = e(b)
		
		gen `touse1' = e(sample)
		gen `touse2'= `touse1'==0
	
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
			mata: st_view(y=., .,"`depvar'", "`touse1'")
			mata: st_view(x1=.,.,"`chosen'", "`touse2'")
			
			mata: yhat1 = quadcross((x , J(rows(x),1,1))',b')
			mata: yhat2 = quadcross((x1,J(rows(x1),1,1))',b')
	
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
				if ("`sorty'"=="") mata: y1 = y1 \	y[_Mpmm(yhat1, yhat2, `knn')]
				else               mata: y1 = y1 \	_randomleo(y,yhat2)			
			}
		}	
	}
	
	if ("`mlong'"=="mlong"){
		use `mydata', clear
		
		
		
	}
	else{
		use `mydata', clear
	}	
}
end

//MATA functions
mata
	//Function will return an index selection vector for Y
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
	//Function for bootstrapping
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



