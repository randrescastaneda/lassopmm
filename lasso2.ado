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
	];
#delimit cr
set more off

qui{
	tokenize `varlist'
	// Local for dependent variable
	local depvar `1'
	
	// obtain the independent variables
	macro shift 
	local _my_x `*'
	
	if ("`sorty'"!="" & "`knn'"!=""){
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
		clonevar `_group' = psu		
	}
	
	tempfile mydata
	save `mydata'
	
	set seed `seed'
	if (`numsim'>1){
		//Get numsim vectors of Y
		
		
		
	}
