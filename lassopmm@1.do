set more off
clear all

mata
//Function will return an index selection vector for Y
function _Mpmm(yh1, yh2, knn){
	//Search distance to yh2
	ry2 = rows(yh2)
	ry1 = rows(yh1)
	if (knn>1){
		for(i=1; i<=ry2; i++){
			myy = order(abs(yh1:-yh2[i]),1)[1..knn]	
			
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


//n is the number of simulations, dim is the number of rows of the output, epsi is the source
function _f_sampleepsi(real scalar n, real scalar dim, real matrix eps) {				  
	sige2 = J(dim,n,0)
	N = rows(eps)
	if (cols(eps)==1) for(i=1; i<=n; i++) sige2[.,i]=eps[ceil(N*runiform(dim,1)),1]
	else              for(i=1; i<=n; i++) sige2[.,i]=eps[ceil(N*runiform(dim,1)),i]
	//for(i=1; i<=n; i++) sige2[.,i]=eps[ceil(rows(eps)*runiform(dim,1)),i]
	return(sige2)	
}
end

sysuse auto, clear

gen _numobs = _n
preserve
	sample 10
	replace price = .
	
	tempfile uno
	save `uno'
restore
append using `uno', gen(samples)

tab samples

//Local with all candidate variables
local _my_x mpg headroom trunk weight length turn displacement gear_ratio foreign
//Local with dependent variable
local depvar price 

//Lasso regress
lassoregress `depvar' `_my_x' [aw=weight]
local myvar =  e(varlist_nonzero)
	tempname _beta BB
	tempvar touse touse1 touse2
	
	mat `_beta' = e(b)
	
gen `touse1' = e(sample)
gen `touse2'= `touse1'==0

	local a=1
	foreach x of local _my_x{
		if (`_beta'[1,`a']!=0){
			local chosen `chosen' `x'
			mat `BB' = nullmat(`BB'),`_beta'[1,`a']
		}
		local ++a
	}
	mat `BB' = `BB',`_beta'[1,`a']
	
matlist `_beta'
matlist `BB'

mata: b = st_matrix("`BB'")
mata: st_view(x=., .,"`chosen'", "`touse1'")
mata: st_view(y=., .,"`depvar'", "`touse1'")
mata: st_view(y1=.,.,"`depvar'", "`touse2'")
mata: st_view(x1=.,.,"`chosen'", "`touse2'")

mata: yhat1 = quadcross((x , J(rows(x),1,1))',b')
mata: yhat2 = quadcross((x1,J(rows(x1),1,1))',b')

mata: y1[.,.]=y[_Mpmm(yhat1, yhat2, 1)]

end
