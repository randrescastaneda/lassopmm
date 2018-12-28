set more off
clear all

run "C:\Users\WB378870\OneDrive - WBG\000.my_ados\lassopmm\lassopmm.ado"

sysuse auto, clear

gen _numobs = _n
preserve
	sample 10
	sum price
	replace price = .
		
	tempfile uno
	save `uno'
restore
append using `uno', gen(samples)


//Local with all candidate variables
local _x mpg headroom trunk weight length turn displacement gear_ratio foreign
//Local with dependent variable
local _y price 

preserve
lassopmm `_y' `_x' [aw=weight], knn(1)
sum `_y' if samples==1
restore 
preserve
lassopmm `_y' `_x' [aw=weight], knn(5)
sum `_y' if samples==1
restore

preserve
lassopmm `_y' `_x' [aw=weight], sort
sum `_y' if samples==1
restore

preserve
lassopmm `_y' `_x' [aw=weight], sort lambda(0)
sum `_y' if samples==1
restore

