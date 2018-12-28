set more off
clear all

run "C:\Users\WB378870\OneDrive - WBG\000.my_ados\lassopmm\lasso2.ado"

sysuse auto, clear

gen psu = 7 if _n<74
replace psu = 6 if _n<60
replace psu = 5 if _n<50
replace psu = 4 if _n<40
replace psu = 3 if _n<30
replace psu = 2 if _n<20
replace psu = 1 if _n<10

gen _numobs = _n
preserve
	sample 10
	sum price
	replace price = .
		
	tempfile uno
	save `uno'
restore
append using `uno', gen(samples)

gen _numobs11 = _n

//Local with all candidate variables
local _x mpg headroom trunk weight length turn displacement gear_ratio foreign
//Local with dependent variable
local _y price 

preserve
set trace on
set traced 1
lassopmm `_y' `_x' [aw=weight], knn(5) numsim(5) psu(psu) seed(12388) uniqid(_numobs11)
sum `_y' if samples==1

sss
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

