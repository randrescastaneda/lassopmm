{smcl}
{* *! version 1.0.0  26Dec2018}{...}
{cmd:help lassopmm}
{hline}

{title:Title}

{p2colset 5 24 26 2}{...}
{p2col :{cmd:lassopmm} {hline 1}} PMM enabled Lasso regression for multiple imputation, command makes use of Wilbur Townsend's lassoregress command {p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 17 2}{opt lassoregress}   {depvar} [{indepvars}] {ifin} {weight} {cmd:,} uniqid(varlist) add(integer) [{it:options}]

{synoptset 20 tabbed}{...}
{synoptline}
{synopthdr}
{synoptline}
{syntab:{title:Required}}

{synopt:{opt uniqid(varlist)}}varlist of unique identifiers - used to ensure replicability{p_end}
{synopt:{opt add(integer)}}number of simulations to be run for multiple imputation{p_end}

{syntab:{title:Optional}}

{synopt:{opt SORTy}}method which takes a random draw of the observed dependent variable, sorts it and matches to the sorted yhat results for the non-observed data. When not specified, the default, pmm is used.{p_end}
{synopt:{opt psu(varlist)}}varlist for psu bootstrapping of data{p_end}
{synopt:{opt lambda(real)}}penalty placed on larger coefficients — by default found by cross-validation {p_end}
{synopt:{opt numfolds(integer)}}number of folds used when cross-validating lambda or alpha — default is 10. {p_end}
{synopt:{opt numlambda(integer)}}number of lambda tested when lambda is found by cross-validation{p_end}
{synopt:{opt knn(integer)}}number of closest observations to draw result from {p_end}
{synopt:{opt seed(integer)}}Seed to ensure replicability, if not specified it uses the current seed's state. {p_end}
{synopt:{opt NOIsily}}Display results of lassoregress{p_end}

{p 2 2 1}
{cmd:aweight}s are allowed; see {help weight}{break} 
{cmd:mi set wide} is allowed; see {help mi set}{break}
{cmd:mi set mlong} is allowed; see {help mi set}{break} 

{title:Example}
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

mi set wide

mi register imputed price 

lassopmm `_y' `_x' [aw=weight], sorty knn(5) add(5) psu(psu) seed(12388) uniqid(_numobs11)
mi estimate: mean price if samples==1 [aw=weight]


{title:Authors:}

{pstd}
Raul Andres Castaneda{break}
The World Bank - Poverty and Equity Global Practice {break}
Washington, DC{break}
acastaneda@worldbank.org{p_end}

{pstd}
Paul Corral{break}
The World Bank - Poverty and Equity Global Practice {break}
Washington, DC{break}
Corresponding author{break} 
pcorralrodas@worldbank.org{p_end}

{pstd}
Leonardo Ramiro Lucchetti{break}
The World Bank - Poverty and Equity Global Practice {break}
Washington, DC{break}
llucchetti@worldbank.org{p_end}

{title:References}

{pstd}
Wilbur Townsend, 2017. "ELASTICREGRESS: Stata module to perform elastic net regression, lasso regression, ridge regression," Statistical Software Components S458397, Boston College Department of Economics, revised 16 Apr 2018.

{title:Disclaimer}

{pstd}
Any error or omission is the author's responsibility alone.

