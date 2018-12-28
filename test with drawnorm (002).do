
clear all
set more off

/*
Author: Paul Corral
Source data from Nobuo Yoshida's SWIFT Summer University course 
*/


*Specify the folder where the base year data (y0) is stored
*==========================================================

global swift "C:\Users\WB378870\OneDrive - WBG\Summer U\swift"


*Open a log file
*===============

capture log close
log using "$swift\log_simulationandestimation", replace



*Define globals for regression
*=============================

global location urban region2 region4 region5 region6 region7 region8 region9 region10 
global demographics hhsize head_male head_age marital1 marital2 marital3 marital5 marital6 depratio 
global education head_schooling2 head_schooling3 head_schlvl2 head_schlvl3 head_schlvl4 head_schlvl5 noschoolingp read write
global employment head_employed2 head_employed3 employedp 

#delimit ;
global dwelling rooms owned tenure2 tenure3 water_drinking2 water_drinking3 lighting2 lighting3 lighting4 fuel2 fuel3 fuel4 
solidwaste2 solidwaste3 solidwaste4 toilet2 toilet3 toilet4 toilet5 wall2 wall3 floor2 floor3 roof2 roof3 roof4;
#delimit cr

global assets furniture fridge fan radio desktop laptop dish tv iron bicycle mcycle house land microwave rcooker


*Impute from y0 to y1
*============================

use "$swift\y0", clear
qui append using "$swift\y1", gen(round)
keep if urban==1		//Use when running urban/rural model
keep if region==8		//Use when running regional model



stepwise, pr(.0200001) pe(.02): reg lnrpcexp $location $demographics $education $employment $dwelling $assets susu [aw=pwt]
matrix X = e(b)
matrix X = X[1,1..`e(df_m)']
global xvars: colnames X	

*Run Multiple Imputation

//mi set mlong
* I commented this out since I want to keep it in wide
reg lnrpcexp $xvars [pw = pwt]
predict therr, r
gen noused=!(e(sample))

	//Send residual vector to mata
	putmata non_n_eps = therr if noused==0
 
	//degrees of freedom
	local df=e(df_r)
	//Root mean squared erro
	local eps =e(rmse)
	
	//Send rmse to mata
	mata: err=strtoreal(st_local("eps"))
	//Send deg. freedom to mata
	mata: df=strtoreal(st_local("df"))
	//Draw Sigma 2 like in MI methods and formulas from Stata Corp.
	mata: sig2 = (err^2)*(df):/rchi2(1,20,df)
 

 /*
 
Linear regression                               Number of obs     =        507
                                                F(11, 495)        =      28.89
                                                Prob > F          =     0.0000
                                                R-squared         =     0.4472
                                                Root MSE          =     .45419

------------------------------------------------------------------------------
             |               Robust
    lnrpcexp |      Coef.   Std. Err.      t    P>|t|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
      hhsize |  -.0823107   .0074718   -11.02   0.000    -.0969911   -.0676304
      mcycle |   .1478375   .0509937     2.90   0.004     .0476466    .2480283
 solidwaste2 |  -.3368129   .1413745    -2.38   0.018     -.614581   -.0590448
 solidwaste4 |  -.7171317   .1089622    -6.58   0.000    -.9312172   -.5030462
         fan |   .2340242   .0487921     4.80   0.000      .138159    .3298895
 solidwaste3 |  -.4749995     .09862    -4.82   0.000    -.6687649   -.2812342
     toilet5 |   .1133193    .049202     2.30   0.022     .0166488    .2099899
        dish |   .2210782   .0864841     2.56   0.011      .051157    .3909993
        land |   .2038748   .0508539     4.01   0.000     .1039588    .3037909
       roof2 |   .1694612   .0766482     2.21   0.027     .0188653    .3200571
      laptop |   .3339867   .1734012     1.93   0.055    -.0067065    .6746799
       _cons |   8.102768   .1125103    72.02   0.000     7.881711    8.323825
------------------------------------------------------------------------------

 */
preserve

	//Indicate the variables from the regression
	foreach x in $xvars _cons{
		local imps `imps' `x'_imp
	}
	clear
	
	
	//SPECIFY NUMBER OF SIMS
	local sims = 20
	//Draw the betas from the OLS
	drawnorm `imps', cov(e(V)) means(e(b)) n(`sims')
	
	//Send the beta vector to mata
	putmata betas = (`imps')
	//Transpose the betas
	mata: betas=betas'
restore

//Send the X matrix to mata
putmata X=($xvars) if noused==1
	//Count the number of observations
	mata: rows(X)
	
//Send weight vector to mata
putmata w=pwt if noused==1
//Get 20 vectors of Yhat
mata: X=(X,J(rows(X),1,1))*betas
// Get 20 simulated Y (includes error terms, like in MI)
mata: X1=rnormal(1,1, X,sqrt(sig2))
//Indicate poverty per simulation
mata: poor=mean((X1:<=ln(1314)),w)
// Simulated poverty rate for exercise
mata: mean(poor')


	//Note the above considers normal errors, lets draw from the pool of errors
	// I need a function to get random draws
	mata
	function _f_sampleeps(real scalar n, real scalar dim, real matrix eps) {				  
		sige2 = J(dim,n,0)
		N = rows(eps)
		if (cols(eps)==1) for(i=1; i<=n; i++) sige2[.,i]=eps[ceil(N*runiform(dim,1)),1]
		else              for(i=1; i<=n; i++) sige2[.,i]=eps[ceil(N*runiform(dim,1)),i]
		return(sige2)	
	}
	//Remember X is our YHAT
	residuals = _f_sampleeps(`sims', rows(X) ,  non_n_eps)
	
	X1 = X+residuals
	poor=mean((X1:<=ln(1314)),w)
	// Simulated poverty rate for exercise
	mean(poor')

end

	
	
