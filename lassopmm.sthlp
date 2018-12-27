{smcl}
{* *! version 1.0.0  26Dec2018}{...}
{cmd:help lassopmm}
{hline}

{title:Title}

{p2colset 5 24 26 2}{...}
{p2col :{cmd:lassopmm} {hline 1}} PMM enabled Lasso regression for multiple imputation, command makes use of Wilbur Townsend's lassoregress command {p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 17 2}{opt lassoregress}   {depvar} [{indepvars}] {ifin} {weight} [{cmd:,} {it:options}]

{synoptset 15 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt lambda}}penalty placed on larger coefficients — by default found by cross-validation; {p_end}
{synopt:{opt numfolds}}number of folds used when cross-validating lambda or alpha — default is 10. {p_end}
{synopt:{opt numlambda}}number of lambda tested when lambda is found by cross-validation;{p_end}
{synopt:{opt numsim}}number of simulations to be run for multiple imputation{p_end}
{synopt:{opt knn}}Number of closest observations to draw result from {p_end}
{synopt:{opt SORTy}}Method which takes a random draw of the observed dependent variable, sorts it and matches to the sorted yhat results for the non-observed data.{p_end}
{synopt:{opt seed}}Seed to ensure replicability, if not specified it uses the current seed's state. {p_end}

{p 4 6 2}
{cmd:aweight}s are allowed; see {help weight}.






