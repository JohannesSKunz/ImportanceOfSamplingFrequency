* Packages 
*ssc install xtspj  

clear all 
set more off

* Insert path to data

***
*** Data cleaning and preparation
***
***

collapse lifesat sah famsat incsat econsat hlthsat jobsat ragedobok gender, ///
			by(prim_key wave)

xtset prim_key wave

g yearly=mod(wave,12)
g quarterly=mod(wave,3)

g year =  2015 * (wave<=4) ///
		+ 2016 * (wave>4  & wave<=16) ///
		+ 2017 * (wave>16 & wave<=28) ///
		+ 2018 * (wave>28 )

g month = 1 * (wave==5  | wave==17 | wave==29) ///
		+ 2 * (wave==6  | wave==18 | wave==30) ///
		+ 3 * (wave==7  | wave==19 | wave==31) ///
		+ 4 * (wave==8  | wave==20 | wave==32) ///
		+ 5 * (wave==9  | wave==21 | wave==33 | wave==-3) ///
		+ 6 * (wave==10 | wave==22 | wave==34 | wave==-2) ///
		+ 7 * (wave==11 | wave==23 | wave==35 | wave==-1) ///
		+ 8 * (wave==12 | wave==24 | wave==36 | wave==0) ///
		+ 9 * (wave==13 | wave==25 | wave==37 | wave==1) ///
		+10 * (wave==14 | wave==26 | wave==38 | wave==2) ///
		+11 * (wave==15 | wave==27 | wave==39 | wave==3) ///
		+12 * (wave==16 | wave==28 | wave==40 | wave==4) 

ren ragedobok age		
replace wave = wave + 4 

sort prim_key wave



***
*** AR(1) Pooled OLS without sample restrictions
***
***
		
qui reg lifesat L.lifesat i.year i.age i.gender	, cluster(prim_key)
	est sto pool_m_bigsample

qui reg lifesat L3.lifesat i.year i.age i.gender, cluster(prim_key)
	est sto pool_q1_bigsample
qui reg lifesat L3.lifesat i.year i.age i.gender ///
		if quarterly==0, cluster(prim_key)
	est sto pool_q2_bigsample
	
qui reg lifesat L12.lifesat i.year i.age i.gender, cluster(prim_key)
	est sto pool_y1_bigsample
qui reg lifesat L12.lifesat i.year i.age i.gender ///
		if yearly==0, cluster(prim_key)
	est sto pool_y2_bigsample
	

	
***
*** Data cleaning and preparation continued...
***
***	
	

drop if wave<5
drop if wave==42 

bys prim_key: g ti = _n
bys prim_key: g Ti = _N
bys prim_key: g gap = wave[_n]-wave[_n-1]		
replace gap = 1 if gap==.
bys prim_key: egen hasgap = max(gap)

keep if hasgap==1
keep if Ti==37

replace lifesat = -1 if lifesat==.
bysort prim_key: egen minls = min(lifesat)
drop if minls==-1
drop if age==.

su age year gender lifesat

***
*** 1. AR(1) Pooled OLS
***
***
		
qui reg lifesat L.lifesat i.year i.age i.gender	, cluster(prim_key)
	est sto pool_m

qui reg lifesat L3.lifesat i.year i.age i.gender, cluster(prim_key)
	est sto pool_q1
qui reg lifesat L3.lifesat i.year i.age i.gender ///
		if quarterly==0, cluster(prim_key)
	est sto pool_q2
	
qui reg lifesat L12.lifesat i.year i.age i.gender, cluster(prim_key)
	est sto pool_y1
qui reg lifesat L12.lifesat i.year i.age i.gender ///
		if yearly==0, cluster(prim_key)
	est sto pool_y2	


	
***
*** 2. AR(1) FE 
***
***
		
qui xtreg lifesat L.lifesat i.year i.age i.gender, fe cluster(prim_key)
	est sto fe_m

qui xtreg lifesat L3.lifesat i.year i.age i.gender, fe cluster(prim_key)
	est sto fe_q1
qui xtreg lifesat L3.lifesat i.year i.age i.gender ///
		if quarterly==0, fe cluster(prim_key)
	est sto fe_q2
	
qui xtreg lifesat L12.lifesat i.year i.age i.gender, fe cluster(prim_key)
	est sto fe_y1
qui xtreg lifesat L12.lifesat i.year i.age i.gender ///
		if yearly==0, fe cluster(prim_key)
	est sto fe_y2	

	
** Results from AR(1) models, pooled and fe	

esttab pool* , b(3) se keep(*lifesat*) ///
	rename(L3.lifesat L.lifesat L12.lifesat L.lifesat) ///
	mtitle("monthly" "quarterly" "quarterly" "yearly" "yearly")	
	
esttab fe* , b(3) se keep(*lifesat*) ///
	rename(L3.lifesat L.lifesat L12.lifesat L.lifesat) ///
	mtitle("monthly" "quarterly" "quarterly" "yearly" "yearly")	
	
	
***	
*** 3. FE without covariates ("no x")
***
***
	
g Llifesat = L.lifesat
global lags Llifesat
forval t=2/12 {
	g L`t'lifesat = L`t'.lifesat
	global lags $lags L`t'lifesat
	}

	
qui xtreg lifesat $lags , fe cluster(prim_key)
	est sto all_nox
	g esample = e(sample)
	mat bhat = e(b)
	mat Vhat = e(V)
	
qui xtspj lifesat $lags , model(regress) method(parm) cluster(prim_key)
	est sto all_nox_bc	
	
			
sort prim_key esample wave
bys prim_key esample (wave): g ti2 = _n 
bys prim_key esample: g Ti2 = _N 

gen odd = mod(Ti2,2)

gen half = Ti2/2 if odd==0
	set seed 38495902 	  
	g u = runiform()  	 
	qui su u			 
	replace u = (u-r(mean))/2  
	bysort prim_key: replace u = u[1]  
	replace half = round(Ti2/2+u,1) if odd==1
gen first_half = esample & ti2<=half
replace first_half=. if esample==0

sort prim_key wave	

qui xtreg lifesat  $lags if first_half, fe cluster(prim_key)
	est sto all_nox_1
	mat b1hat = e(b)

qui xtreg lifesat  $lags if first_half==0, fe cluster(prim_key)
	est sto all_nox_2
	mat b2hat = e(b)	

mat b = 2*bhat - 0.5*(b1hat+b2hat) // spj bias-corrected estimate

est rest all_nox
erepost b = b 
est sto all_nox_manual
	

esttab  all_nox all_nox_bc all_nox_manual, b(3) se keep(*lifesat*) ///
	mtitle("fe"  "spj" "manual spj" )	
	

	
***
*** 4. FE and covariates (age, year, gender)
***
***
	
	
*** Program to automate FE with SPJ bias-correction
	
cap program drop spj_manual
program define spj_manual
        version 15.1
        syntax , yvar(string ) xvar(string) eststo(string ) [bsize(integer 1)]
        
		cap drop esample 
		cap drop ti2 
		cap drop Ti2 
		cap drop odd 
		cap drop half 
		cap drop u 
		cap drop first_half	

		qui xtreg `yvar' `xvar' , fe cluster(prim_key)
			est sto `eststo'
			qui g esample = e(sample)
			mat bhat = e(b)
			mat bhat = bhat[1,1..`bsize']
			mat Vhat = e(V)
			mat Vhat = Vhat[1..`bsize',1..`bsize']
			
		sort prim_key esample wave
		bys prim_key esample (wave): g ti2 = _n 
		bys prim_key esample: g Ti2 = _N 
		qui gen odd = mod(Ti2,2)		 
		qui gen half = Ti2/2 if odd==0	 
		set seed 38495902				 
		qui g u = runiform()			 
		qui su u						 
		qui replace u = (u-r(mean))/2	  
		bysort prim_key: replace u = u[1]  
		qui replace half = round(Ti2/2+u,1) if odd==1  
		qui gen first_half = esample & ti2<=half   
		qui replace first_half=. if esample==0  
		sort prim_key wave				 
		
		qui xtreg `yvar' `xvar' if first_half, fe cluster(prim_key)
			est sto `eststo'1
			mat b1hat = e(b)
			mat b1hat = b1hat[1,1..`bsize']

		qui xtreg `yvar' `xvar' if first_half==0, fe cluster(prim_key)
			est sto `eststo'2
			mat b2hat = e(b)
			mat b2hat = b2hat[1,1..`bsize']

		mat b = 2*bhat - 0.5*(b1hat+b2hat)
		est rest `eststo'
		erepost b = b V = Vhat
		est sto `eststo'_manual
	
		esttab  `eststo' `eststo'_manual, b(3) se keep(*`yvar'*) ///
			mtitle("fe" "spj (manual)" )	
	
    end
	
	
*** 4.1 Monthly AR(12) FE with SPJ bias-correction

spj_manual, yvar(lifesat) xvar($lags i.age i.year i.gender) ///
	eststo(all) bsize(12)
	
*** 4.2 Monthly AR(1) FE with SPJ bias-correction	

spj_manual, yvar(lifesat) xvar(Llifesat i.age i.year i.gender) ///
	eststo(m_ar1_fe)

*** 4.3 Quarterly AR(1) FE with SPJ bias-correction	

spj_manual, yvar(lifesat) xvar(L3lifesat i.age i.year i.gender) ///
	eststo(q_ar1_fe)

*** 4.4 Yearly AR(1) FE with SPJ bias-correction	

spj_manual, yvar(lifesat) xvar(L12lifesat i.age i.year i.gender) ///
	eststo(y_ar1_fe)
	
	
***	
*** TABLES
***
***


** Table 1: AR(1), pooled OLS

esttab pool_y1 pool_y2 pool_q1 pool_q2 pool_m , b(3) se keep(*lifesat*) ///
	rename(L3.lifesat L.lifesat L12.lifesat L.lifesat) ///
	mtitle( "yearly" "yearly" "quarterly" "quarterly" "monthly")	

esttab pool_y1 pool_y2 pool_q1 pool_q2 pool_m ///
	using "~/Desktop/table_pool.tex"  ///
	, tex replace b(3) se keep(*lifesat*) ///
	rename(L3.lifesat L.lifesat L12.lifesat L.lifesat) ///
	mtitle( "yearly" "yearly" "quarterly" "quarterly" "monthly")	
	
 
	
** Table 2: AR(1), FE and bias-corrected FE

esttab 	y_ar1_fe y_ar1_fe_manual ///
		q_ar1_fe q_ar1_fe_manual ///
		m_ar1_fe m_ar1_fe_manual ///
		, b(3) se keep(*lifesat*) ///
		rename(L3lifesat Llifesat L12lifesat Llifesat) ///
		mtitle( "yearly" "yearly bc" "quarterly" "quart., bc" "monthly" "monthly, bc")

esttab 	y_ar1_fe y_ar1_fe_manual ///
		q_ar1_fe q_ar1_fe_manual ///
		m_ar1_fe m_ar1_fe_manual ///
		using "~/Desktop/table_fe.tex"  ///
		, tex replace b(3) se keep(*lifesat*) ///
		rename(L3lifesat Llifesat L12lifesat Llifesat) ///
		mtitle( "yearly" "yearly bc" "quarterly" "quart., bc" "monthly" "monthly, bc")
		
   
	
*** Figure 1: AR(12), FE and bias-corrected FE	

est resto m_ar1_fe_manual
	mat b = e(b)
	sca rho = b[1,1]
	forval i=1/12 {
		nlcom  _b[Llifesat]^`i'
		mat V = r(V)
		sca var`i' = V[1,1]
		}
est resto all_manual
	mat bneu = e(b)
	mat Vneu = e(V)
	forval i=1/12 {
		mat bneu[1,`i'] = rho^`i'
		mat Vneu[`i',`i'] = var`i'
		}
erepost b = bneu V=Vneu 
est sto all_ar1
	
coefplot 	(all_ar1,    label("AR(1)") lw(thick) lp(dash) ) ///
			(all_manual, label("AR(12)") lw(thick)  ), ///
			keep(*lifesat*) vertical yline(0, lcolor(black)) ///
			coeflabel(Llifesat = "1" L2lifesat = "2" L3lifesat = "3" ///
				L4lifesat = "4" L5lifesat = "5" L6lifesat = "6" ///
				L7lifesat = "7" L8lifesat = "8" L9lifesat = "9" ///
				L10lifesat = "10" L11lifesat = "11" L12lifesat = "12" ) ///
		levels(99) cismooth mlabel format(%9.2g) mlabposition(1) mlabgap(*2) ///
		ytitle("well-being") xtitle("months") recast(line) ///
		ysize(4.5) xsize(5) ///
		nooffset ///
		 graphregion(color(white) fcolor(white)) ///
		 legend(ring(0) position(2) cols(1) region(lstyle(none)) bmargin(medlarge)) ///
		 ylabel(, glcolor(gs14) glwidth(thin) )
		 
graph export "~/fig1_r1.pdf", replace	
		
 
esttab  all all_manual all_ar1, b(3) se keep(*lifesat*) ///
	mtitle("fe" "spj (manual)" "based on ar(1)" )	 

esttab  all all_manual all_ar1 ///
	using "~/Desktop/table_ar12.tex", tex b(3) se keep(*lifesat*) ///
	mtitle("fe" "spj (manual)" "based on ar(1)" ) replace	

	
*** Table 3: Domain satisfaction outcomes


cap erase "~/Desktop/table_domain.tex"
foreach yvar in famsat incsat econsat hlthsat jobsat {
	
	* (1) Pooled, yearly
	qui reg `yvar' L12.`yvar' i.year i.age i.gender, cluster(prim_key)
		est sto `yvar'_pool
		
	* (2) Bias-corrected FE, yearly
	qui spj_manual, yvar(`yvar') xvar(L12.`yvar' i.age i.year i.gender) ///
	eststo(`yvar'_y)
	
	* (3) Bias-corrected FE, quarterly
	qui spj_manual, yvar(`yvar') xvar(L3.`yvar' i.age i.year i.gender) ///
	eststo(`yvar'_q)
	
	* (4) Bias-corrected FE, monthly
	qui spj_manual, yvar(`yvar') xvar(L.`yvar' i.age i.year i.gender) ///
	eststo(`yvar'_m)
	
	esttab  `yvar'_pool `yvar'_y `yvar'_q `yvar'_m, b(3) se keep(*`yvar'*) ///
	rename(L12.`yvar' L.`yvar' L3.`yvar'  L.`yvar' ) ///
	mtitle("pooled" "bc fe year" "bc fe quart" "bc fe month")	
	
	esttab  `yvar'_pool `yvar'_y `yvar'_q `yvar'_m ///
	using "~/Desktop/table_domain.tex", tex append ///
	b(3) se keep(*`yvar'*) ///
	rename(L12.`yvar' L.`yvar' L3.`yvar'  L.`yvar' ) ///
	mtitle("pooled" "bc fe year" "bc fe quart" "bc fe month")	
	
	}
	
	
*** Appendix table: No balanced sample restriction:
	
	esttab pool_y1_bigsample pool_y2_bigsample pool_q1_bigsample ///
			pool_q2_bigsample pool_m_bigsample , b(3) se keep(*lifesat*) ///
	rename(L3.lifesat L.lifesat L12.lifesat L.lifesat) ///
	mtitle( "yearly" "yearly" "quarterly" "quarterly" "monthly")
	
	esttab pool_y1_bigsample pool_y2_bigsample pool_q1_bigsample ///
			pool_q2_bigsample pool_m_bigsample ///
			using "~/Desktop/table_unbalanced.tex", replace tex ///
			b(3) se keep(*lifesat*) ///
	rename(L3.lifesat L.lifesat L12.lifesat L.lifesat) ///
	mtitle( "yearly" "yearly" "quarterly" "quarterly" "monthly")

*** Appendix table: Descriptive statistics

	estpost sum lifesat age year gender
	est sto descstat
	
	esttab descstat using "~/Desktop/table_descstat.tex", ///
		cells("mean(fmt(2)) sd(fmt(2)) min(fmt(0)) max(fmt(0))" ) ///
		replace tex
