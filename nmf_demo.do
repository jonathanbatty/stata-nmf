clear all
set seed 12345

import delimited "img.csv"

timer clear
qui {
	nois _dots 0, title(Timing NMF:) reps(1000)
	forvalues i = 1/1000 {
		timer on 1
		nmf v*, 				///
			k(20) 				///
			iter(500) 			///
			initial(randomu) 	///
			stop(1.0e-4) 		///
			method(mu) 			///
			loss(is) 			///
			nograph noframes	
		timer off 1
		nois _dots `i' 0
	}
}
timer list

matrix W = r(W)
matrix H = r(H)
matrix norms = r(norms)
matrix list norms

// Timings (over 100 iterations):
// Command - nmf v*, k(20) iter(500) initial(nndsvd) stop(1.0e-4) method(cd) beta(2) nograph noframes
// Unoptimised - 3.0363
// Mean function put in if statements - 2.8732
// Variable types set - 2.5845
// 'verbose' option with no printfs -
// Use pointers - 
// 


// Test of variable declarations
clear all
set matastrict on
do nmf.ado
set matastrict off







// Run single test
clear all
set seed 12345

//import delimited "sample_data/img.csv"

import delimited "sample_data/nsclc.csv"
rename v1 gene

// Make a percentage of data missing (MCAR)

// gen origorder = _n
// gen randsort = .
// foreach var of varlist v* {
// 	replace randsort = runiform()
// 	sort randsort
// 	replace `var' = . if runiform() < 0.05
// }
// sort origorder
// drop randsort
// drop origorder

// Run NMF
nmf p*, 				///
	k(15) 				///
	iter(2000) 			///
	initial(randomu) 	///
	stop(1.0e-4) 		///
	method(mu) 			///
	loss(eu) 			///
	//nograph noframes	

	
	
	

matrix W = r(W)
matrix H = r(H)
matrix norms = r(norms)
matrix list norms




// matrix list W
// matrix list H



// matrix A = (1, 2, 3 \ 4, 5, 6 \ 7, 8, 9)
// matrix list A
//
// matrix svd U W _V = A
// matrix V = _V'
//
// matrix list U
// matrix list W
// matrix list V


// clear all
// mata:
// components = 5
// permutation = rowshape(range(1, components, 1), 1)
// permutation
// permutation[3]
// end


clear all
mata:
A = (., 2, 3 \ 4, 0, 6 \ ., 8, 9)
A

mask = (A :!= .)
maskedA = A :* mask

mask
maskedA

_editmissing(A, 0)
A

end




