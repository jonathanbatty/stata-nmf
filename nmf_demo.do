clear all
set seed 12345

import delimited "C:\Users\jonat\OneDrive\Desktop\NMF\img.csv"

timer clear
qui {
	nois _dots 0, title(Timing NMF:) reps(100)
	forvalues i = 1/100 {
		timer on 1
		nmf v*, 				///
		    k(20) 				///
			iter(500) 			///
			initial(nndsvd) 	///
			stop(1.0e-4) 		///
			method(cd) 			///
			beta(2) 			///
			nograph 			///
			noframes
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
import delimited "C:\Users\jonat\OneDrive\Desktop\NMF\img.csv"

nmf v*, 				///
	k(20) 				///
	iter(500) 			///
	initial(nndsvd) 	///
	stop(1.0e-4) 		///
	method(cd) 			///
	beta(2) 			///
	nograph noframes	

	
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
A = (-1, 2, -3 \ 4, 0, 6 \ 7, 8, 9)
(A :>= 0) :* A
A= (A :>= 0) :* A
A
end




