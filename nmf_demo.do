// Toy example

clear all
set seed 12345

// set obs 10
// forvalues i = 1/5 {
// 	gen v`i' = runiformint(1, 10)
// }

import delimited "M:\Code\Stata\NMF\Resources\image_test\image_values.csv", clear
// timer on 1
nmf v*, k(40) iter(500) initial(random) method(mu) beta(2) stop(1.0e-4) // nograph noframes
// timer off 1
// timer list

matrix W = r(W)
matrix H = r(H)
matrix norms = r(norms)



// matrix list W
// matrix list H
// matrix list norms







clear all
set seed 12345
import delimited "M:\Code\Stata\NMF\Resources\image_test\image_values.csv", clear

qui {
	forvalues i = 2/200 {
		nmf v*, k(`i') iter(1000) initial(random) beta(2) stop(1.0e-4) nograph // noframes
		
		matrix norms = r(norms)
		
		local row `= rowsof(norms)'

		local norm = norms[`row', 1]

		noi display "`norm'"
		
	}
}










// clear all
// set seed 123
// // matrix M = (-4, -3, -2\-1, 0, 1\2, 3, 4)
// // matrix list M
//
//
// mata:
// M = (-1, 2, -3\ -3, -4, 0\ 0, 7, 8)
// M
//
// (M :== 0)
//
// M = (M :== 0) :* epsilon(1)
//
// M
//
// end




clear all
mata:
A = (8, 1, 6 \ 3, 5, 7 \ 4, 9, 2)
B = (15 \ 15 \ 15)
A
B
lusolve(A, B)

C = (16,2,3,13\5,11,10,8\9,7,6,12\4,14,15,1)
D = (34 \ 34 \ 34 \ 34)
C
D
lusolve(C, D)

end








clear all
use "C:\Users\jonat\Desktop\MB Frailty-Volume Outcome PCI\ACS.dta", clear
keep i10_dx*

quietly {
	// Generate 3-character ICD 10 dummy variables for each diagnostic code
	forvalues x = 1/40 {
		noisily display "`x'"
		gen i10_3_dx`x' = substr(i10_dx`x', 1, 3)
		drop i10_dx`x'
		// get levels of variable
		levelsof i10_3_dx`x'
		foreach lev in `r(levels)'{
			// check if dummy variable name already created; if not, create it
			capture confirm variable icd_dx_`lev'
			if !_rc == 0 {
				gen icd_dx_`lev' = i10_3_dx`x' == "`lev'"
			}
			replace icd_dx_`lev' = 1 if i10_3_dx`x' == "`lev'"
		}
	}
}
drop i10_3_dx*
drop icd_dx_inc icd_dx_inv

compress

nmf icd_dx_*, k(10) iter(5) initial(random)
