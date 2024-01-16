clear all
cls

// Install dependencies (for this demo; none are required for nmf)
// ssc install heatplot, replace
// ssc install palettes, replace
// ssc install colrspace, replace
// ssc install gtools, replace

// Set seed
set seed 12345

// Load Olivetti Faces dataset
use faces

// Select some random faces to preview

// Firstly, create an empty matrix to hold the face preview
mata: preview_matrix = J(4032, 4032, .)
mata: st_matrix("preview_matrix", preview_matrix)

// Select 8 x 8 = 64 faces to display
local counter = 0
forvalues i = 1/8 {
	forvalues j = 1/8 {
		
		local counter = `counter' + 1
		
		// Select face id randomly 
		local face = runiformint(1, 400)
		
		// Create 64 x 64 matrix from the randomly selected column (variable)
		mkmat v`face', matrix("facematrix")
		mata: mmatrix = colshape(st_matrix("facematrix"), 64)
		mata: st_matrix("facematrix_updated", mmatrix)
		
		// Calculate which section of preview matrix to replace
		local topleftr = ((`i' - 1) * 64) + 1
		local topleftc = ((`j' - 1) * 64) + 1
		
		//di "i = " `i' ", j = " `j' ", topleftr = " `topleftr' ", topleftc = " `topleftc'
 		
		// Overwrite preview_matrix at correct co-ords
		matrix preview_matrix[`topleftr', `topleftc'] = facematrix_updated
	}
}

// Plot image
heatplot preview_matrix, colors(HCL reds, gscale) ///
                         plotregion(margin(zero)) ///
			             legend(off) ///
			             xlabel(none, nogrid) xtitle("") xscale(noline) xsize(8) ///
			             ylabel(none, nogrid) ytitle("") yscale(noline) ysize(8) ///
					 	 aspectratio(1) plotregion(margin(0 0 0 0))
						 
// Drop preview matrices
mata: mata clear
matrix drop preview_matrix

// 