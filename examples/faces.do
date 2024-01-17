clear all
cls

// Install nmf
net install nmf, from("https://raw.githubusercontent.com/jonathanbatty/stata-nmf/main/installation/") replace

// Install demo-specific dependencies (nmf has no dependencies)
ssc install heatplot
ssc install palettes
ssc install colrspace
ssc install gtools

// Set seed
set seed 12345

// Load Olivetti Faces dataset
use faces
// Each column of the faces data matrix, W, contains the pixel vaues for an 8 x 8 
// facial image. Applying NMF to these data will give W (the basis matrix) containing
// facial features such as eyes, noses, lips etc., and H (the coeficcient matrix) 
// wite indicates which to what extent each feature is included in each image.

// Select some random faces to preview

// Firstly, create an empty matrix to hold the face preview
mata: preview_matrix = J(512, 512, .)
mata: st_matrix("preview_matrix", preview_matrix)

// Select 8 x 8 = 64 faces to display
forvalues i = 1/8 {
	forvalues j = 1/8 {
		
		// Select face id randomly 
		local face = runiformint(1, 400)
		
		// Create 64 x 64 matrix from the randomly selected column (variable)
		mkmat v`face', matrix("facematrix")
		mata: mmatrix = colshape(st_matrix("facematrix"), 64)
		mata: st_matrix("facematrix_updated", mmatrix)
		
		// Calculate which section of preview matrix to replace
		local topleftr = ((`i' - 1) * 64) + 1
		local topleftc = ((`j' - 1) * 64) + 1
		
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

// Perform NMF
nmf v*, k(10) epoch(1000) initial(randomu) loss(eu) nograph

// Change to frame containing basis matrix (W) and plot these features
frame change W

mata: feature_matrix = J(448, 448, .)
mata: st_matrix("feature_matrix", feature_matrix)

forvalues i = 1/7 {
	forvalues j = 1/7 {
		
		local counter = ((`i' - 1) * 7) + `j'
	
		if `counter' <= 10 {
			
			// Create 64 x 64 matrix from the randomly selected column (basis)
			mkmat W`counter', matrix("feature")
			mata: fmatrix = colshape(st_matrix("feature"), 64)
			mata: st_matrix("feature_updated", fmatrix)

			// Calculate which section of preview matrix to replace
			local topleftr = ((`i' - 1) * 64) + 1
			local topleftc = ((`j' - 1) * 64) + 1
			
			// Overwrite preview_matrix at correct co-ords
			matrix feature_matrix[`topleftr', `topleftc'] = feature_updated	
		} 
	}
}

// Plot features
heatplot feature_matrix, colors(HCL reds, gscale) ///
                         plotregion(margin(zero)) ///
			             legend(off) ///
			             xlabel(none, nogrid) xtitle("") xscale(noline) xsize(8) ///
			             ylabel(none, nogrid) ytitle("") yscale(noline) ysize(8) ///
					 	 aspectratio(1) plotregion(margin(0 0 0 0))
						 

// Plot features for a single person in frame HCL
frame change H

local random = runiformint(1, 400)


// Multiply W and H and recapitulate initial face
Save result to new frame (X_hat)


