/*============================================================================== 
   
	Title:          Demonstration of NMF on the Olivetti Faces dataset
	Author(s):      Dr Jonny Batty, BSc(Hons) MBChB(Hons) MPH
	Contact(s):     J.Batty@leeds.ac.uk
	Date updated:   17th January 2024
	License:		MIT
								
	Overview
	--------
	This script demonstrates how NMF can be used in Stata to extract fetures  
	from images of faces. These features can then be used to 'compress' face 
	images, or as part of a facial recognition algorithm. This uses the Olivetti 
	Faces dataset, which contains a set of 400 images taken between April 1992 
	and April 1994 at AT&T Laboratories, Cambridge. There are ten different 
	images of each of 40 distinct subjects, stored as 64 x 64 pixel images. The 
	images have been quantised to 256 grey levels, converted to floating point 
	values in the range [0, 1].
	
==============================================================================*/

clear all


// Install nmf
net install nmf, from("https://raw.githubusercontent.com/jonathanbatty/stata-nmf/main/installation/") replace

// Install dependencies specific to this demonstration 
ssc install heatplot
ssc install palettes
ssc install colrspace
ssc install gtools

// Set seed
set seed 12345

// Change default frame name
frame rename default X

// Load Olivetti Faces dataset
use "https://raw.githubusercontent.com/jonathanbatty/stata-nmf/main/installation/faces.dta"
// Each column of the faces data matrix, W, contains the pixel vaues for an 8 x 8 
// facial image. Applying NMF to these data will give W (the basis matrix) containing
// facial features such as eyes, noses, lips etc., and H (the coeficcient matrix) 
// wite indicates which to what extent each feature is included in each image.

// Define a program that plots a matrix containing pixel data
capture program drop plotimage
program define plotimage
	args matrix_name plot_name y_size x_size
	
	heatplot `matrix_name', colors(HCL reds, gscale) ///
                            plotregion(margin(zero)) ///
			                legend(off) ///
			                xlabel(none, nogrid) xtitle("") xscale(noline) xsize(`x_size') ///
			                ylabel(none, nogrid) ytitle("") yscale(noline) ysize(`y_size') ///
					 	    plotregion(margin(0 0 0 0)) ///
						    name(`plot_name')
end

// Select a single random person to preview
local randomperson = runiformint(1, 400)
display "Random selected ID: " `randomperson'

mkmat v`randomperson', matrix("randomperson")
mata: rmatrix = colshape(st_matrix("randomperson"), 64)
mata: st_matrix("randomperson", rmatrix)

plotimage randomperson selectedperson 8 8

// Preview one image of each person
// Firstly, create an empty matrix to hold the face preview
mata: preview_matrix = J(320, 512, .)
mata: st_matrix("preview_matrix", preview_matrix)

local counter = 0
// Select 4 x 10 = 40 faces to display (one for each person)
forvalues i = 1/5 {
	forvalues j = 1/8 {
		
		// Select face id randomly (apart from first, which has been selected)
		local face = (`counter' * 10) + 1
		local counter = `counter' + 1
		
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
plotimage preview_matrix previewgrid 5 8
						 
// Drop preview matrices
mata: mata clear
matrix drop preview_matrix

// Perform NMF
local rank = 40
nmf v*, k(`rank') epoch(1000) initial(randomu) loss(eu) nograph

// Change to frame containing basis matrix (W) and plot these features
frame change W

mata: feature_matrix = J(320, 512, .)
mata: st_matrix("feature_matrix", feature_matrix)

local counter = 0
forvalues i = 1/5 {
	forvalues j = 1/8 {
		
		local counter = `counter' + 1
	
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

// Plot features
plotimage feature_matrix basismatrix 5 8			 

// Plot features for a single randomly-selected person in frame
frame change H
mkmat H`randomperson', matrix("coefmatrix")
mata: rmatrix = colshape(st_matrix("coefmatrix"), 8)
mata: st_matrix("coefmatrix", rmatrix)

plotimage coefmatrix coefficientmatrix 8 8
				 
// Multiply W and H and recapitulate initial matrix
frame W: mkmat W*, matrix(W)
frame H: mkmat H*, matrix(H)
matrix X_hat = W * H

frame create X_hat
frame X_hat: svmat X_hat, name(x)

// Plot face of randomly selected person
frame change X_hat
mkmat x`randomperson', matrix("reconstructperson")
mata: reconmatrix = colshape(st_matrix("reconstructperson"), 64)
mata: st_matrix("reconstructperson", reconmatrix)

plotimage reconstructperson reconstructedperson	8 8

// Plot gallery of reconstituted faces (same as initial gallery)
mata: reconstituted_matrix = J(320, 512, .)
mata: st_matrix("reconstituted_matrix", reconstituted_matrix)

local counter = 0
forvalues i = 1/5 {
	forvalues j = 1/8 {
		
		// Select face id randomly (apart from first, which has been selected)
		local face = (`counter' * 10) + 1
		local counter = `counter' + 1
		
		// Create 64 x 64 matrix from the randomly selected column (variable)
		mkmat x`face', matrix("facematrix")
		mata: mmatrix = colshape(st_matrix("facematrix"), 64)
		mata: st_matrix("facematrix_updated", mmatrix)
		
		// Calculate which section of preview matrix to replace
		local topleftr = ((`i' - 1) * 64) + 1
		local topleftc = ((`j' - 1) * 64) + 1
		
		// Overwrite preview_matrix at correct co-ords
		matrix reconstituted_matrix[`topleftr', `topleftc'] = facematrix_updated
	}
}

plotimage reconstituted_matrix recongrid 5 8

// This has demonstrated the compression and reconstitution of a dataset 
// containing 64 x 64 x 400 = 1,638,400 datapoints (pixels) to two smaller 
// matrices, of size 4096 x 40 = 163,840 (W) and 40 * 400 = 16,000 (H). In 
// doing this, the algorithm detected common facial features between people.

// End of file


