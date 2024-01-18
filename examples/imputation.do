/*============================================================================== 
   
	Title:          Demonstration of NMF for use in imputation
	Author(s):      Dr Jonny Batty, BSc(Hons) MBChB(Hons) MPH
	Contact(s):     J.Batty@leeds.ac.uk
	Date updated:   17th January 2024
	License:		MIT
								
	Overview
	--------
	This script demonstrates how NMF can be used in Stata to impute missing data.
	This will be demonstrated using a dataset that contains a subset of open
	microarray data, from: Botling J, et al. Biomarker discovery in non-small 
	cell lung cancer: integrating gene expression profiling, meta-analysis, and 
	tissue microarray validation. Clin Cancer Res. 2013;19(1):194–204. Missing 
	data in this dataset will be simulated. NMF will be used to impute these
	missing data in order to recapitulate values similar to those in the original 
	dataset.
	
==============================================================================*/

clear all

// Install nmf
net install nmf, from("https://raw.githubusercontent.com/jonathanbatty/stata-nmf/main/installation/") replace

// Set seed
set seed 12345

// Load non-small cell lung cancer gene expression dataset
use "https://raw.githubusercontent.com/jonathanbatty/stata-nmf/main/installation/nsclc.dta"
// Each row contains expression data for a given gene; each column contains data
// pertaining to a patient. 

// Create a new frame to store data matrix
frame copy default A
frame change A
drop gene_id


// Simulate missing data across all columns of microarray data:
foreach var of varlist p* {
	
	// Randomly select a proportion of rows to set to missing, from 0 to 0.3 
	// (i.e. up to 30% missingness of each column)
	local missing_proportion = runiform(0, 0.3)
	
	forvalues i = 1 / `= _N' {
		if runiform() <= `missing_proportion' {
			qui replace `var' = . in `i'
		}
	}
}

// View data to demonstrate presence of missing values
browse

// Run NMF to factorise gene expression matrix
nmf p*, k(20) epoch(1000) initial(randomu) loss(eu) stop(1.0e-5) nograph

// Multiply W and H to recapitulate initial data matrix
frame W: mkmat W*, matrix(W)
frame H: mkmat H*, matrix(H)
matrix A_hat = W * H

frame create A_hat
frame A_hat: svmat A_hat, name(p)

// Show columns from each frame for comparison
display "Original (complete) dataset:"
frame default: list p001 - p005 in 1/20

display "Dataset with simulated missing values:"
frame A: list p001 - p005 in 1/20

display "Dataset with missing values imputed by NMF:"
frame A_hat: list p1 - p5 in 1/20

// Note: this approach requires that the data are not missing completely at 
// random (MCAR): correlates of missingness will be required in order for NMF to 
// successfully predict missing values.

// NMF as a tool in multiple imputation:
// NMF could be repeated multiple times. After each iteration, the imputed
// values could be combined with the nonmissing values to give one potential 
// complete dataset. These datasets could be analysed using Rubin's rules in 
// line with best practices in multiple imputation.

// End of file