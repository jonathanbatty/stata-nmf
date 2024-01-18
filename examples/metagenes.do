/*============================================================================== 
   
	Title:          Demonstration of NMF on a NSCLC gene expression dataset
	Author(s):      Dr Jonny Batty, BSc(Hons) MBChB(Hons) MPH
	Contact(s):     J.Batty@leeds.ac.uk
	Date updated:   17th January 2024
	License:		MIT
								
	Overview
	--------
	This script demonstrates how NMF can be used in Stata to identify 'metagenes'
	in microarray gene expression data, from 100 samples of non-small cell lung 
	cancer (NSCLC) tissue. This is a subset of open data from: Botling J, et al. 
	Biomarker discovery in non-small cell lung cancer: integrating gene expression 
	profiling, meta-analysis, and tissue microarray validation. Clin Cancer Res. 
	2013;19(1):194–204.
	
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

// Add another frame to hold matrix data
frame copy default A
frame A: drop gene_id

// Perform NMF
local rank = 5
frame A: nmf p*, k(`rank') epoch(1000) initial(randomu) loss(eu)
// When factorised, the columns of W can be interpreted as 'metagenes' and H contains
// the degree of expression of each metagene for a given individual.

// Display available frames
frames dir

// View first 20 rows of frame W (containing the gene membership of each 'metagene')
frame W: list in 1/20

// View first 5 columns of frame (containing expression profiles of each 
// metagene for the first 5 patient samples.
frame H: list H1 - H5

// Therefore, the dimensionality has been reduced from 200 genes to 5 'metagenes'.
// Further work should focus on (i) how to select the rank of the factorisation 
// (i.e. the number of metagenes), and (ii) how expression of the identified 
// 'metagenes' may influence e.g. survival for patients in NSCLC.

// End of file


