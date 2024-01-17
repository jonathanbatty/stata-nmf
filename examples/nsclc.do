/*============================================================================== 
   
	Title:          Demonstration of NMF on a NSCLC gene expression dataset
	Author(s):      Dr Jonny Batty, BSc(Hons) MBChB(Hons) MPH
	Contact(s):     J.Batty@leeds.ac.uk
	Date updated:   17th January 2024
	License:		MIT
								
	Overview
	--------
	This script demonstrates how NMF can be used in Stata to identify 'metagenes'
	among 100 patients with non-small cell lung cancer (NSCLC).
	
==============================================================================*/

clear all

// Install nmf
net install nmf, from("https://raw.githubusercontent.com/jonathanbatty/stata-nmf/main/installation/") replace

// Set seed
set seed 12345

// Load non-small cell lung cancer gene expression dataset
use "https://raw.githubusercontent.com/jonathanbatty/stata-nmf/main/installation/nsclc.dta"
// Each column 

// Add another frame to hold matrix data
frame copy default X
frame X: drop gene_id

// Perform NMF
local rank = 40
nmf v*, k(`rank') epoch(1000) initial(randomu) loss(eu) nograph

// End of file




net install nmf, from("https://raw.githubusercontent.com/jonathanbatty/stata-nmf/main/installation/") replace

nmf p*, k(5) epoch(500)

nmf p*, k(5)

help nmf