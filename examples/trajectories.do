/*============================================================================== 
   
	Title:          Demonstration of NMF for use in identifying teporal  
					trajectories using clinical data
	Author(s):      Dr Jonny Batty, BSc(Hons) MBChB(Hons) MPH
	Contact(s):     J.Batty@leeds.ac.uk
	Date updated:   17th January 2024
	License:		MIT
								
	Overview
	--------
	This script demonstrates how NMF can be used in Stata to derive trajectories
	of disease accumulation over time. Given a person-time x disease state matrix,
	NMF will decompose this into a matrix that represents the expression of a 
	disease trajectory over time (W), and another that represents the membership of 
	each disease state in a given trajectory (H). Note that the data used in 
	this example (disease_matrix.dta) are entirely synthetic and illustrative.
	The presence or absence of each disease state at a given time point for a 
	given patient was determined by a univariate probability distribution.
	
==============================================================================*/

clear all

// Install nmf
net install nmf, from("https://raw.githubusercontent.com/jonathanbatty/stata-nmf/main/installation/") replace

// Set seed
set seed 12345

// Load synthetic disease matrix dataset
use "https://raw.githubusercontent.com/jonathanbatty/stata-nmf/main/installation/disease_matrix.dta"
// This dataset comprises

// Run NMF to factorise dataset
// As this dataset comprises binary data, the Kullback-Leibler divergence is likely
// to be the correct loss metric.
nmf disease_*, k(5) epoch(500) initial(randomu) loss(kl) nograph

// Display available frames
frames dir

// Link original data frame with frame W to get patient ID and time variables in
// frame W
gen obs_no = _n
frame W: gen obs_no = _n
frame change W
frlink 1:1 obs_no, frame(default)
frget patient_id time, from(default)
drop obs_no default
order patient_id time W*
frame change default

// View person-time x trajectory matrix (W) for the first patient
frame W: list in 1/10, separator(0)

// View trajectory x disease membership matrix (H) for the first 8 disease states
frame H: list H1 - H8, separator(0)

// End of file


