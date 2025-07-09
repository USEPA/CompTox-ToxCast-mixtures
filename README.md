# CompTox-ToxCast-mixtures
Authors: Broughton RA, Feshuk M, Stanfield I, Isaacs K, Paul Friedman K

## Overview
This repository accompanies pilot work to implement mathematical models of mixtures responses from ToxCast bioactivity data. The code performs the following steps:
1.	Fit the single chemical component concentration-response data using three different techniques: the best fit curve from tcplfit2, a bootstrap resampling approach, and a Bayesian statistical framework for estimation.
2.	Simulate modeled mixture concentration-response from the single component concentration-response curves (obtained through each of the three different fitting methods in Step 1) via the CA, IA, and MP mixture models. The full concentration-response curve is modeled as well as just the concentation at the point of departure (activity concentration at the cutoff, ACC).
3.	Evaluate the abilities of the different mixture models to capture the observed mixture concentration-response for each mixture in the test dataset through a fit metric analysis.
4.	Apply the implemented mixture concentration-response simulation capabilities to perform a case study on bioactivity:exposure ratio (BER) comparisons.

## Usage
The code is split into four sections that can be run in sequence or independently. All steps require tcpl and tcplfit2.
1. Fit single component concentration response with best fit model from tcplfit2, bootstrap resampling, and a Bayesian framework. Simulate the CA, IA, and MP mixture models with each fitting approach. 
  * File: 1_mixtures_analysis.Rmd
2.
3.
4.
5.

A vignette is provided that walks through an example of the full procedure for simulating a mixture model from single chemical components for one example binary mixture.

