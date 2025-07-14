# CompTox-ToxCast-mixtures
Development of Mathematical New Approach Methods to Assess Chemical Mixtures
Authors: Broughton RA, Feshuk M, Stanfield Z, Isaacs KK, Paul Friedman K

## Overview
This repository accompanies pilot work to implement mathematical models of mixtures responses from ToxCast bioactivity data presented in the manuscript _Development of Mathematical New Approach Methods to Assess Chemical Mixtures_, submitted for peer review. The code performs the following steps:
1.	Fit the single chemical component concentration-response data using three different techniques: the best fit curve from _tcplfit2_, a bootstrap resampling approach, and a Bayesian statistical framework for estimation.
2.	Simulate modeled mixture concentration-response from the single component concentration-response curves (obtained through each of the three different fitting methods in Step 1) via the Concentration Addition (CA), Independent Action (IA), and a model that treats the mixture as the most potent single chemical component (MP). The full concentration-response curve is modeled as well as just the concentation at the point of departure (activity concentration at the cutoff, ACC).
3.	Evaluate the abilities of the different mixture models to capture the observed mixture concentration-response for each mixture in the test dataset through a fit metric analysis.
4.	Apply the implemented mixture concentration-response simulation capabilities to perform a case study on bioactivity:exposure ratio (BER) comparisons.

## Usage
The code is split into four sections that can be run in sequence or independently. All steps require tcpl and tcplfit2.
1. Fit single component and observed mixture concentration-response data with best fit model from _tcplfit2_, bootstrap resampling, and a Bayesian framework. Simulate the CA, IA, and MP mixture models from single component concentration-response obtained with each fitting approach. Single component concentration-reponse is obtained from tested chemicals screened in this work and from existing legacy data within ToxCast.
    * File: [1_mixtures_analysis.Rmd](1_mixtures_analysis.Rmd)
2. Perform analysis of results for manuscript. Filter mixture and single component activity rates based on hitcalls. Assess qulaity-of-fit of simulated mixtures from the mathematical models compared to the observed mixture by evaluating the accuracy of ACC predictions, conservative ACC predictive ability, and full concentration-response curves.
    * File: [2_results_analysis.Rmd](2_results_analysis.Rmd)
3. Generate figures for manuscript, including plotting all observed mixtures and the CA, IA, and MP model simulations. 
    * File: [3_mixture_figures.Rmd](3_mixture_figures.Rmd)
4. Conduct bioactivty:exposure ratio (BER) case study to compare simulated mixtures point of departures to [CDC NHANES](https://www.cdc.gov/nchs/nhanes/index.html) human blood serum exposure concentrations. 
    * File: [4_mixtures_ber_case_study.Rmd](4_mixtures_ber_case_study.Rmd)

A vignette is provided that walks through an example of the full procedure for simulating a mixture model from single chemical components for one example binary mixture.

![Graphical Abstract.](/figures/GraphicalAbstract.png)
