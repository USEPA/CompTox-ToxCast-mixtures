# CompTox-ToxCast-mixtures
This repository contains the code needed to produce the results described in Broughton et al., "Development of Mathematical New Approach Methods to Assess Chemical Mixtures." In this work, 3 models of chemical mixture activity (concentration addition, independent action, and considering the mixture as the activity of the single most potent mixture constituent) are evaluated using a test dataset with the goal of informing future simulation of mixtures using ToxCast data to prioritize mixtures for testing.
# How to use this repository
All inputs are provided in the "input" folder.
All output data are provided in the "output" folder.
Mixtures models functions are provided in the "scripts" folder.
Analysis of the mixtures test dataset is performed in 1_mixtures_analysis.Rmd.
Generation of results from this analysis is performed in 2_results_analysis.Rmd.
Figures for Broughton et al. are generated in 3_mixtures_figures.Rmd.
The analysis for the bioactivity:exposure ratio case study, using mixtures modeling to prioritize mixtures for screening, is performed in 4_mixtures_ber_case_study.Rmd.

