# Targeting-Social-Safety-Net-Programs-on-Human-Capabilities
Depends on package CBART, which can be installed from binary with code:

library(devtools);
install_github("LendieFollett/CBART");
library(CBART)

## General

plotting_functions.R
  - contains functions to create plots for both the simulation experiments and Indonesian data analysis
  - is sourced in targeting_simulation_study_results.R and 

functions.R
  - contains functions to create plots for both the simulation experiments and Indonesian data analysis
  - sourced in alatas_analysis.R
  
## Simulation experiments:

one_x_sim.R
  - sources plotting_functions.R

targeting_simulation_study.R
  - runs multivariate simulation study
  - saves plots in Simulation Plots subfolder
  - saves RDS output of predictions, summaries
  
targeting_simulation_study_results.R
  - creates additional plots, tables from above results


## Indonesian data analysis:
alatas_analysis.R
  - Fits OLR, OLR (AE), RF, BART, CBART to all response variables
  - Saves output, makes plots
