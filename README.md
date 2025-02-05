# RumenMegasphaera2025
This repository is for the analysis code from the Strachan et al. 2025 paper "Distinct lactate utilization strategies drive niche differentiation between two co-existing Megasphaera species in the rumen microbiome"

## Statistical analysis and figure generation
The code to reproduce the analysis of the various raw data or pipeline output tables can be found in the folder 'figures'. The input data for these analyses are either found in the output folders within 'data', with the exception of the parameters for the ECM modeling, which is found in 'ECM_models'. All this code can be run from this repository to reproduce any of the figures (after ensuring that any packages being loaded are installed). 

## Bioinformatic processing
The code used to process sequence data to create bins and various output tables (ex. count tables) is found in the folder 'processing'. This is where the command line settings for various bioinformatic tools can be found. The commands are run within pipelines that are written in python or bash. 
