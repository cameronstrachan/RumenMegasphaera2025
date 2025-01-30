# RumenMegasphaera2025
This repository is for the analysis code from the Strachan et al. 2025 paper "Distinct lactate utilization strategies drive niche differentiation between two co-existing Megasphaera species in the rumen microbiome"

## Statistical analysis and figure generation
The code for the analysis of the various output tables can be found in the folder 'figures'. The analysis is organized by the figures that it was used to generate. The input data for these analysis are either found in the output folders within 'data'. All this code can be run from this repository to reproduce any of the figures.

## Bioinformatic processing
The code used to process sequence data to various output tables (ex. count tables) is found in the folder 'processing'. This is where the command line settings for various bioinformatic tools can be found. Most time, commands are run within pipelines that are written in python or bash. Due to the file sizes, the raw input data needs to be downloaded and then the locations and names of the input data would need to be changed.
