#!/bin/bash -e
#SBATCH --account=ga03488
#SBATCH --job-name=ITS_tents
#SBATCH --time=24:00:00
#SBATCH --mem=120G
#SBATCH --output ITS__%j.out  
#SBATCH --error ITS__%j.err   

##### unload modules ####
module purge

##### Load required ones ####
module load R

##### Code to run ####
Rscript ITS_full.r