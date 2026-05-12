#!/bin/bash -e
#SBATCH --account=ga03488
#SBATCH --job-name=16S_tents
#SBATCH --time=24:00:00
#SBATCH --mem=120G
#SBATCH --output 16S__%j.out  
#SBATCH --error 16S__%j.err 

##### unload modules ####
module purge

##### Load required ones ####
module load R/4.3.2-foss-2023a

##### Code to run ####
Rscript 16S_full.r