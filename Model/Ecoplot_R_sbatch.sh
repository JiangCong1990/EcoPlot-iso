#!/bin/bash

#SBATCH --job-name=EPr_forA       # Specify job name
#SBATCH --qos=short               # Quality of service
#SBATCH --account=vewa            # Charge resources to the 'vewa' account
#SBATCH --ntasks=1               # Specify max number of tasks
#SBATCH --cpus-per-task=20        # Request 16 CPUs per task
#SBATCH --mem-per-cpu=6g           # Memory request
#SBATCH --time=10:00:00            # Set time limit to 3 hours
#SBATCH --partition=computehm      # Specify partition for compute nodes
#SBATCH --output=my_job.o%j       # Standard output file (%j will be replaced by job ID)
#SBATCH --error=my_job.e%j        # Standard error file (%j will be replaced by job ID)


# Load the necessary module for Conda

source  ~/.bashrc

module load miniconda/2024.05

#  Activate the conda environment with R

conda activate my_r_env

# Execute the R script

Rscript Script_SWBiso_Forest.R

Rscript Script_SWBiso_Scenario.R

# Job completion
exit

