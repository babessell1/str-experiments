#!/usr/bin/bash

#SBATCH --job-name=strling_sum
#SBATCH --partition=standard
#SBATCH --account=remills99
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=36:00:00
#SBATCH --output=logs/strling.log
#SBATCH --error=logs/strling.err

eval "$(conda shell.bash hook)"

# Set up any environment variables or activate a virtual environment if needed
conda activate str-downstream

python run_strling.py