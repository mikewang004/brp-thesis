#!/bin/sh
#
#SBATCH --job-name="JRunAnalysis"
#SBATCH --partition=htc
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G

./download_data.sh
