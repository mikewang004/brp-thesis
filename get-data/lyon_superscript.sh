#!bin/sh
#
#SBATCH --job-name="JRunAnalysis"
#SBATCH --partition=compute
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G

./download_data.sh
