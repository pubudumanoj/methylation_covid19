#!/bin/bash
#SBATCH --account=rrg-bourqueg-ad
#SBATCH -N 1 -n 10
#SBATCH --mem=250GB
#SBATCH -t 16:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pubudu.nawarathna@mail.mcgill.ca
#SBATCH -J score.job.0
#SBATCH -o covid_methylation_all%j.out
#SBATCH -e covid_methylation_all%j.err
set -eu -o pipefail

module load mugqic/R_Bioconductor/3.6.0_3.9

Rscript day0.2021.oct.R
