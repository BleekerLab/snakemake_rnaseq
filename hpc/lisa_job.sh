#!/bin/bash
#SBATCH -t 02:00:00 
#SBATCH -p normal 
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=m.galland@uva.nl
#SBATCH --nodes 1

conda activate rnaseq

snakemake --cores 20