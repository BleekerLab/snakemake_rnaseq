#!/bin/bash
#
#SBATCH --job-name=snakemake_rnaseq     # job name
#SBATCH  --time=24:00:00                # mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-type=END,FAIL            # When to send emails
#SBATCH --mail-user=m.galland@uva.nl    # email address
#SBATCH --mem-per-cpu=8G                # RAM requested per job (in Snakemake, one rule = one job)
#SBATCH --output=parallel_%j.log        # standard output and error log (%j substitutes the JOB ID)


#SBATCH --nodes=1                       # Run all processes on a single node	
#SBATCH --cpus-per-task=30              # Number of CPUs per task

source activate rnaseq

srun snakemake -j $SLURM_CPUS_PER_TASK