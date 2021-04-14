#!/bin/sh
#
# Reserve 8 CPUs for this job
#$ -pe parallel 6
# Request 8G of RAM
#$ -l h_vmem=8G
# The path used for the standard output stream of the job
#$ -o /ebio/abt5_projects/Sex_biased_genes_Phaeophyceae/code/
# Merge stdout and stderr. The job will create only one output file which
# contains both the real output and the error messages.
#$ -j y
# Use /bin/bash to execute this script
#$ -S /bin/bash
#Run job from current working directory
#$ -cwd
# Send email when the job begins, ends, aborts, or is suspended
#$ -m beas


conda activate geneOntology_venv

emapper.py -m diamond --data_dir /ebio/abt5_projects/Sex_biased_genes_Phaeophyceae/data/Dictyota/eggnog-mapper/data -i /ebio/abt5_projects/Sex_biased_genes_Phaeophyceae/data/Dictyota/Dictyota-dichotoma_MALE_proteins.fa -o emapper-output --output_dir /ebio/abt5_projects/Sex_biased_genes_Phaeophyceae/data/Dictyota/eggnog-mapper
