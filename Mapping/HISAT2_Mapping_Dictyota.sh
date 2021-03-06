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

conda activate map_venv
input="/ebio/abt5_projects/Sex_biased_genes_Phaeophyceae/data/Dictyota/catTrimmomatic_output"
ID=${1}  #will be replaced with sampleID
FORWARD=$input/${ID}_forward_paired.fq.gz
REVERSE=$input/${ID}_reverse_paired.fq.gz
UNPAIR1=$input/${ID}_forward_unpaired.fq.gz
UNPAIR2=$input/${ID}_reverse_unpaired.fq.gz
OUT=${ID}_sort.bam

hisat2 -p 6 -x /ebio/abt5_projects/Sex_biased_genes_Phaeophyceae/data/Dictyota/Ddi_MALE -1 $FORWARD -2 $REVERSE -U $UNPAIR1,$UNPAIR2 | samtools view -@ 6 -bh | samtools sort -@ 6 >$OUT
samtools index -b $OUT $OUT.bai

