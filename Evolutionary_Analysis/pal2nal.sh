#!/bin/sh
#
# Reserve 8 CPUs for this job
#$ -pe parallel 4
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



input_dir1="/ebio/abt5_projects/Sex_biased_genes_Phaeophyceae/data/Comparative_analysis/proteomes/Results_Aug12/muscle_output"
input_dir2="/ebio/abt5_projects/Sex_biased_genes_Phaeophyceae/data/Comparative_analysis/proteomes/Results_Aug12/evo_Nucleotides"
output_dir="/ebio/abt5_projects/Sex_biased_genes_Phaeophyceae/data/Comparative_analysis/proteomes/Results_Aug12/pal2nal_output"
conda activate seqtk-1.3

for f in $(ls $input_dir2/*.fa | cut -d/ -f10)

do pal2nal.pl -output fasta ${input_dir1}/$f.output ${input_dir2}/$f > $output_dir/aligned_nucl_$f

done
