#!/bin/bash
#SBATCH -J samtools-sort
#SBATCH -c 4                               # Request one core
#SBATCH -t 0-0:10                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=10G                         # Memory total in MiB (for all cores)
#SBATCH -o logs/%x_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e logs/%x_%j.err                 # File to which STDERR will be written, including job ID (%j)


set -eu
module load gcc/9.2.0 samtools/1.15.1

set -x

args=("$@")
for i in "${args[@]}"
do
FILE_PREFIX="${i%.*}"

samtools sort -m 2G -@ 4 -o "${FILE_PREFIX}_sorted.bam" "$i"
samtools index -@ 4 "${FILE_PREFIX}_sorted.bam"

done
echo "done"
