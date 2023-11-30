#!/bin/bash
#SBATCH -J sra-to-fastq
#SBATCH -c 4                               # Request one core
#SBATCH -t 0-1:30                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=20G                         # Memory total in MiB (for all cores)
#SBATCH -o sra-to-fastq_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e sra-to-fastq_%j.err                 # File to which STDERR will be written, including job ID (%j)

# USAGE:
# comm -1 -3 <(find . -name '*.done' | sed s/.done//g | sort | uniq) <(find . -name '*.sra' | sed s/.sra//g | sort | uniq) |
# xargs echo sbatch sra-to-fastq.sh

module load gcc/9.2.0
module load sratoolkit/3.0.2
module load pigz/2.3.4

set -eux

for i in $@
do
dir=${i%/*}
file="${i##*/}"
echo "$SLURM_JOB_ID" > "${i}.lock"
fasterq-dump -3 -x -e 4 -m 10G -O "$dir" -o "$file" "${i}.sra"
if [ -f "${i}" ]; then
    mv "${i}" "${i}_1.fastq"
fi
if [ -f "${i}_1.fastq" ]; then
    pigz -p 4 "${i}_1.fastq"
fi
if [ -f "${i}_2.fastq" ]; then
    pigz -p 4 "${i}_2.fastq"
fi
rm "${i}.lock"
echo "$SLURM_JOB_ID" > "${i}.done"
done
echo "done"
