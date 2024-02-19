#!/bin/bash
#SBATCH -J bam_to_fastq
#SBATCH -c 8                               # Request one core
#SBATCH -t 0-0:25                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=40G                         # Memory total in MiB (for all cores)
#SBATCH -o logs/%x_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e logs/%x_%j.err                 # File to which STDERR will be written, including job ID (%j)

# USAGE:
# comm -1 -3 <(ls -1 fastq | grep _1.fastq.gz$ | sed s/_1.fastq.gz//g | sort | uniq) <(ls -1 raw | grep .bam$ | grep -v sorted | sed s/.bam//g | sort | uniq) |
#  parallel sbatch bam_to_fastq.sh {}

module load gcc/9.2.0
module load samtools/1.15.1
module load pigz/2.3.4

set -eux
OUTPUT_PATH="$1"

# argument 1 contains path to output folder
# 2 onwards contains sample paths
args=("$@")
for i in "${args[@]:1}"
do
# Remove extension from input i
path_no_ext=${i%.*}
dir=${i%/*}
file="${i##*/}"
file_no_ext=${file%.*}
output_prefix="${OUTPUT_PATH}/${file_no_ext}"
echo "$SLURM_JOB_ID" > "${output_prefix}.lock"
sorted_bam_path="${output_prefix}_sorted.bam"
if [ ! -f "$sorted_bam_path" ]; then
  samtools sort -n -@ 8 -m 4G -o "$sorted_bam_path" "$i"
fi
if [ ! -f "${output_prefix}_1.fastq.gz" ]; then
  samtools fastq -1 >(pigz -p 3 -c > "${output_prefix}_1.fastq.gz") \
    -2 >(pigz -p 3 -c > "${output_prefix}_2.fastq.gz") \
    -0 >(pigz -p 3 -c > "${output_prefix}_unpaired.fastq.gz") \
    -@ 4 "$sorted_bam_path"
  rm "$sorted_bam_path"
fi
rm "${output_prefix}.lock"
echo "$SLURM_JOB_ID" > "${output_prefix}.done"
done

echo "done"
