#!/bin/bash
#SBATCH -J star
#SBATCH -c 8                               # Request one core
#SBATCH -t 0-0:15                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=40G                         # Memory total in MiB (for all cores)
#SBATCH -o logs/%x_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e logs/%x_%j.err                 # File to which STDERR will be written, including job ID (%j)

#
# comm -1 -3 <(ls -1 star | grep -e .done -e .lock | sed s/.done//g | sed s/.lock//g | sort | uniq) <(ls -1 trimmed | grep _1.fastq.gz | sed s/_1.fastq.gz//g | sort | uniq) |
#  parallel sbatch star.sh /n/scratch3/users/c/ch305/indices/hg38_107/star trimmed/{} star_first/*SJ.out.tab

# USAGE EXAMPLE
# sbatch star.sh /n/scratch3/users/c/ch305/indices/hg38_107/star trimmed/SRR123 (tab)

module load gcc/6.2.0 star/2.7.9a

OUTPUT_PATH=star

set -eux

# argument 1 contains path to salmon index
# 2 onwards contains sample paths
args=("$@")
for i in "${args[@]:1}"
do
dir=${i%/*}
file="${i##*/}"
echo "$SLURM_JOB_ID" > "${OUTPUT_PATH}/${file}.lock"
# --seqBias --gcBias --posBias --recoverOrphans --softclip -p 4 \
if [ -f "${i}_2.fastq.gz" ]; then
  STAR --runThreadN 8 \
    --genomeDir "$1" \
    --readFilesIn "${i}_1.fastq.gz" "${i}_2.fastq.gz" \
    --readFilesCommand gzip -cd \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix "${OUTPUT_PATH}/${file}_" \
    --chimOutType WithinBAM
else
  STAR --runThreadN 8 \
    --genomeDir "$1" \
    --readFilesIn "${i}_1.fastq.gz" \
    --readFilesCommand gzip -cd \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix "${OUTPUT_PATH}/${file}_" \
    --chimOutType WithinBAM
fi
rm "${OUTPUT_PATH}/${file}.lock"
echo "$SLURM_JOB_ID" > "${OUTPUT_PATH}/${file}.done"
done
echo "done"

