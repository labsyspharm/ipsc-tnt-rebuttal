#!/bin/bash
#SBATCH -J salmon
#SBATCH -c 4                               # Request one core
#SBATCH -t 0-0:20                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=25G                         # Memory total in MiB (for all cores)
#SBATCH -o logs/%x_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e logs/%x_%j.err                 # File to which STDERR will be written, including job ID (%j)

# USAGE EXAMPLE
# sbatch salmon.sh /n/scratch3/users/c/ch305/rna-editing/rna-seq/fastq/Homo_sapiens.GRCh38.gentrome_including_variants_index 962_TCX_2

source ~/miniconda3/etc/profile.d/conda.sh
conda activate salmon

OUTPUT_PATH=quants

set -eux

# argument 1 contains path to salmon index
# 2 onwards contains sample paths
args=("$@")
for i in "${args[@]:1}"
do
dir=${i%/*}
file="${i##*/}"
echo "$SLURM_JOB_ID" > "${OUTPUT_PATH}/${file}.lock"
# --recoverOrphans causes segfaults at the moment
if [ -f "${i}_2.fastq.gz" ]; then
  salmon quant -i "$1" -l A \
    --seqBias --gcBias --posBias --softclip -p 4 \
    -1 "${i}_1.fastq.gz" \
    -2 "${i}_2.fastq.gz" \
    -o "${OUTPUT_PATH}/${file}"
else
  salmon quant -i "$1" -l A \
    --seqBias --gcBias --posBias --softclip -p 4 \
    -r "${i}_1.fastq.gz" \
    -o "${OUTPUT_PATH}/${file}"
fi
rm "${OUTPUT_PATH}/${file}.lock"
echo "$SLURM_JOB_ID" > "${OUTPUT_PATH}/${file}.done"
done
echo "done"

