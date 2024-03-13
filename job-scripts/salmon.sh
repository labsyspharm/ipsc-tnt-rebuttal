#!/bin/bash
#SBATCH -J salmon
#SBATCH -c 4                               # Request one core
#SBATCH -t 0-0:20                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=25G                         # Memory total in MiB (for all cores)
#SBATCH -o logs/%x_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e logs/%x_%j.err                 # File to which STDERR will be written, including job ID (%j)

# comm -1 -3 <(ls -1 quants_yamanaka | grep -e .done | sed s/.done//g | sort | uniq) <(find raw/SRP* -name '*_1.fastq.gz' | grep -v SRP223953 | grep -v SRP262162 | sed s/_1.fastq.gz//g | sed 's/raw\/SRP.*\///g' | sort | uniq)
# USAGE EXAMPLE
# sbatch salmon.sh /n/scratch3/users/c/ch305/rna-editing/rna-seq/fastq/Homo_sapiens.GRCh38.gentrome_including_variants_index 962_TCX_2

source ~/miniconda3/etc/profile.d/conda.sh
conda activate salmon

set -eux
INDEX="$1"
OUTPUT_PATH="$2"

# argument 1 contains path to salmon index
# argument 1 contains path to output folder
# 3 onwards contains sample paths
args=("$@")
for i in "${args[@]:2}"
do
dir=${i%/*}
file="${i##*/}"
echo "$SLURM_JOB_ID" > "${OUTPUT_PATH}/${file}.lock"
# --recoverOrphans causes segfaults at the moment
if [ -f "${i}_2.fastq.gz" ]; then
  # Paired-end
  salmon quant -i "$INDEX" -l A \
    --seqBias --gcBias --posBias --softclip -p 4 \
    -1 "${i}_1.fastq.gz" \
    -2 "${i}_2.fastq.gz" \
    -o "${OUTPUT_PATH}/${file}"
else
  # single-end
  salmon quant -i "$INDEX" -l A \
    --seqBias --gcBias --posBias --softclip -p 4 \
    -r "${i}_1.fastq.gz" \
    -o "${OUTPUT_PATH}/${file}"
fi
rm "${OUTPUT_PATH}/${file}.lock"
echo "$SLURM_JOB_ID" > "${OUTPUT_PATH}/${file}.done"
done
echo "done"

