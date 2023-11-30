#!/bin/bash
#SBATCH -J fastp
#SBATCH -c 4                             # Request one core
#SBATCH -t 0-0:20                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=2G                         # Memory total in MiB (for all cores)
#SBATCH -o logs/%x_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e logs/%x_%j.err                 # File to which STDERR will be written, including job ID (%j)

source ~/miniconda3/etc/profile.d/conda.sh
conda activate fastp

set -eux

OUTPUT_PATH=trimmed

for i in $@
do
dir=${i%/*}
file="${i##*/}"
echo "$SLURM_JOB_ID" > "${OUTPUT_PATH}/${file}.lock"
if [ -f "${i}_2.fastq.gz" ]; then
  fastp -i "${i}_1.fastq.gz" -I "${i}_2.fastq.gz" \
    -o "${OUTPUT_PATH}/${file}_1.fastq.gz" -O "${OUTPUT_PATH}/${file}_2.fastq.gz" \
    --unpaired1 "${OUTPUT_PATH}/${file}_1.fastq.gz" --unpaired2 "${OUTPUT_PATH}/${file}_2.fastq.gz" \
    --cut_right \
    --detect_adapter_for_pe \
    -w 4 \
    --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --json "${OUTPUT_PATH}/${file}_report.json" --html "${OUTPUT_PATH}/${file}_report.html" \
    --overrepresentation_analysis \
    --report_title "$file"
else
  fastp -i "${i}_1.fastq.gz" \
    -o "${OUTPUT_PATH}/${file}_1.fastq.gz" \
    --cut_right \
    -w 4 \
    --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --json "${OUTPUT_PATH}/${file}_report.json" --html "${OUTPUT_PATH}/${file}_report.html" \
    --overrepresentation_analysis \
    --report_title "$file"
fi
rm "${OUTPUT_PATH}/${file}.lock"
echo "$SLURM_JOB_ID" > "${OUTPUT_PATH}/${file}.done"
done
echo "done"



