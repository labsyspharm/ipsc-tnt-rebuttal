# pysradb download --out-dir . --srp SRP259918
# pysradb download --out-dir . --srp SRP286549

fastq-dl -a SRP259918 -o SRP259918
fastq-dl -a SRP286549 -o SRP286549

pysradb metadata SRP259918 > SRP259918/meta.tsv
pysradb metadata SRP286549 > SRP286549/meta.tsv

python merge_reads.py meta.tsv
