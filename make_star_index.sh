# https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

wget ftp://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

gzip -cd Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > Homo_sapiens.GRCh38.dna.primary_assembly.fa
gzip -cd Homo_sapiens.GRCh38.110.gtf.gz > Homo_sapiens.GRCh38.110.gtf

# Include whole Sendai genome
cat Homo_sapiens.GRCh38.dna.primary_assembly.fa ../raw/sequences/pSeV-idSOX2-Pmut.fa > \
  Homo_sapiens.GRCh38.dna.primary_assembly_with_sendai.fa

STAR --runThreadN 12 \
  --runMode genomeGenerate \
  --genomeDir Homo_sapiens.GRCh38.dna.primary_assembly_with_sendai_star \
  --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly_with_sendai.fa \
  --sjdbGTFfile Homo_sapiens.GRCh38.110.gtf \
  --sjdbOverhang 100

