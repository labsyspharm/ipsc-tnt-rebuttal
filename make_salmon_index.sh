# https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/

wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Salmon authors recommend including the whole genome as decoy sequences, which
# are taken into account for finding alignments but are not actually quantified

grep "^>" <(gunzip -c Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz) | cut -d " " -f 1 > decoys.txt

sed -i.bak -e 's/>//g' decoys.txt

# Include viral genes and concatenate cDNA and whole genome
cat Homo_sapiens.GRCh38.cdna.all.fa.gz \
  <(gzip -c "../raw/sequences/Sendai genes.TXT") \
  Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > \
  Homo_sapiens.GRCh38.gentrome_including_sendai.fa.gz

salmon index -t Homo_sapiens.GRCh38.gentrome_including_sendai.fa.gz \
  -d decoys.txt -p 12 \
  -i Homo_sapiens.GRCh38.gentrome_including_sendai_index

# Include variants of Yamanaka transcripts including and exluding UTR
# Filter out the UTR containing variant from the Ensembl cDNA
cat <(gzip -c "../raw/sequences/yamanaka_utr_variants.fa") \
  <(gzip -cd Homo_sapiens.GRCh38.cdna.all.fa.gz | python ../filter_fasta.py - - 'ENST00000259915.13|ENST00000374672.5|ENST00000621592.8|ENST00000325404.3' | gzip -c) \
  Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > \
  Homo_sapiens.GRCh38.yamanaka_utr_variants.fa.gz

salmon index -t Homo_sapiens.GRCh38.yamanaka_utr_variants.fa.gz \
  -d decoys.txt -p 12 \
  -i Homo_sapiens.GRCh38.yamanaka_utr_variants_index
