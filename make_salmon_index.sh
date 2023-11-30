# https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/

wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

wget ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Salmon authors recommend including the whole genome as decoy sequences, which
# are taken into account for finding alignments but are not actually quantified

grep "^>" <(gunzip -c Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz) | cut -d " " -f 1 > decoys.txt

sed -i.bak -e 's/>//g' decoys.txt

# Include viral genes and concatenate cDNA and whole genome

cat Homo_sapiens.GRCh38.cdna.all.fa.gz \
  <(gzip "../raw/sequences/Sendai genes.TXT") \
  Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > \
  Homo_sapiens.GRCh38.gentrome_including_sendai.fa.gz



salmon index -t Homo_sapiens.GRCh38.gentrome_including_sendai.fa.gz \
  -d decoys.txt -p 12 \
  -i Homo_sapiens.GRCh38.gentrome_including_sendai_index

