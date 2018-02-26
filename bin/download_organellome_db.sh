mkdir organellome_db
echo "Downloading mitochondrial genomes"
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.1.genomic.fna.gz
mv mitochondrion.1.1.genomic.fna.gz organellome_db/mito1.fa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.2.1.genomic.fna.gz
mv mitochondrion.2.1.genomic.fna.gz organellome_db/mito2.fa.gz
echo "Downloading chloroplastic genomes"
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.1.1.genomic.fna.gz
mv plastid.1.1.genomic.fna.gz organellome_db/plast1.fa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.2.1.genomic.fna.gz
mv plastid.2.1.genomic.fna.gz organellome_db/plast2.fa.gz
echo "Extracting files"
gunzip organellome_db/mito1.fa.gz
gunzip organellome_db/mito2.fa.gz
gunzip organellome_db/plast1.fa.gz
gunzip organellome_db/plast2.fa.gz
cat organellome_db/*.fa > organellome_db/organellome.fa
rm organellome_db/plast*
rm organellome_db/mito*
