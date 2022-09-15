# Chicken
ref=data/galGal3/galGal3_all.fa
./JupiterPlot/jupiter2 name=todChl-to-galGal3 ref=$ref fa=genomes/todChl/todChl.fasta

# Hornbill
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/769/605/GCA_009769605.1_bBucAby1.pri/GCA_009769605.1_bBucAby1.pri_genomic.fna.gz
gunzip GCA_009769605.1_bBucAby1.pri_genomic.fna.gz
ref=GCA_009769605.1_bBucAby1.pri_genomic.fna
./JupiterPlot/jupiter2 name=todChl-to-bucAby ref=$ref fa=genomes/todChl/todChl.fasta

# Bee eater
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/819/595/GCA_009819595.1_bMerNub1.pri/GCA_009819595.1_bMerNub1.pri_genomic.fna.gz
gunzip GCA_009819595.1_bMerNub1.pri_genomic.fna.gz
ref=GCA_009819595.1_bMerNub1.pri_genomic.fna
./JupiterPlot/jupiter2 name=todChl-to-merNub ref=$ref fa=genomes/todChl/todChl.fasta
