cd ~/uce-alcedinidae

# head data/galGal3/galGal3_all.fa

# Convert psl to sam-
# https://github.com/lh3/samtools-legacy/blob/master/misc/psl2sam.pl
# psl2sam.pl todChl.psl > todChl-to-galGal3.sam

bedtools makewindows -i src -g cnee/galGal3.chr_length.txt -w 1000000 > cnee/test.bed
head cnee/test.bed

liftOver -multiple cnee/test.bed cnee/galGal3-to-todChl.over.chain todChlJupiter-agp.bed unMapped
head todChlJupiter-agp.bed

# ./JupiterPlot/jupiter name='todChl-to-galGal3' ref=data/galGal3/galGal3_all.fa fa=genomes/todChl/todChl.fasta


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
