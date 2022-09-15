#------------------------------------------------------------------------------
# Annotate with gemoma
#------------------------------------------------------------------------------

# Get genomes/annotations-

mkdir data; cd data

# calAnn
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/957/555/GCA_003957555.1_bCalAnn1_v1.p/GCA_003957555.1_bCalAnn1_v1.p_genomic.fna.gz

# Budgie
wget ftp://ftp.ensembl.org/pub/rapid-release/gff3/melopsittacus_undulatus//Melopsittacus_undulatus.bMelUnd1.mat.Z.101.gff3.gz
wget ftp://ftp.ensembl.org/pub/rapid-release/fasta/melopsittacus_undulatus/dna//Melopsittacus_undulatus.bMelUnd1.mat.Z.dna.toplevel.fa.gz

# aquChr
wget ftp://ftp.ensembl.org/pub/release-100/fasta/aquila_chrysaetos_chrysaetos/dna/Aquila_chrysaetos_chrysaetos.bAquChr1.2.dna_sm.toplevel.fa.gz

#chicken
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.gff.gz

# round 2 = repeat-masked genome
cd ..
mkdir annotations; cd annotations

# Reference genome fast file
infa='../data/todChl.scaffolds.full_mask.fa'

# Hard-mask the reference genome to work with gemoma
sed -e '/^>/! s/[[:lower:]]/N/g' $infa > masked.fasta

# Do the annotation
java -jar -Xms80G -XX:ParallelGCThreads=48 ../GeMoMa-1.7.1.jar CLI GeMoMaPipeline \
	restart=true \
	threads=48 \
	GeMoMa.Score=ReAlign \
	AnnotationFinalizer.r=NO \
	o=true \
	p=false \
	tblastn=false \
	t=masked.fasta \
	s=own \
	i=melUnd \
	g=data/Melopsittacus_undulatus.bMelUnd1.mat.Z.dna.toplevel.fa \
	a=data/Melopsittacus_undulatus.bMelUnd1.mat.Z.101.gff3 \
	s=own \
	i=galGal \
	g=data/GCF_000002315.6_GRCg6a_genomic.fna \
	a=data/GCF_000002315.6_GRCg6a_genomic.gff \
	s=own \
	i=aquChr \
	g=data/Aquila_chrysaetos_chrysaetos.bAquChr1.2.dna_sm.toplevel.fa \
	a=data/Aquila_chrysaetos_chrysaetos.bAquChr1.2.101.gff3 \
	s=own \
	i=calAnn \
	g=data/GCF_003957555.1_bCalAnn1_v1.p_genomic.fna \
	a=data/GCF_003957555.1_bCalAnn1_v1.p_genomic.gff

# --mask-lower-case

gff=final_annotation.gff

# Get number of genes, lengths-
cat $gff | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 - $4) } END { print NR, sum / NR }'
# 16539 genes, 17744.6 bp avg.

# Run gemoma annotation filter (GAF)
# Keep only best (top 1, m=1) matches-
java -jar ../GeMoMa-1.7.1.jar CLI GAF m=1 g=final_annotation.gff

# Sort
gt gff3 -sort filtered_predictions.gff > filtered_predictions_sorted.gff

#------------------------------------------------------------------------------
# Functional genome annotation using interproscan
#------------------------------------------------------------------------------

# Download and (protip) unzip tar .gz file with (SO much faster than gunzip + tar):
# tar -xf FILENAME.tar.gz
mkdir interproscan
cd interproscan
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.47-82.0/interproscan-5.47-82.0-64-bit.tar.gz
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.47-82.0/interproscan-5.47-82.0-64-bit.tar.gz.md5

# Recommended checksum to confirm the download was successful:
md5sum -c interproscan-5.47-82.0-64-bit.tar.gz.md5
# Must return *interproscan-5.47-82.0-64-bit.tar.gz: OK*
# If not - try downloading the file again as it may be a corrupted copy.
tar -pxvzf interproscan-5.47-82.0-*-bit.tar.gz

alias ips=/home/FM/celiason/uce-alcedinidae/interproscan-5.47-82.0/interproscan.sh

# Run
ips -i proteins_all.faa -goterms -f tsv

# Output protein sequences
gffread filtered_predictions.gff -g $ref -y proteins.fa

