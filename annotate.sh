# Annotate with gemoma


# Get genomes/annotations-

cd /home/FM/celiason/uce-alcedinidae/annotations/data

# calAnn
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/957/555/GCA_003957555.1_bCalAnn1_v1.p/GCA_003957555.1_bCalAnn1_v1.p_genomic.fna.gz

# Budgie
# wget ftp://ftp.ensembl.org/pub/rapid-release/gff3/melopsittacus_undulatus//Melopsittacus_undulatus.bMelUnd1.mat.Z.101.gff3.gz
# wget ftp://ftp.ensembl.org/pub/rapid-release/fasta/melopsittacus_undulatus/dna//Melopsittacus_undulatus.bMelUnd1.mat.Z.dna.toplevel.fa.gz

# aquChr
# wget ftp://ftp.ensembl.org/pub/release-100/fasta/aquila_chrysaetos_chrysaetos/dna/Aquila_chrysaetos_chrysaetos.bAquChr1.2.dna_sm.toplevel.fa.gz

#chicken
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.fna.gz
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.gff.gz



# round 1 = unmasked genome
# t=/home/FM/celiason/uce-alcedinidae/genomes/todChl/todChl.fasta \

# round 2 = repeat-masked genome
mkdir ~/uce-alcedinidae/annotations/gemoma/round2

cd ~/uce-alcedinidae/annotations/gemoma/round2

java -jar -Xms20G -XX:ParallelGCThreads=48 GeMoMa-1.7.1.jar CLI GeMoMaPipeline

infa=/home/FM/celiason/uce-alcedinidae/repeats/Full_mask/todChl.scaffolds.full_mask.fa
# head $infa

# Hard-mask to work with gemoma-
sed -e '/^>/! s/[[:lower:]]/N/g' $infa > masked.fasta
# head masked.fasta -n100

# java -jar -Xms20G -XX:ParallelGCThreads=48 GeMoMa-1.7.1.jar CLI GeMoMaPipeline

# Do the annotation
java -jar -Xms80G -XX:ParallelGCThreads=48 ../GeMoMa-1.7.1.jar CLI GeMoMaPipeline \
	restart=true \
	threads=48 \
	GeMoMa.Score=ReAlign \
	AnnotationFinalizer.r=NO \
	o=true \
	p=false \
	tblastn=false \
	t=/home/FM/celiason/uce-alcedinidae/annotations/gemoma/masked.fasta \
	s=own \
	i=melUnd \
	g=/home/FM/celiason/uce-alcedinidae/annotations/data/Melopsittacus_undulatus.bMelUnd1.mat.Z.dna.toplevel.fa \
	a=/home/FM/celiason/uce-alcedinidae/annotations/data/Melopsittacus_undulatus.bMelUnd1.mat.Z.101.gff3 \
	s=own \
	i=galGal \
	g=/home/FM/celiason/uce-alcedinidae/annotations/data/GCF_000002315.6_GRCg6a_genomic.fna \
	a=/home/FM/celiason/uce-alcedinidae/annotations/data/GCF_000002315.6_GRCg6a_genomic.gff \
	s=own \
	i=aquChr \
	g=/home/FM/celiason/uce-alcedinidae/annotations/data/Aquila_chrysaetos_chrysaetos.bAquChr1.2.dna_sm.toplevel.fa \
	a=/home/FM/celiason/uce-alcedinidae/annotations/data/Aquila_chrysaetos_chrysaetos.bAquChr1.2.101.gff3 \
	s=own \
	i=calAnn \
	g=/home/FM/celiason/uce-alcedinidae/annotations/data/GCF_003957555.1_bCalAnn1_v1.p_genomic.fna \
	a=/home/FM/celiason/uce-alcedinidae/annotations/data/GCF_003957555.1_bCalAnn1_v1.p_genomic.gff

# --mask-lower-case

gff=final_annotation.gff

cat $gff | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 - $4) } END { print NR, sum / NR }'

# 17001 genes, 17522.4 bp avg. (round1)
# 16539 genes, 17744.6 bp avg. (round2)

# java -jar GeMoMa-1.7.1.jar CLI Extractor c=true a=final_annotation.gff g=$ref

# Run gemoma annotation filter (GAF)
# Keep only best (top 1, m=1) matches-
java -jar ../GeMoMa-1.7.1.jar CLI GAF m=1 g=final_annotation.gff

# head round1/filtered_predictions.gff
# head round2/filtered_predictions.gff

# Sort
gt gff3 -sort round1/filtered_predictions.gff > round1/filtered_predictions_sorted.gff
gt gff3 -sort round2/filtered_predictions.gff > round2/filtered_predictions_sorted.gff

# Not sure what this does
gt eval round2/filtered_predictions_sorted.gff round1/filtered_predictions_sorted.gff

gffcompare -r round1/filtered_predictions.gff round2/filtered_predictions.gff

# wc -l round1/filtered_predictions.gff

#------------------------------------------------------------------------------
# Functional genome annotation
#------------------------------------------------------------------------------

# Using interproscan

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



# Test run-
./interproscan.sh -i test_proteins.fasta -f tsv
./interproscan.sh -i test_proteins.fasta -f tsv -dp

alias ips=/home/FM/celiason/uce-alcedinidae/interproscan-5.47-82.0/interproscan.sh
# ipsJar=/home/FM/celiason/uce-alcedinidae/interproscan/interproscan-5.47-82.0/interproscan-5.jar
# ips
# Real run-
# proteins.fasta file was produced by MAKER as standard output
ips -i proteins.fa -goterms -f tsv
# OK it's going. I think 1000 genes takes a day, so 17000 might take a while........
# NO! Much quicker- I started 11/26 3:06 PM; ended 11/27 12:36 AM (so <24 h!)

conda create --name MAKER
conda activate MAKER
conda install python=2.7
conda install -c bioconda maker
# conda create --name MAKER python=2.7.6
conda deactivate

# When you want to use maker, just activate the MAKER environment by
# source activate MAKER

# Add functional annotations to transcripts + GFF files
# on VORTEX server

# cd ~/kingfisher_annotation/todChl_rnd3.maker.output

ref="/home/FM/celiason/uce-alcedinidae/genomes/todChl/todChl.fasta"


cd /home/FM/celiason/uce-alcedinidae/annotations/gemoma/round2

alias ips=/home/FM/celiason/uce-alcedinidae/interproscan-5.47-82.0/interproscan.sh

ips -i proteins_all.faa -goterms -f tsv


# Output protein sequences
gffread filtered_predictions.gff -g $ref -y proteins.fa

# This once again is an example command line for running BLASTP against UniProt/Swiss-Prot:

# wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
# gunzip uniprot_sprot.fasta
# makeblastdb -dbtype prot -in uniprot_sprot.fasta
# blastp -query proteins.fa -db uniprot_sprot.fasta -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out output.blastp

# Using uniprot
# Filtered transcripts (N = 16540)
mmseqs easy-rbh --threads 48 proteins.fa ../uniprot_sprot.fasta uniprot.mmseqs tmp
# All transcripts (N = 21562)
mmseqs easy-rbh --threads 48 proteins_all.fa ../uniprot_sprot.fasta uniprot.all.mmseqs tmp

mmseqs easy-rbh --threads 48 todChl_augustus_proteins.faa ../uniprot_sprot.fasta todChl_augustus_proteins.faa.uniprot.mmseqs tmp


# Using chicken
mmseqs easy-rbh --threads 48 proteins.fa ../UP000000539_9031.fasta chicken.mmseqs tmp

# We will use the prefix 'GEMOMA' for our gene names, and an eight digit identifier.
source activate MAKER
maker_map_ids --prefix GEMOMA_ --justify 8 filtered_predictions.gff > filtered_predictions.map

# The output is a two column file translating old gene and mRNA names to new more standardized names.
# less $basename.map
#  maker-NT_010783.15-snap-gene-0.0	GMOD_00000001

# uniprot="/home/celiason/data/uniprot_sprot.fasta"

# These script do in-place replacement of names, so lets copy the files before running the scripts.
cp filtered_predictions.gff filtered_predictions.renamed.gff
# cp filtered_predictions.proteins.fasta filtered_predictions.proteins.renamed.fasta
# cp filtered_predictions.transcripts.fasta filtered_predictions.transcripts.renamed.fasta
# cp proteins.fa.tsv proteins.fa.renamed.tsv
# cp output.blastp output.renamed.blastp
cp uniprot.mmseqs uniprot.renamed.mmseqs
map_gff_ids filtered_predictions.map filtered_predictions.renamed.gff
# map_gff_ids filtered_predictions.map filtered_predictions.noseq.renamed.gff
# map_fasta_ids filtered_predictions.map filtered_predictions.proteins.renamed.fasta
# map_fasta_ids filtered_predictions.map filtered_predictions.transcripts.renamed.fasta
# map_data_ids filtered_predictions.map proteins.fa.renamed.tsv
# map_data_ids filtered_predictions.map output.renamed.blastp
map_data_ids filtered_predictions.map uniprot.renamed.mmseqs
head uniprot.renamed.mmseqs

# You can see names have changed by looking at the files.
# less $basename.noseq.gff
# less $basename.noseq.renamed.gff
# less $basename.map
# less output.renamed.blastp


# Use these commands to update your annotations with information from the BLAST report:
maker_functional_gff $uniprot output.renamed.blastp $basename.noseq.renamed.gff > $basename.noseq.renamed.putative_function.gff
maker_functional_fasta $uniprot output.renamed.blastp $basename.proteins.renamed.fasta > $basename.proteins.renamed.putative_function.fasta
maker_functional_fasta $uniprot output.renamed.blastp $basename.transcripts.renamed.fasta > $basename.transcripts.renamed.putative_function.fasta

# Look at the files to see that putative functions were added.

# less $basename.noseq.renamed.putative_function.gff
less $basename.transcripts.renamed.putative_function.fasta


# Finally we will add protein domain information to the final annotations using a report from InterProScan.
# This is done using the following scripts:

# ipr_update_gff - adds searchable tags to the gene and mRNA features in the GFF3 files
# iprscan2gff3 - adds physical viewable features for daomains that can be displayed in JBrowse, Gbrowse, and Web Apollo.

# This once again is an example command line for running InterProScan:

# interproscan.sh -appl pfam -dp -f TSV -goterms -iprlookup -pa -t p -i hsap_contig.maker.proteins.fasta -o output.iprscan

# Use these commands to update your annotations with information from the InterProScan report:

# head $basename.noseq.renamed.putative_function.gff

ipr_update_gff $basename.noseq.renamed.putative_function.gff output.renamed.iprscan > $basename.noseq.renamed.putative_function.domain_added.gff
cat $basename.noseq.renamed.putative_function.domain_added.gff | less

# iprscan2gff3 output.renamed.iprscan hsap_contig.renamed.gff > visible_iprscan_domains.gff 

conda deactivate

