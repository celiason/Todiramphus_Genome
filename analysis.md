# Collared kingfisher draft genome assembly and analysis

## Repeat masking

```sh
mkdir repeats
mkdir repeats/todChl_mask
mkdir repeats/galGal_mask
mkdir repeats/Full_mask

# De Novo Repeat Identification
RepeatModeler/BuildDatabase -name todChl -engine ncbi genomes/todChl.fasta

# ./RepeatModeler --help
RepeatModeler -pa 15 -engine ncbi -database todChl 2>&1 | tee repeatmodeler.log

# Soft-mask using de novo library-
RepeatMasker -pa 48 -xsmall -e ncbi -lib repeats/todChl-families.fa -dir repeats/todChl_mask genomes/todChl.fasta

# Then the maksed FASTA from this search can be used as input for the next search, using the chicken library from Repbase. I also normally rename the outputs after each round so they are more representative of what they contain.
RepeatMasker -pa 48 -xsmall -e ncbi -species chicken -dir repeats/galGal_mask repeats/todChl_mask/todChl.fasta.masked

# Finally, results from each round must be analyzed together to produce the final repeat annotation
mkdir repeats/Full_mask
cp repeats/galGal_mask/todChl.fasta.masked.masked genomes/todChl.scaffolds.full_mask.fa

# Output full-masked for annotation in genoma
gzip -k repeats/Full_mask/todChl.scaffolds.full_mask.fa

# Create BED files
./RMout_to_bed.pl repeats/galGal_mask/todChl.fasta.masked.out base0
./RMout_to_bed.pl repeats/todChl_mask/todChl.fasta.out base0

# Find intersection between 2 masks
bedtools intersect -a repeats/galGal_mask/todChl.fasta.masked.out.bed -b repeats/todChl_mask/todChl.fasta.out.bed > repeats/Full_mask/intersect.bed

# Get statistics on overlap

wc -l repeats/galGal_mask/todChl.fasta.masked.out.bed  # 549083 - 98% overlapping with chicken repeats
wc -l repeats/todChl_mask/todChl.fasta.out.bed  # 680766 - 79% de novo
wc -l repeats/Full_mask/intersect.bed  # 536414 common to both
```

## Structural annotation with GeMoMa

```sh
# Anna's hummingbird genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/957/555/GCA_003957555.1_bCalAnn1_v1.p/GCA_003957555.1_bCalAnn1_v1.p_genomic.fna.gz -o genomes/

# Budgerigar genome and annotation
wget ftp://ftp.ensembl.org/pub/rapid-release/gff3/melopsittacus_undulatus//Melopsittacus_undulatus.bMelUnd1.mat.Z.101.gff3.gz -o genomes/
wget ftp://ftp.ensembl.org/pub/rapid-release/fasta/melopsittacus_undulatus/dna//Melopsittacus_undulatus.bMelUnd1.mat.Z.dna.toplevel.fa.gz -o annotation/

# Golden eagle genome
wget ftp://ftp.ensembl.org/pub/release-100/fasta/aquila_chrysaetos_chrysaetos/dna/Aquila_chrysaetos_chrysaetos.bAquChr1.2.dna_sm.toplevel.fa.gz -o genomes/

# Chicken genome and annotation
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.fna.gz -o genomes/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.gff.gz -o annotation/

infa=genomes/todChl.scaffolds.full_mask.fa

# Hard-mask to work with gemoma-
sed -e '/^>/! s/[[:lower:]]/N/g' $infa > genomes/todChl_hardmasked.fasta

java -jar -Xms20G -XX:ParallelGCThreads=48 GeMoMa-1.7.1.jar CLI GeMoMaPipeline

# Working-
java -jar -Xms80G -XX:ParallelGCThreads=48 ../GeMoMa-1.7.1.jar CLI GeMoMaPipeline \
	restart=true \
	threads=48 \
	GeMoMa.Score=ReAlign \
	AnnotationFinalizer.r=NO \
	o=true \
	p=false \
	tblastn=false \
	t=genomes/todChl_hardmasked.fasta \
	s=own \
	i=melUnd \
	g=genomes/Melopsittacus_undulatus.bMelUnd1.mat.Z.dna.toplevel.fa \
	a=annotation/Melopsittacus_undulatus.bMelUnd1.mat.Z.101.gff3 \
	s=own \
	i=galGal \
	g=genomes/GCF_000002315.6_GRCg6a_genomic.fna \
	a=annotation/GCF_000002315.6_GRCg6a_genomic.gff \
	s=own \
	i=aquChr \
	g=genomes/Aquila_chrysaetos_chrysaetos.bAquChr1.2.dna_sm.toplevel.fa \
	a=annotation/Aquila_chrysaetos_chrysaetos.bAquChr1.2.101.gff3 \
	s=own \
	i=calAnn \
	g=genomes/GCF_003957555.1_bCalAnn1_v1.p_genomic.fna \
	a=annotation/GCF_003957555.1_bCalAnn1_v1.p_genomic.gff

# --mask-lower-case

# 16539 genes, 17744.6 bp avg. in resulting annotation-
gff=final_annotation.gff
cat $gff | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 - $4) } END { print NR, sum / NR }'

# Run gemoma annotation filter (GAF) to keep only best (top 1, m=1) matches
java -jar ../GeMoMa-1.7.1.jar CLI GAF m=1 g=final_annotation.gff
gt gff3 -sort filtered_predictions.gff > filtered_predictions_sorted.gff

# Output protein sequences
gffread final_annotation.gff -g $ref -y proteins_all.faa
gffread filtered_predictions.gff -g $ref -y proteins.faa

# We will use the prefix 'GEMOMA' for our gene names, and an eight digit identifier.
maker_map_ids --prefix GEMOMA_ --justify 8 filtered_predictions.gff > filtered_predictions.map
# The output is a two column file translating old gene and mRNA names to new more standardized names.
# head $basename.map
# maker-NT_010783.15-snap-gene-0.0	GMOD_00000001

# These script do in-place replacement of names, so lets copy the files before running the scripts.
cp filtered_predictions.gff filtered_predictions.renamed.gff
# cp todChl_uniprot.mmseqs todChl_uniprot.renamed.mmseqs
map_gff_ids filtered_predictions.map filtered_predictions.renamed.gff
# map_data_ids filtered_predictions.map uniprot.renamed.mmseqs
```

## Functional annotation

```r
mkdir interproscan
cd interproscan

# Using interproscan
# Download and (protip) unzip tar .gz file with (SO much faster than gunzip + tar):
# tar -xf FILENAME.tar.gz
# wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.47-82.0/interproscan-5.47-82.0-64-bit.tar.gz
# wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.47-82.0/interproscan-5.47-82.0-64-bit.tar.gz.md5

# Recommended checksum to confirm the download was successful:
# md5sum -c interproscan-5.47-82.0-64-bit.tar.gz.md5
# Must return *interproscan-5.47-82.0-64-bit.tar.gz: OK*
# If not - try downloading the file again as it may be a corrupted copy.

# tar -pxvzf interproscan-5.47-82.0-*-bit.tar.gz

alias ips=/home/FM/celiason/uce-alcedinidae/interproscan-5.47-82.0/interproscan.sh

# Run
ref=genomes/todChl.fasta
ips -i annotation/proteins_all.faa -goterms -f tsv

# Using uniprot

# Filtered transcripts (N = 16540)
mmseqs easy-rbh --threads 48 proteins.fa ../uniprot_sprot.fasta uniprot.mmseqs tmp

# All transcripts (N = 21562)
mmseqs easy-rbh --threads 48 proteins_all.fa ../uniprot_sprot.fasta uniprot.all.mmseqs tmp

# Using chicken
mmseqs easy-rbh --threads 48 proteins.fa ../UP000000539_9031.fasta chicken.mmseqs tmp
```

## Genome consistency plots

```sh
# Chicken
ref=genomes/galGal3_all.fa
./JupiterPlot/jupiter2 name=todChl-to-galGal3 ref=$ref fa=genomes/todChl.fasta

# Hornbill
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/769/605/GCA_009769605.1_bBucAby1.pri/GCA_009769605.1_bBucAby1.pri_genomic.fna.gz
# gunzip GCA_009769605.1_bBucAby1.pri_genomic.fna.gz
ref=genomes/GCA_009769605.1_bBucAby1.pri_genomic.fna
./JupiterPlot/jupiter2 name=todChl-to-bucAby ref=$ref fa=genomes/todChl.fasta

# Bee eater
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/819/595/GCA_009819595.1_bMerNub1.pri/GCA_009819595.1_bMerNub1.pri_genomic.fna.gz
# gunzip GCA_009819595.1_bMerNub1.pri_genomic.fna.gz
ref=GCA_009819595.1_bMerNub1.pri_genomic.fna
./JupiterPlot/jupiter2 name=todChl-to-merNub ref=$ref fa=genomes/todChl.fasta
```

## PSMC analysis

### Align raw reads to reference

```sh
# Combine raw reads (from SRA) into single files-
cat King*R1*fastq.gz > todChl_reads1.fastq.gz
cat King*R2*fastq.gz > todChl_reads2.fastq.gz

# Setup analysis
reads1="todChl_reads1.fastq.gz"
reads2="todChl_reads2.fastq.gz"
ref="ref/todChl.fasta"

# Run
alias bwa=./bwa/bwa
alias fastp=./fastp
fastp -i $reads1 -I $reads2 --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --stdout -h todChl-to-todChl.html |\
# align
bwa mem -p -t 48 genomes/kingfisher - |\
# make BAM
samtools view -Sb - |\
# sort BAM and save (using 100 GB of RAM with -m argument)
samtools sort -m 100G > todChl-to-todChl.bam

# Index BAM file
samtools index todChl-to-todChl.bam

# Get statistics
samtools flagstat 

# View alignment
samtools tview todChl-to-todChl.bam genomes/todChl.fasta
```

### Call SNVs in R

```r
library(parallel)

ref="genomes/todChl.fasta"
bam="todChl-to-todChl.bam"
# regions <- system(paste0("grep -Po '(?<=>)[0-9_]+' ", ref), intern=TRUE)

# run in shell-
# cat reference_genomes/kingfisher.fasta | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > kingfisher.fasta.seqlengths

# Create regions for parallelizing
nms <- read.delim("kingfisher.fasta.seqlengths", sep="\t", head=F)
nms$seqid <- stringr::str_extract(nms[,1], "^\\d+")
res <- lapply(1:nrow(nms), function(x) {
	if (nms$V2[x] > 5e6) {
		starts <- seq(1, nms$V2[x], by=5e6)
		ends <- c(starts[-1]-1, nms$V2[x])
	} else {
		starts <- 1
		ends <- nms$V2[x]
	}
	data.frame(seq=nms$seqid[x], start=starts, end=ends, id=paste0(nms$seqid[x], "_", 1:length(starts)))
})
options(scipen=999)
res2 <- do.call(rbind, res)
regions <- setNames(paste0(res2$seq, ":", res2$start, "-", res2$end), res2$id)

mclapply(1:length(regions), mc.cores=48, function(x) {
	reg <- unname(regions[x])
	filename <- names(regions[x])
	run <- paste0("samtools mpileup -Q 30 -q 30 -r ", reg, " -uf ", ref, " ", bam, " | bcftools call -c | vcfutils.pl vcf2fq -d 10 -D 100 | gzip > temp/", filename, ".fq.gz")
	system(run)
})
```

### Cleanup and merge fq files

```sh
# Delete files <100 bytes
find *.gz -type f -size -100b -delete

# Combine
cat temp/*.gz > todChl_diploid.fq.gz
```

### Prepare for PSMC

```r
# Split kingfisher10X genome into smaller fragments/contigs:
system("./psmc/utils/splitfa genomes/todChl.fasta > genomes/todChl_split.fa")
```

### Run PSMC

```sh
# Make psmcfa file - quick
./psmc/utils/fq2psmcfa todChl_diploid.fq.gz > todChl_diploid.psmcfa
# cat todChl_diploid.psmcfa

# Run PSMC analyses
# Using settings in nadachowska-bryska et al 2015: "Based on the results from the pilot runs, we chose the final settings for the PSMC to be ‘‘N30 –t5 –r5 –p 4+30*2+4+6+10’’ for all species.""
./psmc/psmc -N30 -t5 -r5 -p"4+30*2+4+6+10" -o todChl_diploid.psmc todChl_diploid.psmcfa
```

### PSMC bootstraps

```r
# To perform bootstrapping, one has to run splitfa first to split long chromosome
# sequences to shorter segments. When the `-b' option is applied, psmc will then
# randomly sample with replacement from these segments. As an example, the
# following command lines perform 100 rounds of bootstrapping:
library(parallel)
system("./psmc/utils/splitfa todChl_diploid.psmcfa > todChl_split.psmcfa")
# Run
mclapply(1:100, mc.cores=48, function(x) {
	run <- paste0("./psmc/psmc -N30 -t5 -r5 -b -p'4+30*2+4+6+10' -o temp/todChl_round-", x, ".psmc todChl_split.psmcfa")
	system(run)
})

# Combine
system("cat todChl_diploid.psmc temp/todChl_round-*.psmc > todChl_combined.psmc")

# Cleanup
system("rm temp/todChl_round-*.psmc")
```

## Orthofinder analysis

```sh
ref=genomes/todChl.scaffolds.full_mask.fa

gffread -J -V -y proteins.faa -g $ref filtered_predictions.gff

grep -c ">" proteins.faa


# Download proteomes (see manuscript, supp. tables for details/NCBI links)
# mkdir proteomes; cd proteomes

ls proteomes # calAnn, melUnd, strHab, falRus, taeGut, corMon + todChl

# Extract primary transcripts only (cuts down on size and runtime)
# wget https://raw.githubusercontent.com/davidemms/OrthoFinder/master/tools/primary_transcript.py primary_transcript2.py 
for f in proteomes/*faa
	do
		python primary_transcript2.py $f
	done

python primary_transcript2.py annotation/proteins.faa

# Run with coraciiform-heavy species set
orthofinder -t 24 -a 24 -M msa -og -S mmseqs -T iqtree -f proteomes/primary_transcripts/
```

## CAFE

### Run

```sh
# Estimate error to account for incomplete assemblies (e.g., for B10k genomes)
python cafe/caferror.py -i kingfishers/cafe_error1.sh -d kingfishers/error -f 1 -v 0 -s 1

# Now run accounting for error
cafe cafe_king1.sh
```

### Analysis

```r
library(ape)
library(stringr)
library(phytools)
library(readxl)
library(stringr)
library(UpSetR)

# CAFE3- "An error value of ε = 0.1 means that in 90% of gene families, the observed size is equal to the true size, whereas in 10% of gene families, the observed size is either too large or too small (fig. 1A, C, and E)."

# Plot tree
tr <- read.tree("coraci192-beast-fixed.tre")  # downloaded from dryad (Eliason et al. 2021 Am Nat)
tr <- drop.tip(tr, which(!(tr$tip %in% c('Todiramphus_chloris','Todus_mexicanus','Halcyon_senegalensis','Chloroceryle_aenea','Ceyx_cyanopectus'))))
source("sppabbrev.R")
tr$tip.label <- sppAbbrev(tr$tip)
tr$root.edge <- 25
maxage <- max(branching.times(tr))
tr <- bind.tip(tr, tip.label='taeGut', position=62-maxage, edge.length=62)
maxage <- max(branching.times(tr))
tr <- bind.tip(tr, tip.label='calAnn', position=67.4-maxage, edge.length=67.4)
plot(tr)
is.ultrametric(tr)


resfile <- "cafe/cafe_final_report.cafe"
# resfile <- "/Users/chadeliason/Downloads/CAFE/kingfishers/cafe_final_report.cafe"

raw <- readLines(resfile)

# avg expansion
# avgexp <- as.numeric(str_match_all(raw[6], "[0-9\\.\\-]+")[[1]][, 1])

# caford <- rank(as.numeric(str_match_all(raw[4], "[0-9]+")[[1]][,1]))
caford <- as.character(as.numeric(str_match_all(raw[4], "[0-9]+")[[1]][,1]))
# names(avgexp) <- caford

tr <- read.tree(text=paste0(gsub("Tree:", "", raw[1]), ";"))

phyord <- c('0'=1,'2'=2,'4'=3,'6'=4,'8'= 5,'7'=13,'5'=12,'3'=11,'1'=10,'10'=6,'9'=9,'12'=7,'11'=8)

# Significant changes only
out <- read.delim(resfile, head=T, skip=9, sep="\t")
out$padj <- p.adjust(out[, 3], method="fdr")
# table(out$padj < 0.05)
keep <- out[, 3] < 0.05  # significant ones
table(keep)  # N = 1443 orthogroups

# sig <- as.data.frame(do.call(rbind, strsplit(gsub("\\(|\\)", "", out[keep, 4]), split=",")))
# plot(sig[, 1])
# abline(h=0.05, lty=2)

# Get family-wide P values of gene families expanded/contracted
res <- t(sapply(1:nrow(out), function(x) str_match_all(out[x, 2], "_(\\d+)\\:")[[1]][, 2] ))
res <- apply(res, 2, as.numeric)
colnames(res) <- str_match_all(raw[3], "<(\\d+)>")[[1]][, 2]
colnames(res) <- phyord[colnames(res)]  # 
res <- res[, order(as.numeric(colnames(res)))]
rownames(res) <- paste0("og", 1:nrow(res))

# Ids for tips and nodes to use later
tid <- as.character(1:Ntip(tr))
nid <- as.character((Ntip(tr)+1) : (Ntip(tr) + Nnode(tr)))

# Node-wise P values-
nodesig <- as.data.frame(do.call(rbind, strsplit(gsub("\\(|\\)", "", out[keep, 4]), split=",")))
colnames(nodesig) <- phyord[caford]

# Only keep certain OGs
res.signif <- res[keep, ]

# Calculate changes along branches of a phylogeny
changes <- t(sapply(1:nrow(res.signif), function(x) {
	states <- res.signif[x, ]
	changes <- states[tr$edge[, 2]] - states[tr$edge[, 1]]
	pvals <- nodesig[x, match(names(changes), colnames(nodesig))]  # Viterbi P values from cafe
	ifelse(pvals < 0.05, changes, 0)
}))
rownames(changes) <- rownames(res.signif)
colnames(changes) <- tr$edge[, 2]

sum_nochange <- apply(changes, 2, function(x) sum(x == 0))  # no change
sum_expand <- apply(changes, 2, function(x) sum(x > 0))  # expand
sum_contract <- apply(changes, 2, function(x) sum(x < 0))  # contract

tot_expand <- apply(changes, 2, function(x) sum(x[x > 0]))  # expand
tot_contract <- apply(changes, 2, function(x) sum(x[x < 0]))  # contract

change_per_gene <- apply(changes, 2, sum) / nrow(res.signif)

pdf("figs/gene_expansion_err_signif_viter.pdf", width=3, height=6)
plot(tr, edge.col = ifelse(change_per_gene > 0, "red", "blue"), no.margin = TRUE, cex = 0.75)
edgelabels(text = round(change_per_gene, 3), frame = "none", adj = c(0.5, -0.5), cex = 0.7)
edgelabels(text = paste0("+", sum_expand, "/ -", sum_contract), frame = "none", adj = c(0.5, 1.5), cex = 0.7)
dev.off()

# Load orthogroup output
ogs <- read.delim("cafe/Orthogroups.tsv")

# Significantly expanded in todChl
picks <- names(which((changes[, "5"] > 0)))
picks <- as.numeric(str_extract(picks, "\\d+"))
ogs[picks,]
transcripts_expanded <- unlist(strsplit(ogs[picks, 8], ", "))  # Convert to list of transcript names
picks <- names(which((changes[, "5"] < 0)))
picks <- as.numeric(str_extract(picks, "\\d+"))
transcripts_contracted <- unlist(strsplit(ogs[picks, 8], ", "))

# Load positive selection/gene name annotations
nms <- read.csv("nms.csv", row=1)
nms$file <- gsub(":", "_", nms$file)

# Get gene name list
genes <- nms$finalgene[match(transcripts_expanded, nms$file)]
genes <- genes[!is.na(genes)]
genes <- sort(unique(genes))
clipr::write_clip(genes)  # put this into STRING network website for analysis

# Load interproscan annotation
iprdat <- read.delim("annotation/todChl_proteins_all_gemoma_faa.tsv", head=FALSE)
iprdat$V1 <- gsub(":", "_", iprdat$V1)  # fix transcript names

# todiramphus green color (RGB) in map- 72	170	129	

# Get list of GO terms

goterms1 <- iprdat$V14[match(transcripts_expanded, iprdat$V1)]
goterms1 <- goterms1[!is.na(goterms1) & goterms1!=""]
goterms1 <- strsplit(goterms1, "\\|")
goterms1 <- sort(unique(unlist(goterms1)))
clipr::write_clip(unlist(goterms1)) # for use in REVIGO (Fig. 3B panel)
clipr::write_clip(unlist(goterms1), breaks=",") # for use in REVIGO (Fig. 3B panel)

# Genes associated with each term-
picks <- iprdat$V1 %in% transcripts_expanded & grepl("GO:0007186", iprdat$V14) # GPCR
iprdat[picks, ]

picks <- iprdat$V1 %in% transcripts_expanded & grepl("GO:0007186", iprdat$V14) & grepl("[Oo]psin", iprdat$V13)# GPCR
iprdat[picks, ]
sort(unique(iprdat$V1[picks]))  # 34 transcripts

gene_lookup$finalgene[match(unique(iprdat$V1[picks]), gene_lookup$file)]

goterms2 <- iprdat$V14[match(transcripts_contracted, iprdat$V1)]
goterms2 <- goterms2[!is.na(goterms2) & goterms2!=""]
goterms2 <- strsplit(goterms2, "\\|")
goterms2 <- sort(unique(unlist(goterms2)))
clipr::write_clip(unlist(goterms2)) # for use in REVIGO (Fig. 3B panel)

# interproscan results 
iprnums <- iprdat[match(transcripts, iprdat$V1), "V12"]
iprnums <- iprnums[!is.na(iprnums) & iprnums!=""]
```

### Load gene families

```r
x <- strsplit(ogs[,8], ",")

# Find families that have expanded in todChl
res <- read.delim("cafe/report.txt.cafe", skip=10)
trees <- res[,2]

# Find todChl-halSen branch in tree
m <- str_match(trees, "halSen_(\\d+).*?todChl_(\\d+).*?\\)_(\\d+)")

# Get number of unique gene families expanded in todChl
x1 <- as.numeric(m[,2])  # halSen
x2 <- as.numeric(m[,3])  # todChl
x3 <- as.numeric(m[,4])  # ancestor

# Get number of unique gene families expanded in ceyCya
m2 <- str_match(trees, "ceyCya\\_(\\d+)\\:.*?(\\d+)\\:16\\.7776")

ogstats <- read.delim("Orthogroups.GeneCount.tsv")
names(ogstats)[2:8] <- c('todMex','halSen','chlAen','ceyCya','calAnn','taeGut','todChl')

idx <- which(x2>x3 & ogstats[,"ceyCya"]!=0)

# Find todChl genes in orthogroups gained in todChl relative to MRCA with halSen
transcripts <- unlist(strsplit(ogs[which(x2 > x3), 8], ", "))

setdat <- ogstats[,2:8]
setdat <- cbind(set='a',setdat)
setdat[,-1] <- ifelse(setdat[,-1] >= 1, 1, 0)

# Upset plot for manuscript
pdf("figs/todChl_upset.pdf", width=6.5, height=4.5)
upset(setdat, nsets=7, sets=c('todMex','ceyCya','chlAen','halSen','todChl','taeGut','calAnn'), nintersects=35, order.by="freq", mainbar.y.label="Number of orthogroups", keep.order=TRUE,
		text.scale=c(1,1,1,1,1,.75))
dev.off()

# Lookup genes from GFF files used in gemoma-
annot <- read.delim("~/Downloads/uniprot_todChl_annot_gemoma_round2.txt", head=F)
annot$V1 <- gsub(":", "_", annot$V1)

# Get gene names
genes <- nms$finalgene[match(transcripts, nms$file)]
genes <- genes[!is.na(genes)]
genes <- sort(unique(genes))
clipr::write_clip(genes)  # put this into STRING network website for analysis

# Functional annotation of genes
dat <- read.delim("Orthogroups.GeneCount.tsv")

# Orthogroups specific to todChl-
picks <- as.character(dat[apply(dat[,2:7], 1, sum)==0, "Orthogroup"])
picks <- paste0("Orthogroup_Sequences/", picks, ".fa")
seqs <- lapply(picks, readLines)
# tmp <- read.csv("~/uce-alcedinidae/paml/pos_selection_results_repeatmasked_20210504.csv")
# tmp$finalgene[match(gsub(">", "", grep(">", unlist(seqs), value=T)), tmp$file)]
cat(unlist(seqs), sep="\n", file="~/uce-alcedinidae/ogs_todChl.faa")

# 1) Functional annotations of this set of 404 genes using MMSEQS-
# mmseqs easy-rbh --threads 12 ogs_todChl.faa annotation/gemoma/uniprot_sprot.fasta ogs_todChl_uniprot.mmseqs tmp
# head ogs_todChl_uniprot.mmseqs

# Gene list

# 2) Functional annotation using InterProScan-
alias ips=/home/FM/celiason/uce-alcedinidae/interproscan-5.47-82.0/interproscan.sh
# ipsJar=/home/FM/celiason/uce-alcedinidae/interproscan/interproscan-5.47-82.0/interproscan-5.jar
# Real run-
# proteins.fasta file was produced by MAKER as standard output
# ips -i ogs_todChl.faa -goterms -f tsv
# head ogs_todChl.faa_1.tsv

# Load IPR annotated genes
dat <- read.delim("output/ogs_todChl.faa_1.tsv", head=F)

# Run enrichment analysis
library(dcGOR)
res <- dcEnrichment(dat$V12, domain="InterPro", ontology="GOBP")

# Visualize as a graph
g <- visEnrichment(res, node.info="full_term_name")
g2 <- visEnrichment(res, node.info="full_term_name", num_top_nodes=10)

# Export terms for REVIGO analysis
tops <- view(res, 125)
cat(tops$term_id, sep="\n")
```

### R code to plot in REVIGO-

```r
# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0008150","biological_process",100,1,0,"biological_process"),
c("GO:0008152","metabolic process",62.693911296527,1,0,"metabolic process"),
c("GO:0009628","response to abiotic stimulus",0.603634500970981,0.925217811316376,0,"response to abiotic stimulus"),
c("GO:0009416","response to light stimulus",0.188810302859614,0.843634660246174,0.30112543,"response to abiotic stimulus"),
c("GO:0051716","cellular response to stimulus",11.0866663899183,0.89014938374801,0.46891116,"response to abiotic stimulus"),
c("GO:0009987","cellular process",78.0376878191047,1,0,"cellular process"),
c("GO:0043933","protein-containing complex subunit organization",1.6681875395806,0.935840927896943,0,"protein-containing complex subunit organization"),
c("GO:0006996","organelle organization",4.10954278773561,0.929874097335854,0.66494786,"protein-containing complex subunit organization"),
c("GO:0050794","regulation of cellular process",19.8740720541611,0.97184669100809,0,"regulation of cellular process"),
c("GO:0050896","response to stimulus",13.7208133053304,1,0,"response to stimulus"),
c("GO:0051179","localization",18.7941114253392,1,0,"localization"),
c("GO:0065007","biological regulation",23.3881336449155,1,0,"biological regulation"),
c("GO:1901135","carbohydrate derivative metabolic process",7.03889068001952,0.911362892805519,0,"carbohydrate derivative metabolic process"),
c("GO:1901360","organic cyclic compound metabolic process",22.4511996998524,0.878150523573222,0.12633896,"carbohydrate derivative metabolic process"),
c("GO:0043170","macromolecule metabolic process",34.4366477094063,0.85874161483324,0.20469617,"carbohydrate derivative metabolic process"),
c("GO:1901576","organic substance biosynthetic process",23.3073421911027,0.837907478511037,0.20773315,"carbohydrate derivative metabolic process"),
c("GO:1901564","organonitrogen compound metabolic process",32.1987066241535,0.841282333392547,0.23825687,"carbohydrate derivative metabolic process"),
c("GO:0044267","cellular protein metabolic process",15.254359469727,0.823898453789706,0.2866446,"carbohydrate derivative metabolic process"),
c("GO:0006139","nucleobase-containing compound metabolic process",18.7768975137146,0.83816886145711,0.30786996,"carbohydrate derivative metabolic process"),
c("GO:0034641","cellular nitrogen compound metabolic process",26.2084518110319,0.83172830309132,0.34939378,"carbohydrate derivative metabolic process"),
c("GO:0019538","protein metabolic process",19.8225173210066,0.837878750066921,0.43922804,"carbohydrate derivative metabolic process"),
c("GO:0044260","cellular macromolecule metabolic process",25.5354976543922,0.808874814264069,0.46854268,"carbohydrate derivative metabolic process"),
c("GO:0009059","macromolecule biosynthetic process",11.1039548744453,0.840549362958524,0.58291612,"carbohydrate derivative metabolic process"),
c("GO:0016192","vesicle-mediated transport",1.41412180423006,0.868094089715752,0.01304146,"vesicle-mediated transport"),
c("GO:0051649","establishment of localization in cell",2.04442026072841,0.861975425874124,0.34310922,"vesicle-mediated transport"),
c("GO:0051641","cellular localization",2.42716153907533,0.861405482004253,0.35682686,"vesicle-mediated transport"),
c("GO:0006812","cation transport",3.53439099252141,0.850009956472001,0.38656461,"vesicle-mediated transport"),
c("GO:0006810","transport",17.9971715326728,0.827491774029602,0.5528133,"vesicle-mediated transport"),
c("GO:0006811","ion transport",5.1712786278122,0.856725097409823,0.59779702,"vesicle-mediated transport"),
c("GO:0071702","organic substance transport",5.81485106092797,0.85425991410083,0.61317491,"vesicle-mediated transport"),
c("GO:0071840","cellular component organization or biogenesis",9.05104158872763,0.981371738760136,0.0167685,"cellular component organization or biogenesis"),
c("GO:0044281","small molecule metabolic process",15.7742278866694,0.9081532110822,0.07740644,"small molecule metabolic process"),
c("GO:0009058","biosynthetic process",24.4296768002127,0.894418865175645,0.10698742,"small molecule metabolic process"),
c("GO:0006807","nitrogen compound metabolic process",45.3621081991498,0.866258630750214,0.15835912,"small molecule metabolic process"),
c("GO:0044238","primary metabolic process",49.9247269409607,0.860620178906147,0.23456847,"small molecule metabolic process"),
c("GO:0071704","organic substance metabolic process",55.5008618556042,0.853950970497847,0.27143659,"small molecule metabolic process"),
c("GO:0044237","cellular metabolic process",51.5163612533451,0.848707235179424,0.27824028,"small molecule metabolic process"),
c("GO:0006796","phosphate-containing compound metabolic process",13.8765960983292,0.85407087781059,0.091154,"phosphate-containing compound metabolic process"),
c("GO:0006793","phosphorus metabolic process",14.1931042105766,0.875515013458844,0.16376197,"phosphate-containing compound metabolic process"),
c("GO:0046483","heterocycle metabolic process",21.6400578122282,0.859605123703167,0.18465362,"phosphate-containing compound metabolic process"),
c("GO:0006725","cellular aromatic compound metabolic process",21.6625912575535,0.859560879108118,0.21016732,"phosphate-containing compound metabolic process"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );
```

### Output figures

```r
# by default, outputs to a PDF file
pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches
# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "uniqueness",
  type = "categorical",
  vColor = "representative",
  title = "Revigo TreeMap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)
dev.off()

jpg("figs/network_ogs.jpg", width=7, height=7)
plot(g, font=16)
dev.off()
```
