setwd("~/uce-alcedinidae/paml/allbirds")

# 1. figure out which Shultz genes correspond to our annotations 

# head ~/uce-alcedinidae/data/shultz/raw_results_genetrees_all_annotations.csv

meta <- read.csv("~/uce-alcedinidae/data/shultz/raw_results_genetrees_all_annotations.csv")

# /Users/chadeliason/Dropbox (The Field Museum)/Projects/King_genome/data/Shultz Sackton eLife R_data_and_scripts/compgen_alignments

# meta$hog
# meta$external_gene_name.x

# our genes

dat <- read.csv("~/uce-alcedinidae/paml/pos_selection_results_repeatmasked_20210504.csv")

hogs <- meta$hog[match(toupper(dat$finalgene), toupper(as.character(meta$external_gene_name.x)))]
table(is.na(hogs))  # 7633 good

tmp <- data.frame(file=as.character(dat$file), gene=dat$finalgene, hog=hogs, stringsAsFactors=F)
tmp <- tmp[!is.na(tmp$gene) & !is.na(tmp$hog), ]
nrow(tmp) # 6198 genes good to go

# 2. align sequences

i=1

f <- paste0("~/uce-alcedinidae/output/cds_FINAL/", tmp$file[i], ".fa")
raw <- readLines(f)
idx <- grep("todChl", raw)
part1 <- raw[c(idx,idx+1)]
part2 <- readLines(paste0("~/uce-alcedinidae/data/shultz/compgen_alignments/", tmp$hog[i], ".phy"))
part2 <- gsub("\\s+", "\n", part2)
part2 <- part2[-1]
part2 <- part2[-length(part2)]
part2 <- paste0(">", part2)
cat(part1,part2,sep="\n",file="test.fa")

fa <- "test.fa"

system(paste0("prank -d=", fa, " -codon -o=aln.fa"))  # took ~2 minutes

hyphy absrel --alignment aln.fa.best.fas --tree aves_foc_hyphy.tre --branches Foreground


# 3. fit M0 models

runCodeml(fasta="/home/FM/celiason/uce-alcedinidae/paml/allbirds/aln.fa.best.fas", model="M0", phy="/home/FM/celiason/uce-alcedinidae/paml/allbirds/aves.tre", force=TRUE)


# 4. fit branch-site model alternative



# 5. fit Anull model



# 6. do likelihood comparisons


# 7. run enrichment test to find novel genes in todiramphus




# 8. Orthofinder??
cd annotations/gemoma/round2/

ref="/home/FM/celiason/uce-alcedinidae/repeats/Full_mask/todChl.scaffolds.full_mask.fa"
echo $ref
gffread -J -V -y proteins.faa -g $ref filtered_predictions.gff
# gffread -J -V -y todChl_augustus_proteins.faa -g $ref todChl_augustus_final.gff3

grep -c ">" proteins.faa
# head proteins_all.faa
# cp proteins_all.faa /home/FM/celiason/uce-alcedinidae/ortho/proteomes/todChl_gemoma_proteins.faa
cp proteins.faa /home/FM/celiason/uce-alcedinidae/ortho/proteomes/todChl_gemoma_proteins.faa

python primary_transcript2.py

# cd ~/uce-alcedinidae
# mkdir ortho; mkdir ortho/proteomes
cd ~/uce-alcedinidae/ortho
ls proteomes # calAnn, melUnd, strHab, falRus, taeGut, corMon + todChl
# orthofinder -t 24 -a 24 -M msa -S mmseqs -T iqtree -f proteomes/

# download  script to extract primary transcripts only (cuts down on size and runtime)
# wget https://raw.githubusercontent.com/davidemms/OrthoFinder/master/tools/primary_transcript.py primary_transcript2.py 
for f in proteomes/*faa ; do python primary_transcript2.py $f ; done
python primary_transcript2.py proteomes/todChl_augustus_proteins.faa

# Run with coraciiform-heavy species set
orthofinder -t 24 -a 24 -M msa -og -S mmseqs -T iqtree -f proteomes/primary_transcripts/

# Not sure why I had to do this?
# orthofinder -t 24 -a 24 -og -S mmseqs -b proteomes/primary_transcripts/OrthoFinder/Results_May06/WorkingDirectory -f redo

# Final files are here:
# ~/uce-alcedinidae/ortho/proteomes/primary_transcripts/OrthoFinder/Results_May06/WorkingDirectory/OrthoFinder/Results_May06/

scp celiason@10.100.111.202:~/uce-alcedinidae/ortho/proteomes/primary_transcripts/OrthoFinder/Results_May06/WorkingDirectory/OrthoFinder/Results_May06/Comparative_Genomics_Statistics/Statistics_Overall.tsv ~/Dropbox/Projects/King_genome/output/orthofinder_may06
# Feb 4 2022 - Rerun with only high-quality VGP genomes
orthofinder -t 24 -a 24 -M msa -og -S mmseqs -T iqtree -f proteomes/primary_transcripts/new_vgp


# Venn diagrams
library(VennDiagram)
library(scales)
library(RColorBrewer)

x <- read.delim("~/Downloads/Orthogroups.GeneCount.tsv")
head(x)

tmp <- list(todMex=ifelse(x[,2]>0, x[,1], NA),
			halSen=ifelse(x[,3]>0, x[,1], NA),
			chlAen=ifelse(x[,4]>0, x[,1], NA),
			ceyCya=ifelse(x[,5]>0, x[,1], NA),
			calAnn=ifelse(x[,6]>0, x[,1], NA),
			taeGut=ifelse(x[,7]>0, x[,1], NA),
			todChl=ifelse(x[,8]>0, x[,1], NA))
tmp <- lapply(tmp, function(x) x[!is.na(x)])

picks <- c(1,2,3,6,7)
pal <- brewer.pal(length(picks), "Set2")
venn.diagram(tmp[picks], filename="poop.jpg", col=pal, fill=alpha(pal, .25))

library(ape)
tr <- read.tree(text="((((taeGut,melUnd),falRus),todChl), calAnn);")
plot(tr)

