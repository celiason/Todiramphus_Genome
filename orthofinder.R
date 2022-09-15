# Orthofinder run

annot="data/filtered_predictions.gff"
ref="data/todChl.scaffolds.full_mask.fa"

# Output protein sequences
gffread -J -V -y todChl_gemoma_proteins.faa -g $ref filtered_predictions.gff

# Setup orthofinder analysis
mkdir ortho; mkdir ortho/proteomes
cd ortho
ls proteomes # calAnn, melUnd, strHab, falRus, taeGut, corMon + todChl

# orthofinder -t 24 -a 24 -M msa -S mmseqs -T iqtree -f proteomes/

# Download  script to extract primary transcripts only (cuts down on size and runtime)
# wget https://raw.githubusercontent.com/davidemms/OrthoFinder/master/tools/primary_transcript.py primary_transcript2.py 

# Filter only primary transcripts
for f in proteomes/*faa ; do python primary_transcript2.py $f ; done
python primary_transcript2.py proteomes/todChl_augustus_proteins.faa

# Run orthofinder with coraciiform-heavy species set
orthofinder -t 24 -a 24 -M msa -og -S mmseqs -T iqtree -f proteomes/primary_transcripts/
