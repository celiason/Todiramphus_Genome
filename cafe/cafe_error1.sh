# !/Users/chadeliason/Downloads/CAFE/release/cafe
# date
# version

#specify data file, p-value threshold, # of threads to use, and log file
load -i Orthogroups.GeneCount.tsv -p 0.01 -t 4 -l log.txt

#the phylogenetic tree structure with branch lengths
tree (((todMex:54.33589237,(ceyCya:37.55829174,(chlAen:35.59055317,(halSen:23.2026038,todChl:23.2026038):12.38794937):1.967738571):16.77760063):7.66410763,taeGut:62):5.4,calAnn:67.4):11.93589237

#search for 1 parameter model
lambda -s (((1,(1,(1,(1,1)1)1)1)1,1)1,1)

# generate a report
# report tmp/
