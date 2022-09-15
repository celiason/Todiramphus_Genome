#!cafe

#specify data file, p-value threshold, # of threads to use, and log file
load -i kingfishers/Orthogroups.GeneCount.tsv -p 0.01 -t 4 -l kingfishers/reports/log_run2.txt

#the phylogenetic tree structure with branch lengths
tree (((todMex:54.33589237,(ceyCya:37.55829174,(chlAen:35.59055317,(halSen:23.2026038,todChl:23.2026038):12.38794937):1.967738571):16.77760063):7.66410763,taeGut:62):5.4,calAnn:67.4):11.93589237

# the best error model from fitting
errormodel -model kingfishers/cafe_errormodel_0.068603515625.txt -all

#search for 1 parameter model
lambda -s (((1,(1,(1,(1,1)1)1)1)1,1)1,1)

#
report kingfishers/reports/report_run2
