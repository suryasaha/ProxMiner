ProxMiner
=========
PIPELINE
1. miner.gff/out2noise.v1.pl : Create 5 noise OUT files using the real OUT file of RM
2. miner.gff/out2f_itemsets.vX.pl : Reads in the OUT file of RM wt RS library and prints out the frequent itemsets and the copy/cluster information for each relationship. Do this for the real OUT fle and the 5 noise OUT files
3a. miner.merge_ALL_f_itemsets.v3.pl : Merge all the f_itemset files to create merged f_itemsets file and a stats file
3b. miner.merge_ALL_f_itemsets.stage1.v1.multi.pl : Merge the original f_itemset and any number of noise f_itemset files to create merged f_itemsets file and a stats file
4. miner.graph.stage1.v5.pl : Produce connected components using the merged f_itemsets file and RMRB annotations.

STAND ALONE SCRIPTS
1. miner.out2gff.v2.pl : Reads in .out file which is sorted on the start position and translates it to GFF file for TIGR genome browser
2. miner.fam-members.v1.pl : Writes out mfasta file for each family with all its members
3. miner.graph.v5.pl : Reads in merged f_itemsets file which is sorted on confidence and writes out clusters of related repeats using graphs
4. miner.libannotator.v3.pl : Reads in RM .cat file for the novel repeat library and finds the repeat families that are annotated for > X% of their total length and prints out their annotation in .annot file
5. miner.unannotated-fams.v1.pl : Reads in fasta file of families , .annot file and f_itemsets file and writing out all the f_itemsets with unannotated families 
6. replicate.sh : Creates any number of noise files and finds the relationships
7. miner.ucsctable2gff.v1.pl : converts ucsc table to gff2 format