<h1><b>ProxMiner</b></h1>
=========
<p><h4>PIPELINE</h4>
<ol>
<li>miner.gff/out2noise.v1.pl : Create 5 noise OUT files using the real OUT file of RM</li>
<li>miner.gff/out2f_itemsets.vX.pl : Reads in the OUT file of RM wt RS library and prints out the frequent itemsets and the copy/cluster information for each relationship. Do this for the real OUT file and the 5 noise OUT files
<ul>
<li>miner.merge_ALL_f_itemsets.v3.pl : Merge all the f_itemset files to create merged f_itemsets file and a stats file</li>
<li>miner.merge_ALL_f_itemsets.stage1.v1.multi.pl : Merge the original f_itemset and any number of noise f_itemset files to create merged f_itemsets file and a stats file</li>
</ul>
<li>miner.graph.stage1.v5.pl : Produce connected components using the merged f_itemsets file and RMRB annotations.</li>
</ol>
</p>
<p><h4>STAND ALONE SCRIPTS</h4>
<ol>
<li>miner.out2gff.v2.pl : Reads in .out file which is sorted on the start position and translates it to GFF file for TIGR genome browser</li>
<li>miner.fam-members.v1.pl : Writes out mfasta file for each family with all its members</li>
<li>miner.graph.v5.pl : Reads in merged f_itemsets file which is sorted on confidence and writes out clusters of related repeats using graphs</li>
<li>miner.libannotator.v3.pl : Reads in RM .cat file for the novel repeat library and finds the repeat families that are annotated for > X% of their total length and prints out their annotation in .annot file</li>
<li>miner.unannotated-fams.v1.pl : Reads in fasta file of families , .annot file and f_itemsets file and writing out all the f_itemsets with unannotated families </li>
<li>replicate.sh : Creates any number of noise files and finds the relationships</li>
<li>miner.ucsctable2gff.v1.pl : converts ucsc table to gff2 format</li>
</ol></p>
<p><h5>Glossary</h5>
<ul>
<li>RM: RepeatMasker</li>
<li>RMRB: RepeatMasker with RepBase library</li>
<li>OUT: RepeatMasker OUT file</li>
</ul>