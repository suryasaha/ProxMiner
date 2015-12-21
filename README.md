<h1><b>ProxMiner</b></h1>

A association rule based and graph mining pipeline to identify frequently co-occuring repeat regions in chromosomal length genome sequences.

=========
<p>
<h4>PIPELINE</h4>
<ol>
    <li>Run RS on query sequence. Run RM with repeat library constructed by RS on query sequence.</li>
    <li>miner.gff/out2noise.v1.pl: Create 5 noise OUT files using the query OUT file of RM</li>
    <li>miner.gff/out2f_itemsets.stage1.vX.pl: Reads in the OUT file of RM wt RS library and prints out the frequent itemsets and the copy/cluster information for each relationship. Do this for the query OUT file and the 5 noise OUT files</li>
    <li>
        Merge f_itemsets from query sequence and background noise sets 
        <ul>
            <li>miner.merge_ALL_f_itemsets.stage1.v1.pl: Merge all the 5 noise and 1 real f_itemset files to create merged f_itemsets file and a stats file</li>
            <li>miner.merge_ALL_f_itemsets.stage1.v1.multi.pl: Merge the query f_itemset and any number of noise f_itemset files to create merged f_itemsets file and a stats file</li>
        </ul>
    <li>miner.graph.stage1.v8.pl: Produce connected components using the merged f_itemsets file and RMRB annotations of both the RS repeat library and the query sequence.</li>
    </li>
</ol>
</p>
<p>
<h4>STAND ALONE SCRIPTS</h4>
<ol>
    <li>miner.out2gff.v2.pl: Reads in .OUT file which is sorted on the start position and translates it to GFF file</li>
    <li>miner.fam-members.v1.pl: Writes out mfasta file for each repeat family with all its members</li>
    <li>miner.graph.v5.pl: Reads in merged f_itemsets file sorted on confidence and writes out clusters of related repeats using graphs</li>
    <li>miner.libannotator.v3.pl: Reads in RM .CAT file for the RS repeat library and finds the repeat families that are annotated for > X% of their total length and prints out their annotation in .ANNOT file</li>
    <li>miner.unannotated-fams.v1.pl: Reads in fasta file of families, .ANNOT file and f_itemsets file and writing out all the f_itemsets with unannotated families </li>
    <li>replicate.sh: Creates any number of noise files and finds the relationships</li>
    <li>miner.ucsctable2gff.v1.pl: converts ucsc table to gff2 format</li>
</ol>
</p>
<p>
<h5>Glossary</h5>
<ul>
    <li>RM: RepeatMasker http://bix.ucsd.edu/repeatscout/</li>
    <li>RS: RepeatScout http://bix.ucsd.edu/repeatscout/</li>
    <li>RMRB: RepeatMasker with RepBase library http://www.girinst.org/</li>
    <li>OUT: RepeatMasker OUT file http://www.repeatmasker.org/webrepeatmaskerhelp.html</li>
</ul>
