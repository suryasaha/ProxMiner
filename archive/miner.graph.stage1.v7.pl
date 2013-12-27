#!/usr/bin/perl -w
# MGEL
# Surya Saha 04/30/07
# reading cmd line input merged f_itemsets file which is sorted on confidence 
# writing out clusters of related repeats using graphs


# v2: directed graph for U and D relations

# v3: undirected graph for Sum (overall) relations
# v3: printing out family counts and avg lengths
# v3: Printing out the edge weights for all conn. comps.

# Stage 1
# v1: Now using input from stage1 pipeline(only U and D rels), no strand
# v1: Printing the rel info (AvgDist and Std Dev)

# v2: Adding a edge for the reciprocal relation with the same weight

# v3: Adding in RMRB annotation information for each family

# v4:
#   : Fixed bug where conf, stdev and avgdist were assigned wrongly to relations
#   : when opp edges were being added
#   : Fixed bug where low conf rel from a conncomp was reported
#   : Printing out nof members, repeat family members and non repeat family members info

# v5: Fixed the bug where incorrect rel was reported since edges go in both dir's
#   : added in code to handle annot DB name

# v6: Creating a .html file with all coordinates. Printing coordinates in the conn-comps file. Getting elements
#   : from OUT file
#   : Writing out FASTA file for each cluster for every conn-comp with flanking seqs

# v7: optional summary file instead of annotation file.
#   : Creating Graphviz files
#   : Creating flanking sequence files
#   : New statistics for graph
#   : Using GFF file instead of OUT for coordinates

# TODO: To read in optional summary file instead of annotation file. 
# TODO: Getting the conseq lengths for non-annotated fams
# TODO: READ IN JOHN'S SUMMARY FILE (Repeats and CD's)
# TODO: HTML: Generate a table of components at the top to avoid searching manually
# TODO: COUNTS OF A RELATIONSHIP

#fam1, Category, fam2, Occ, Occ conf, Avg Rand Coin, Avg Rand Coin conf, Diff, Diff conf, Avg dist, Std Dev, fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen
# R=501	u2	R=468	8	0.8	0	0	8	0.8	1000 406.45 R=501	10	166	R=468	262	149
# R=810	u3	R=394	8	0.72	0	0	8	0.72	3854 746.43 R=810	11	224	R=394	15	289
# R=748	u3	R=225	7	0.7	0	0	7	0.7	634 37.22 R=748	10	241	R=225	22	457
         
# R=1010	717	6961.43	118	1154.5	34	15	1	1	140
# 
# R=1016	765	3828.25	117	568.94	17	13	1	0	131
# R=1016  131  2 100%
# 	R=1016 1022 8.40 0.00 0.00 R=1016 1 131 (0) F118#DNA/hAT 19 149 (30) 5 mips.rs.chr12.con.lib.cat 
# 	R=1016 1022 8.40 0.00 0.00 R=1016 1 131 (0) F118#DNA/hAT 19 149 (30) 5 rmrb.rs.chr12.con.lib.cat 

my $ver=7.0;

use strict;
use warnings;
use POSIX;
use Graph;
use Graph::Directed;

unless (@ARGV == 9){
	print "USAGE: $0 <input merged.f_itemsets file> <summ USE ANNOT FILE FOR NOW file or -nosumm> <chr fasta file> <Conf threshold> <GFF file> <-clusters or -noclusters> <-graphviz or nographviz> <-flankseq or -noflankseq> <range>\n";
	exit;
}

# Supporting functions
sub round2{
	my ($num);
	$num=$_[0];
	$num=$num*100;
	$num=int($num);
	$num=$num/100;
	return $num;
}
# get the complement
sub comp{
	my $DNA;
	$DNA=$_[0];
	$DNA=~ s/\s*//g;# clean it
	$DNA=~ tr/ACGTacgt/TGCAtgca/;
	return $DNA;
}

my ($ifname,$rec,@temp,@m_f_itemsets_table,@summ_table,$ctr,$flag,$clusterflag,@edges,$fullyannotated,$summflag, $partiallyannotated,$unannotated,$i,$j,$k,$l,$tot_recs,$tot_annots,%annots,$UDgraph,$UD_conf_lim,@UD_conn_comp,
$range,%size_counts,@vertices,%conn_counts,$fam,$length,$percent_annot,%rels,$graphvizflag,$flankseqflag,
@unannotated_ccomps,@fullyannotated_ccomps,@partiallyannotated_ccomps,@gff_table_pos,@gff_table_comp,$chr_seq,
$temp,$user_t,$system_t,$cuser_t,$csystem_t);

$ifname=$ARGV[0];
chomp $ifname;
unless(open(INFILEFITEMSETS,$ifname)){print "not able to open ".$ifname."\n\n";exit;}
$UD_conf_lim=$ARGV[3];
chomp $UD_conf_lim;
unless(open(OUTFILEDATA,">$ifname.$UD_conf_lim.conn-comps.stg1.tab")){print "not able to open $ifname.$UD_conf_lim.conn-comps.stg1.tab\n\n";exit;}

$ifname=$ARGV[1];
chomp $ifname;
if ($ifname eq '-nosumm'){ $summflag=0;}
else{
	$summflag=1;
	unless(open(INFILESUMM,$ifname)){print "not able to open ".$ifname."\n\n";exit;}
}

$ifname=$ARGV[2];
chomp $ifname;
unless(open(INFILESEQ,$ifname)){print "not able to open ".$ifname."\n\n";exit;}

$ifname=$ARGV[4];
chomp $ifname;
unless(open(INFILEGFF,$ifname)){print "not able to open ".$ifname."\n\n";exit;}
unless(open(OUTFILEHTML,">$ifname.ccomp-coords.html")){print "not able to open $ifname.ccomp-coords.html \n\n";exit;}

$clusterflag=$ARGV[5];
chomp $clusterflag;
if ($clusterflag eq "-clusters"){ $clusterflag=1;}
elsif ($clusterflag eq "-noclusters"){ $clusterflag=0;}

$graphvizflag=$ARGV[6];
chomp $graphvizflag;
if ($graphvizflag eq "-graphviz"){ $graphvizflag=1;}
elsif ($graphvizflag eq "-nographviz"){ $graphvizflag=0;}

$flankseqflag=$ARGV[7];
chomp $flankseqflag;
if ($flankseqflag eq "-flankseq"){ $flankseqflag=1;}
elsif ($flankseqflag eq "-noflankseq"){ $flankseqflag=0;}

$range=$ARGV[8];
chomp $range;

print OUTFILEHTML '<html><head><title>TIGR Rice Gbrowse links </title></head><br>';

print OUTFILEHTML '<center><br><h1> TIGR Rice Gbrowse links </h1><br></center><br><body><br>';

#fam1, Category, fam2, Occ, Occ conf, Avg Rand Coin, Avg Rand Coin conf, Diff, Diff conf, Avg dist, Std Dev, fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen
# R=501	U	R=468	8	0.8	0	0	8	0.8	1000 406.45 R=501	10	166	R=468	262	149
# R=810	D	R=394	8	0.72	0	0	8	0.72	3854 746.43 R=810	11	224	R=394	15	289

# R=1010	717	6961.43	118	1154.5	34	15	1	1	140
# 
# R=1016	765	3828.25	117	568.94	17	13	1	0	131
# R=1016  131  2 100%
# 	R=1016 1022 8.40 0.00 0.00 R=1016 1 131 (0) F118#DNA/hAT 19 149 (30) 5 mips.rs.chr12.con.lib.cat 
# 	R=1016 1022 8.40 0.00 0.00 R=1016 1 131 (0) F118#DNA/hAT 19 149 (30) 5 rmrb.rs.chr12.con.lib.cat 

# slurping in the whole freq itemsets file
$ctr=0;
while($rec=<INFILEFITEMSETS>){
	if($rec =~ /^#/){next;}
	if(length ($rec) < 10){next;}#for avoiding last line
	
	push @m_f_itemsets_table,[split(' ',$rec)];
	# record relation
	# rels:{fam1 rel fam2}=[Diff Conf, Avg Dist, Std Dev]
	@temp=();
	$temp[0]=$m_f_itemsets_table[$ctr][8];
	$temp[1]=$m_f_itemsets_table[$ctr][9];
	$temp[2]=$m_f_itemsets_table[$ctr][10];
	$rels{"$m_f_itemsets_table[$ctr][0] $m_f_itemsets_table[$ctr][1] $m_f_itemsets_table[$ctr][2]"}=[@temp];
	$ctr++;
}
close (INFILEFITEMSETS);
# record tot recs
$tot_recs = $ctr;

# slurping in the whole annotation file
if ($summflag){
	$ctr=0;
	while($rec=<INFILESUMM>){
		if($rec =~ /^#/){next;}
		if(length ($rec) < 10){next;}#for empty lines
		push @summ_table,[split(' ',$rec)];
		$ctr++;
	}
	close (INFILESUMM);
	# record tot recs
	$tot_annots = $ctr;
}


# slurping in the whole chromosome file
if ($flankseqflag || $clusterflag){
	while($rec=<INFILESEQ>){
		if($rec =~ /^>/){next;}
		if(length ($rec) < 10){next;}#for avoiding last line
		$rec=~ s/\s*//g;
		$chr_seq=$chr_seq.$rec;
	}
	close (INFILESEQ);
	print STDERR "Read in chromosome file, length : ",length($chr_seq),"\n";
}
# reading in the coordinates from the GFF file
# chr12	RepBase	repeat_region	760	1067	285	+	.	LINE-7_OS
# chr12	RepBase	repeat_region	8935	9222	1966	-	.	SUSU
$i=$j=0;
while($rec=<INFILEGFF>){
	if($rec =~ /^#/){next;}
	if(length ($rec) < 10){next;}#for avoiding last line
	@temp=();
	push @temp,[split(' ',$rec)];
	
	if($temp[0][6] eq "+"){
		$gff_table_pos[$i][0]=$temp[0][8];#fam name
		$gff_table_pos[$i][1]=$temp[0][3];#start
		$gff_table_pos[$i][2]=$temp[0][4];#end
		$i++;
	}
	elsif ($temp[0][6] eq '-'){
		$gff_table_comp[$j][0]=$temp[0][8];#fam name
		$gff_table_comp[$j][1]=$temp[0][3];#start
		$gff_table_comp[$j][2]=$temp[0][4];#end
		$j++;
	}
	else{
		print STDERR "\nBad record read\n"; exit;
	}
}
close (INFILEGFF);


# 1 graphs: One for U and D relstionships 
# Edges will be weighted by the confidence
$UDgraph = Graph::Directed->new;

#fam1, Category, fam2, Occ, Occ conf, Avg Rand Coin, Avg Rand Coin conf, Diff, Diff conf, Avg dist, Std Dev, fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen
# R=501	U	R=468	8	0.8	0	0	8	0.8	1000 406.45 R=501	10	166	R=468	262	149
# R=810	D	R=394	8	0.72	0	0	8	0.72	3854 746.43 R=810	11	224	R=394	15	289

# adding the edges to the graphs
foreach $i (@m_f_itemsets_table){
	if(($i->[8] > $UD_conf_lim) && ($i->[1] =~ /^U/)){ # for U relationships
		# if edge already exists, change the weight if the current one is higher
		if(!$UDgraph->has_edge($i->[0],$i->[2])) {
			# add edge for rel
			$UDgraph->add_weighted_edge($i->[0],$i->[2],$i->[8]);
			#add attributes to edge for avg dist and std dev
			$UDgraph->set_edge_attribute($i->[0],$i->[2],"AvgDist",$i->[9]);
			$UDgraph->set_edge_attribute($i->[0],$i->[2],"StDev",$i->[10]);
			$UDgraph->set_edge_attribute($i->[0],$i->[2],"Valid",1);
			
			# edge in opp. direction
			# check if fam2 U fam1 or fam1 D fam2 exists
			if (exists $rels{"$i->[2] U $i->[0]"} && exists $rels{"$i->[0] D $i->[2]"}){
				$j=$rels{"$i->[2] U $i->[0]"};
				$k=$rels{"$i->[0] D $i->[2]"};
				if(@$j[0] > @$k[0]){
					$UDgraph->add_weighted_edge($i->[2],$i->[0],@$j[0]);
					$UDgraph->set_edge_attribute($i->[2],$i->[0],"AvgDist",@$j[1]);
					$UDgraph->set_edge_attribute($i->[2],$i->[0],"StDev",@$j[2]);
					$UDgraph->set_edge_attribute($i->[2],$i->[0],"Valid",1);
				}
				else{
					$UDgraph->add_weighted_edge($i->[2],$i->[0],@$k[0]);
					$UDgraph->set_edge_attribute($i->[2],$i->[0],"AvgDist",@$k[1]);
					$UDgraph->set_edge_attribute($i->[2],$i->[0],"StDev",@$k[2]);
					$UDgraph->set_edge_attribute($i->[2],$i->[0],"Valid",1);
				}
			}
			# check if fam2 U fam1 exists
			elsif (exists $rels{"$i->[2] U $i->[0]"}){
				$j=$rels{"$i->[2] U $i->[0]"};
				$UDgraph->add_weighted_edge($i->[2],$i->[0],@$j[0]);
				$UDgraph->set_edge_attribute($i->[2],$i->[0],"AvgDist",@$j[1]);
				$UDgraph->set_edge_attribute($i->[2],$i->[0],"StDev",@$j[2]);
				$UDgraph->set_edge_attribute($i->[2],$i->[0],"Valid",1);
			}
			# check if fam1 D fam2 exists
			elsif (exists $rels{"$i->[0] D $i->[2]"}){
				$k=$rels{"$i->[0] D $i->[2]"};
				$UDgraph->add_weighted_edge($i->[2],$i->[0],@$k[0]);
				$UDgraph->set_edge_attribute($i->[2],$i->[0],"AvgDist",@$k[1]);
				$UDgraph->set_edge_attribute($i->[2],$i->[0],"StDev",@$k[2]);
				$UDgraph->set_edge_attribute($i->[2],$i->[0],"Valid",1);
			}
			else{
				$UDgraph->add_weighted_edge($i->[2],$i->[0],$i->[8]);
				$UDgraph->set_edge_attribute($i->[2],$i->[0],"AvgDist",$i->[9]);
				$UDgraph->set_edge_attribute($i->[2],$i->[0],"StDev",$i->[10]);
				$UDgraph->set_edge_attribute($i->[2],$i->[0],"Valid",0);
			}
			
			#add attributes to vertex
			$UDgraph->set_vertex_attribute($i->[0],"Count",$i->[12]);
			$UDgraph->set_vertex_attribute($i->[0],"AvgLen",$i->[13]);
			$UDgraph->set_vertex_attribute($i->[2],"Count",$i->[15]);
			$UDgraph->set_vertex_attribute($i->[2],"AvgLen",$i->[16]);
		}
	}
	elsif(($i->[8] > $UD_conf_lim) && ($i->[1] =~ /^D/)){ # for D relationships
		# if edge already exists, change the weight if the current one is higher
		if(!$UDgraph->has_edge($i->[2],$i->[0])) {
			# add edge for rel
			$UDgraph->add_weighted_edge($i->[2],$i->[0],$i->[8]);
			#add attributes to edge for avg dist and std dev
			$UDgraph->set_edge_attribute($i->[2],$i->[0],"AvgDist",$i->[9]);
			$UDgraph->set_edge_attribute($i->[2],$i->[0],"StDev",$i->[10]);
			$UDgraph->set_edge_attribute($i->[2],$i->[0],"Valid",1);
	
			# edge in opp. direction
			# check if fam2 D fam1 or fam1 U fam2 exists
			if (exists $rels{"$i->[2] D $i->[0]"} && exists $rels{"$i->[0] U $i->[2]"}){
				$j=$rels{"$i->[2] D $i->[0]"};
				$k=$rels{"$i->[0] U $i->[2]"};
				if(@$j[0] > @$k[0]){
					$UDgraph->add_weighted_edge($i->[0],$i->[2],@$j[0]);
					$UDgraph->set_edge_attribute($i->[0],$i->[2],"AvgDist",@$j[1]);
					$UDgraph->set_edge_attribute($i->[0],$i->[2],"StDev",@$j[2]);
					$UDgraph->set_edge_attribute($i->[0],$i->[2],"Valid",1);
				}
				else{
					$UDgraph->add_weighted_edge($i->[0],$i->[2],@$k[0]);
					$UDgraph->set_edge_attribute($i->[0],$i->[2],"AvgDist",@$k[1]);
					$UDgraph->set_edge_attribute($i->[0],$i->[2],"StDev",@$k[2]);
					$UDgraph->set_edge_attribute($i->[0],$i->[2],"Valid",1);
				}
			}
			# check if fam2 D fam1 exists
			elsif (exists $rels{"$i->[2] D $i->[0]"}){#fix 5
				$j=$rels{"$i->[2] D $i->[0]"};
				$UDgraph->add_weighted_edge($i->[0],$i->[2],@$j[0]);
				$UDgraph->set_edge_attribute($i->[0],$i->[2],"AvgDist",@$j[1]);
				$UDgraph->set_edge_attribute($i->[0],$i->[2],"StDev",@$j[2]);
				$UDgraph->set_edge_attribute($i->[0],$i->[2],"Valid",1);
			}
			# check if fam1 U fam2 exists
			elsif (exists $rels{"$i->[0] U $i->[2]"}){# fix 5
				$k=$rels{"$i->[0] U $i->[2]"};
				$UDgraph->add_weighted_edge($i->[0],$i->[2],@$k[0]);
				$UDgraph->set_edge_attribute($i->[0],$i->[2],"AvgDist",@$k[1]);
				$UDgraph->set_edge_attribute($i->[0],$i->[2],"StDev",@$k[2]);
				$UDgraph->set_edge_attribute($i->[0],$i->[2],"Valid",1);
			}
			else{
				$UDgraph->add_weighted_edge($i->[0],$i->[2],$i->[8]);
				$UDgraph->set_edge_attribute($i->[0],$i->[2],"AvgDist",$i->[9]);
				$UDgraph->set_edge_attribute($i->[0],$i->[2],"StDev",$i->[10]);
				$UDgraph->set_edge_attribute($i->[0],$i->[2],"Valid",0);
			}
	
			#add attributes to vertex
			$UDgraph->set_vertex_attribute($i->[0],"Count",$i->[12]);
			$UDgraph->set_vertex_attribute($i->[0],"AvgLen",$i->[13]);
			$UDgraph->set_vertex_attribute($i->[2],"Count",$i->[15]);
			$UDgraph->set_vertex_attribute($i->[2],"AvgLen",$i->[16]);
		}
	}
}

# Parsing the annotations to store annotations for each family
# in a array referenced by a hash
# R=1  10670  4 99.1752577319588
# R=1 16543 3.02 0.20 0.15 R=1 1 1984 (8686) C RIRE3A_I#LTR/Gypsy (2261) 1985 1 5
# R=1 6653 16.11 3.38 3.66 R=1 2010 3517 (7153) TRUNCATOR#LTR/Gypsy 27 1530 (1388) 5
# R=1 8859 10.77 0.71 1.58 R=1 3581 4985 (5685) TRUNCATOR#LTR/Gypsy 1523 2915 (3) 5
# R=1 46605 4.39 0.18 0.16 R=1 4984 10670 (0) C RIRE8B_I#LTR/Gypsy (0) 5981 294 5
# R=5  9188  3 99.8367435785808

# R=1010	717	6961.43	118	1154.5	34	15	1	1	140
# 
# R=1016	765	3828.25	117	568.94	17	13	1	0	131
# R=1016  131  2 100%
# 	R=1016 1022 8.40 0.00 0.00 R=1016 1 131 (0) F118#DNA/hAT 19 149 (30) 5 mips.rs.chr12.con.lib.cat 
# 	R=1016 1022 8.40 0.00 0.00 R=1016 1 131 (0) F118#DNA/hAT 19 149 (30) 5 rmrb.rs.chr12.con.lib.cat 
if ($summflag){
	$flag=0;
	foreach $i (@summ_table){
		if($flag==0){# for first record
			@temp=();
			# get fam info
			$fam=$i->[0];
			$length=$i->[1];
			$percent_annot=$i->[3];
			$flag=1;# for next iteration
			$ctr=0;# for @temp subscript
		}
		elsif ($flag == 1){
			if ($fam ne $i->[0]){
				# if the family exists in the graph
				if($UDgraph->has_vertex($fam)){
					# assign the annots to the hash
					# every family will have at least 1 annotation
					$annots{$fam}=[@temp];
					# assign attributes to the vertex
					$UDgraph->set_vertex_attribute($fam,"Length",$length);
					$UDgraph->set_vertex_attribute($fam,"Percent_annot",$percent_annot);
				}
				@temp=();
				#reset fam info and fam name
				$fam=$i->[0];
				$length=$i->[1];
				$percent_annot=$i->[3];
				# reset counter
				$ctr=0;
			}
			elsif ($fam eq $i->[0]){
				# add record to the temp array
				$temp[$ctr++]=$i
			}
		}
	}
}

# Get the connected components out
@UD_conn_comp=$UDgraph->strongly_connected_components();

# print the connected components
$i=localtime();
# calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;

print OUTFILEDATA "\# Version: $ver\n";
print OUTFILEDATA "\# Time: $i\n\n";
print OUTFILEDATA "\# Runtime details after preparing graphs and getting the conn components:\n";
print OUTFILEDATA "\# System time for process: $system_t\n";
print OUTFILEDATA "\# User time for process: $user_t\n\n";

print OUTFILEHTML "\# Version: $ver<br>\n";
print OUTFILEHTML "\# Time: $i<br><br>\n\n";
print OUTFILEHTML "\# Runtime details after preparing graphs and getting the conn components:<br>\n";
print OUTFILEHTML "\# System time for process: $system_t<br>\n";
print OUTFILEHTML "\# User time for process: $user_t<br><br>\n\n";


print OUTFILEDATA "\n**********************************************************************\n\n";

print OUTFILEDATA "\# Strongly connected components for Upstream and Downstream relationships using directed graph\n";
print OUTFILEHTML "\# Strongly connected components for Upstream and Downstream relationships using directed graph<br>\n";

print OUTFILEDATA "\# Merged f_itemsets file used : $ARGV[0]\n";
print OUTFILEHTML "\# Merged f_itemsets file used : $ARGV[0]<br>\n";
if ($summflag){
	print OUTFILEDATA "\# Summary/Annot file used : $ARGV[1]\n";
	print OUTFILEHTML "\# Summary/Annot file used : $ARGV[1]<br>\n";
	print OUTFILEDATA "\# Total annotations : $tot_annots\n";
	print OUTFILEHTML "\# Total annotations : $tot_annots\n";
}
print OUTFILEDATA "\# Chr file used : $ARGV[2]\n";
print OUTFILEHTML "\# Chr file used : $ARGV[2]<br>\n";
if ($flankseqflag || $clusterflag){
	print OUTFILEDATA "\# Chr length : ",length($chr_seq),"\n";
	print OUTFILEHTML "\# Chr length : ",length($chr_seq),"<br>\n";
}
print OUTFILEDATA "\# RM OUT file used : $ARGV[4]\n";
print OUTFILEHTML "\# RM OUT file used : $ARGV[4]<br>\n";

$i=$UDgraph->vertices;
print OUTFILEDATA "\# Number of vertices (families) :$i\n";
print OUTFILEHTML "\# Number of vertices (families) :$i<br>\n";

$i=$UDgraph->edges;
print OUTFILEDATA "\# Number of edges (relationships) :$i\n";
print OUTFILEHTML "\# Number of edges (relationships) :$i<br>\n";
# $i=$UDgraph->average_degree;
# print OUTFILEDATA "\# Average degree :$i\n";
print OUTFILEDATA "\# Upstream/Downstream Confidence threshold is: $UD_conf_lim\n\n";
print OUTFILEHTML "\# Upstream/Downstream Confidence threshold is: $UD_conf_lim<br><br>\n\n";

if ($graphvizflag){
	unlink glob "GraphVizFiles/* GraphVizFiles/.*";
	rmdir 'GraphVizFiles';
	mkdir "GraphVizFiles",0755 or warn "Cannot create GraphVizFiles directory: $!";
}

if ($clusterflag){
	unlink glob "Clusters/* Clusters/.*";
	rmdir 'Clusters';
	mkdir 'Clusters',0755 or warn "Cannot create Clusters directory: $!";
}

if ($flankseqflag){
	unlink glob "FlankingSeqs/* FlankingSeqs/.*";
	rmdir 'FlankingSeqs';
	mkdir 'FlankingSeqs',0755 or warn "Cannot create FlankingSeqs directory: $!";
}

$fullyannotated=$partiallyannotated=$unannotated=$ctr=0;
foreach $i(@UD_conn_comp){
	if ($#$i > 0){
		print OUTFILEDATA "\n\# COMPONENT :",$ctr,"\n";
		print OUTFILEHTML "<h2>COMPONENT &nbsp : &nbsp ",$ctr++,"</h2>\n"; 
		print OUTFILEDATA "\# Members :",scalar @$i,"\n";
		my ($rfam_count,$nonrfam_count); $rfam_count=$nonrfam_count=0;
		
		foreach $j (@$i){# if we use repeats and genes
			if ($j =~ /^R/){ $rfam_count++;}
			else { $nonrfam_count++;}
		}
		if ($nonrfam_count> 0){
			print OUTFILEDATA "\# Repeat families :",$rfam_count,"\n";
			print OUTFILEDATA "\# Other members :",$nonrfam_count,"\n";
		}
		
		if ($graphvizflag){
			#getting file name
			$temp=$ctr-1;
			$temp="CComp-".$temp.'_';
			foreach (@$i){$temp=$temp.'_'.$_;}
			
			unless(open(OUTFILEGVIZ,'>',"GraphVizFiles/$temp.gv")){print "not able to open  $temp.gv\n\n";exit;}
			print OUTFILEGVIZ "digraph G{\n";
		}
		
		print OUTFILEDATA "\# Families (Count, AvgLen, ConSeqLen, PercentAnnot):\n";
		
		if ($summflag){$l=0;}# ctr for annot count
		
		foreach $j (@$i){
			if ($graphvizflag){
				@edges=$UDgraph->edges_from($j);
				foreach $k (@edges){ print OUTFILEGVIZ "\"$k->[0]\"->\"$k->[1]\"\;\n"; }
				@edges=();
			}
			
			print OUTFILEDATA "\t",$j,"\t",$UDgraph->get_vertex_attribute($j,"Count"),"\t", $UDgraph->get_vertex_attribute($j,"AvgLen");
			if($UDgraph->has_vertex_attribute($j,'Length')){
				print OUTFILEDATA "\t",$UDgraph->get_vertex_attribute($j,'Length'),"\t", $UDgraph->get_vertex_attribute($j,'Percent_annot');
			}
			print OUTFILEDATA "\n";
			if ($summflag){
				# print annotations, if any
				if (exists $annots{$j}){
					$l++;
					$k=$annots{$j};
					@temp=@$k;
					print OUTFILEDATA "\#\tAnnot(Per_Divg, St, End, Rpt_Name\#Class, Strand(def:\+)), AnnotDB\n";
					#R=1 21917 5.13 0.18 0.25 R=1 1 2771 (7899) C Os1_07_2L#LTR (5427) 5151 2383 5 mips.rs.chr12.con.lib.cat 
					#R=1 18455 6.87 3.74 2.22 R=1 1 2752 (7918) spip (5841) 5495 2703 5 retrorz.rs.chr12.con.lib.cat
					foreach $k (@temp){
						if (defined ($k->[15])){ print OUTFILEDATA "\t\t$k->[2]\t$k->[6]\t$k->[7]\t$k->[10]\t$k->[9]\t$k->[15]\n";}
						else{ print OUTFILEDATA "\t\t$k->[2]\t$k->[6]\t$k->[7]\t$k->[9]\t$k->[14]\n";}
					}
				}
			}
		}
		
		if ($summflag){
			# update counters
			if($l==0){
				$unannotated++;
				if(@unannotated_ccomps==0){ $unannotated_ccomps[0]=$i;}
				else{ $unannotated_ccomps[$#unannotated_ccomps+1]=$i;}
			}
			elsif ($l>0 && $l<@$i){ 
				$partiallyannotated++;
				if(@partiallyannotated_ccomps==0){ $partiallyannotated_ccomps[0]=$i;}
				else{ $partiallyannotated_ccomps[$#partiallyannotated_ccomps+1]=$i;}
			}
			elsif ($l==@$i){ 
				$fullyannotated++;
				if(@fullyannotated_ccomps==0){ $fullyannotated_ccomps[0]=$i;}
				else{ $fullyannotated_ccomps[$#fullyannotated_ccomps+1]=$i;}
			}
		}
		
		if ($graphvizflag){
			print OUTFILEGVIZ "}\n";
		}
		
		print OUTFILEDATA "\n\# Relationships (AvgDist, StdDev):\n";
		# print edge weights 
		for($j=0;$j<$#$i;$j++){
			for($k=$j+1;$k<($#$i+1);$k++){#get wt for all vertices to the right
				if($UDgraph->has_edge($i->[$j],$i->[$k])){
					# NOTE: Presuming that all relations are Upstream
					#       which is correct with all elements sorted by the pos strand.
					# print the rel with higher conf if both edges are valid/real (fix 5)
					if ($UDgraph->get_edge_attribute($i->[$j],$i->[$k],"Valid")== 1 && $UDgraph->get_edge_attribute($i->[$k],$i->[$j],"Valid")== 1) {
						if ($UDgraph->get_edge_weight($i->[$j],$i->[$k]) > $UDgraph->get_edge_weight($i->[$k],$i->[$j])) {
							print OUTFILEDATA "\t",$i->[$j]," U(",$UDgraph->get_edge_weight($i->[$j],$i->[$k]),") ",$i->[$k]," ",$UDgraph->get_edge_attribute($i->[$j],$i->[$k],"AvgDist")," ",$UDgraph->get_edge_attribute($i->[$j],$i->[$k],"StDev"),"\n";
						}
						else {
							print OUTFILEDATA "\t",$i->[$k]," U(",$UDgraph->get_edge_weight($i->[$k],$i->[$j]),") ",$i->[$j]," ",$UDgraph->get_edge_attribute($i->[$k],$i->[$j],"AvgDist")," ",$UDgraph->get_edge_attribute($i->[$k],$i->[$j],"StDev"),"\n";
						}
					}
					elsif ($UDgraph->get_edge_attribute($i->[$j],$i->[$k],"Valid")== 1) {
						print OUTFILEDATA "\t",$i->[$j]," U(",$UDgraph->get_edge_weight($i->[$j],$i->[$k]),") ",$i->[$k]," ",$UDgraph->get_edge_attribute($i->[$j],$i->[$k],"AvgDist")," ",$UDgraph->get_edge_attribute($i->[$j],$i->[$k],"StDev"),"\n";
					}
					elsif ($UDgraph->get_edge_attribute($i->[$k],$i->[$j],"Valid")== 1) {
						print OUTFILEDATA "\t",$i->[$k]," U(",$UDgraph->get_edge_weight($i->[$k],$i->[$j]),") ",$i->[$j]," ",$UDgraph->get_edge_attribute($i->[$k],$i->[$j],"AvgDist")," ",$UDgraph->get_edge_attribute($i->[$k],$i->[$j],"StDev"),"\n";
					}
				}
			}
		}
# 		print OUTFILEDATA "\n\# Coordinates:\n";
		
		# get the pos and comp coordinates for all families in the ccomp
		my (@elems_pos,@elems_comp);# MOVE UP TOP AND REINIT FOR EACH ITERATION
		foreach $j (@$i){
			#get elements from @gff_table_pos
			foreach $k (@gff_table_pos){
				if ($k->[0] eq $j){ push @elems_pos,$k;}
			}
			
			#get elements from @gff_table_comp
			foreach $k (@gff_table_comp){
				if ($k->[0] eq $j){ push @elems_comp,$k;}
			}
		}
		# sort both tables on starting location
		@temp = sort {$a->[1] <=> $b->[1]} @elems_pos;
		@elems_pos=@temp;
		
		@temp = sort {$a->[1] <=> $b->[1]} @elems_comp;
		@elems_comp=@temp;
		
		# MOVE UP TOP AND REINIT FOR EACH ITERATION
		my ($ccomp_start,$ccomp_end,@ccomp_coords_pos,@ccomp_coords_comp,@members,$c_ctr);
		# get the ccomps on pos strand
		$j=0;
		@members=();
		# seed start
		$ccomp_start=$elems_pos[0][1];
		$ccomp_end=$elems_pos[0][2];
		$members[0]=$elems_pos[0][0];
		foreach $k (1..$#elems_pos){
			if ($ccomp_end+$range >= $elems_pos[$k][1]){
				$ccomp_end=$elems_pos[$k][2];
				$members[$#members+1]=$elems_pos[$k][0];
			}
			else{
				# record ccomp coords
				# only record coordinates for ccomps with >1 member
				# should this be a % of the number of families in the ccomp
				if($#members>0){
					$ccomp_coords_pos[$j][0]=$ccomp_start;
					$ccomp_coords_pos[$j][1]=$ccomp_end;
					$ccomp_coords_pos[$j++][2]=[@members];
				}
				
				# intialize for new ccomp 
				$ccomp_start=$elems_pos[$k][1];
				$ccomp_end=$elems_pos[$k][2];
				@members=(); $members[0]=$elems_pos[$k][0];
			}
		}
		# the last ccomp
		# only record coordinates for ccomps with >1 member
		# SHOULD THIS BE A % OF THE NUMBER OF FAMILIES IN THE CCOMP
		if($#members>0){
			$ccomp_coords_pos[$j][0]=$ccomp_start;
			$ccomp_coords_pos[$j][1]=$ccomp_end;
			$ccomp_coords_pos[$j][2]=[@members];
		}
		
# 		print OUTFILEDATA "Families :"; foreach (@fams){ print OUTFILEDATA $_,' ';} print OUTFILEDATA "\n";
		print OUTFILEHTML "<h2>Families &nbsp : &nbsp "; 
		foreach (@$i){ print OUTFILEHTML $_,"&nbsp ";} print OUTFILEHTML "</h2>";
		
		#getting file name
		$temp=$ctr-1;
		$temp="CComp-".$temp.'_';
		foreach (@$i){$temp=$temp.'_'.$_;}
		
		if ($clusterflag){
			unless(open(OUTFILECLUSTER,'>',"Clusters/$temp.fas")){print "not able to open  $temp.fas\n\n";exit;}
			unless(open(OUTFILECLUSTER1K,'>',"Clusters/$temp.1k.fas")){print "not able to open $temp.1k.fas\n\n";exit;}
			unless(open(OUTFILECLUSTER2K,'>',"Clusters/$temp.2k.fas")){print "not able to open $temp.2k.fas\n\n";exit;}
			unless(open(OUTFILECLUSTER4K,'>',"Clusters/$temp.4k.fas")){print "not able to open $temp.4k.fas\n\n";exit;}
			unless(open(OUTFILECLUSTER5K,'>',"Clusters/$temp.5K.fas")){print "not able to open $temp.5k.fas\n\n";exit;}
		}
		
		if ($flankseqflag){
			unless(open(OUTFILEFLANK5,'>',"FlankingSeqs/$temp.5prime.fas")){print "not able to open  $temp.5prime.fas\n\n";exit;}
			unless(open(OUTFILEFLANK3,'>',"FlankingSeqs/$temp.3prime.fas")){print "not able to open  $temp.3prime.fas\n\n";exit;}
		}
		
		# printing clusters from the pos strand
		$c_ctr=1;
# 		print OUTFILEDATA "Coordinates and members are listed from 5' to 3' on pos strand\n";
		print OUTFILEHTML "\n<h3>Coordinates &nbsp and &nbsp members &nbsp are &nbsp listed &nbsp from &nbsp 5' &nbsp to &nbsp 3' &nbsp on &nbsp pos &nbsp strand</h3>";
		foreach $j (@ccomp_coords_pos){
# 			print OUTFILEDATA 'Cluster ',$c_ctr,' '; 
			print OUTFILEHTML '<h4>Cluster &nbsp',$c_ctr++,'</h4>';
			$k=$j->[2]; 
# 			print OUTFILEDATA "Members(", scalar @$k,"): "; 
			print OUTFILEHTML "\nMembers(", scalar @$k,"): "; 
			foreach (@$k){ 
# 				print OUTFILEDATA "$_ "; 
				print OUTFILEHTML "$_ &nbsp";
			} 
# 			print OUTFILEDATA "\n"; 
			print OUTFILEHTML '<br>';
			print OUTFILEHTML 'Length : &nbsp',($j->[1]-$j->[0]),' bp <br>';
# 			print OUTFILEDATA "Exact\t\tStart:$j->[0]\tEnd:$j->[1]\n";
			print OUTFILEHTML "\nExact &nbsp &nbsp &nbsp &nbsp &nbsp Start:$j->[0] &nbsp End:$j->[1] &nbsp \t",'<a href="http://rice.plantbiology.msu.edu/cgi-bin/gbrowse/rice?name=Chr12:',$j->[0],'..',$j->[1],'">link</a><br>';
			
			# with 0 flanking
			print OUTFILEHTML "\n1k flanking &nbsp Start:",$j->[0]-1000,"&nbsp End:",$j->[1]+1000,"&nbsp &nbsp\t",'<a href="http://rice.plantbiology.msu.edu/cgi-bin/gbrowse/rice?name=Chr12:',$j->[0]-1000,'..',$j->[1]+1000,'">link</a><br>';
			
# 			print OUTFILEDATA "2k flanking\tStart:",$j->[0]-2000,"\tEnd:",$j->[1]+2000,"\n";
			print OUTFILEHTML "\n2k flanking &nbsp Start:",$j->[0]-2000,"&nbsp End:",$j->[1]+2000,"&nbsp &nbsp\t",'<a href="http://rice.plantbiology.msu.edu/cgi-bin/gbrowse/rice?name=Chr12:',$j->[0]-2000,'..',$j->[1]+2000,'">link</a><br>';
			
# 			print OUTFILEDATA "4k flanking\tStart:",$j->[0]-4000,"\tEnd:",$j->[1]+4000,"\n";
			print OUTFILEHTML "\n4k flanking &nbsp Start:",$j->[0]-4000,"&nbsp End:",$j->[1]+4000,"&nbsp &nbsp\t",'<a href="http://rice.plantbiology.msu.edu/cgi-bin/gbrowse/rice?name=Chr12:',$j->[0]-4000,'..',$j->[1]+4000,'">link</a><br>';
			
# 			print OUTFILEDATA "5k flanking\tStart:",$j->[0]-5000,"\tEnd:",$j->[1]+5000,"\n";
			print OUTFILEHTML "\n5k flanking &nbsp Start:",$j->[0]-5000,"&nbsp End:",$j->[1]+5000,"&nbsp &nbsp\t",'<a href="http://rice.plantbiology.msu.edu/cgi-bin/gbrowse/rice?name=Chr12:',$j->[0]-5000,'..',$j->[1]+5000,'">link</a><br>';
			
			if ($flankseqflag){
				print OUTFILEFLANK5 ">Cluster_pos_flanking500:",($c_ctr-1),"_Start:",$j->[0]-500,"_End:",$j->[0],"\n",substr($chr_seq,($j->[0]-501),500),"\n";
				print OUTFILEFLANK3 ">Cluster_pos_flanking500:",($c_ctr-1),"_Start:",$j->[1],"_End:",$j->[1]+500,"\n",substr($chr_seq,$j->[1],500),"\n";
			}
			
			if ($clusterflag){
				print OUTFILECLUSTER ">Cluster_pos:",($c_ctr-1),"_Start:",$j->[0],"_End:",$j->[1],"_Len:",(($j->[1]-$j->[0])+1),"\n",substr($chr_seq,($j->[0]-1),($j->[1] - $j->[0])),"\n";
				print OUTFILECLUSTER1K ">Cluster_pos:",($c_ctr-1),"_Start:",$j->[0]-1000,"_End:",$j->[1]+1000,"_Len:",(($j->[1]-$j->[0])+2001),"\n",substr($chr_seq,(($j->[0]-1)-1000),(($j->[1] - $j->[0])+2000)),"\n";
				print OUTFILECLUSTER2K ">Cluster_pos:",($c_ctr-1),"_Start:",$j->[0]-2000,"_End:",$j->[1]+2000,"_Len:",(($j->[1]-$j->[0])+4001),"\n",substr($chr_seq,(($j->[0]-1)-2000),(($j->[1] - $j->[0])+4000)),"\n";
				print OUTFILECLUSTER4K ">Cluster_pos:",($c_ctr-1),"_Start:",$j->[0]-4000,"_End:",$j->[1]+4000,"_Len:",(($j->[1]-$j->[0])+8001),"\n",substr($chr_seq,(($j->[0]-1)-4000),(($j->[1] - $j->[0])+8000)),"\n";
				print OUTFILECLUSTER5K ">Cluster_pos:",($c_ctr-1),"_Start:",$j->[0]-5000,"_End:",$j->[1]+5000,"_Len:",(($j->[1]-$j->[0])+10001),"\n",substr($chr_seq,(($j->[0]-1)-5000),(($j->[1] - $j->[0])+10000)),"\n";
			}
# 			print OUTFILEDATA "\n";
			print OUTFILEHTML "\n";
		}
		
		$j=0;
		@members=();
		# seed start
		$ccomp_start=$elems_comp[0][1];
		$ccomp_end=$elems_comp[0][2];
		$members[0]=$elems_comp[0][0];
		
		foreach $k (1..$#elems_comp){
			if ($ccomp_end+$range >= $elems_comp[$k][1]){
				$ccomp_end=$elems_comp[$k][2];
				$members[$#members+1]=$elems_comp[$k][0];
			}
			else{
				# record ccomp coords
				# only record coordinates for ccomps with >1 member
				# should this be a % of the number of families in the ccomp
				if($#members>0){
					$ccomp_coords_comp[$j][0]=$ccomp_start;
					$ccomp_coords_comp[$j][1]=$ccomp_end;
					$ccomp_coords_comp[$j++][2]=[@members];
				}
				
				# intialize for new ccomp 
				$ccomp_start=$elems_comp[$k][1];
				$ccomp_end=$elems_comp[$k][2];
				@members=(); $members[0]=$elems_comp[$k][0];
			}
		}
		# the last ccomp
		# only record coordinates for ccomps with >1 member
		# should this be a % of the number of families in the ccomp
		if($#members>0){
			$ccomp_coords_comp[$j][0]=$ccomp_start;
			$ccomp_coords_comp[$j][1]=$ccomp_end;
			$ccomp_coords_comp[$j][2]=[@members];
		}

		# printing clusters from the comp strand
		$c_ctr=1;
# 		print OUTFILEDATA "Coordinates and members are listed from 3' to 5' on comp strand\n";
		print OUTFILEHTML "\n<h3>Coordinates &nbsp and &nbsp members &nbsp are &nbsp listed &nbsp from &nbsp 3' &nbsp to &nbsp 5' &nbsp on &nbsp comp &nbsp strand</h3>";
		foreach $j (@ccomp_coords_comp){
# 			print OUTFILEDATA 'Cluster ',$c_ctr,' '; 
			print OUTFILEHTML '<h4>Cluster &nbsp',$c_ctr++,'</h4>';
			$k=$j->[2]; 
# 			print OUTFILEDATA "Members(", scalar @$k,"): "; 
			print OUTFILEHTML "\nMembers(", scalar @$k,"): "; 
			foreach (@$k){ 
# 				print OUTFILEDATA "$_ "; 
				print OUTFILEHTML "$_ &nbsp";
			}
# 			print OUTFILEDATA "\n"; 
			print OUTFILEHTML '<br>';
			print OUTFILEHTML 'Length : &nbsp',($j->[1]-$j->[0]),' bp <br>';
# 			
# 			print OUTFILEDATA "Exact\t\tStart:$j->[0]\tEnd:$j->[1]\n";
			print OUTFILEHTML "\nExact &nbsp &nbsp &nbsp &nbsp &nbsp Start:$j->[0] &nbsp End:$j->[1] &nbsp \t",'<a href="http://rice.plantbiology.msu.edu/cgi-bin/gbrowse/rice?name=Chr12:',$j->[0],'..',$j->[1],'">link</a><br>';
						
# 			print OUTFILEDATA "1k flanking\tStart:",$j->[0]-1000,"\tEnd:",$j->[1]+1000,"\n";
			print OUTFILEHTML "\n1k flanking &nbsp Start:",$j->[0]-1000,"&nbsp End:",$j->[1]+1000,"&nbsp &nbsp\t",'<a href="http://rice.plantbiology.msu.edu/cgi-bin/gbrowse/rice?name=Chr12:',$j->[0]-1000,'..',$j->[1]+1000,'">link</a><br>';
			
# 			print OUTFILEDATA "2k flanking\tStart:",$j->[0]-2000,"\tEnd:",$j->[1]+2000,"\n";
			print OUTFILEHTML "\n2k flanking &nbsp Start:",$j->[0]-2000,"&nbsp End:",$j->[1]+2000,"&nbsp &nbsp\t",'<a href="http://rice.plantbiology.msu.edu/cgi-bin/gbrowse/rice?name=Chr12:',$j->[0]-2000,'..',$j->[1]+2000,'">link</a><br>';
			
# 			print OUTFILEDATA "4k flanking\tStart:",$j->[0]-4000,"\tEnd:",$j->[1]+4000,"\n";
			print OUTFILEHTML "\n4k flanking &nbsp Start:",$j->[0]-4000,"&nbsp End:",$j->[1]+4000,"&nbsp &nbsp\t",'<a href="http://rice.plantbiology.msu.edu/cgi-bin/gbrowse/rice?name=Chr12:',$j->[0]-4000,'..',$j->[1]+4000,'">link</a><br>';
			
# 			print OUTFILEDATA "5k flanking\tStart:",$j->[0]-5000,"\tEnd:",$j->[1]+5000,"\n";
			print OUTFILEHTML "\n5k flanking &nbsp Start:",$j->[0]-5000,"&nbsp End:",$j->[1]+5000,"&nbsp &nbsp\t",'<a href="http://rice.plantbiology.msu.edu/cgi-bin/gbrowse/rice?name=Chr12:',$j->[0]-5000,'..',$j->[1]+5000,'">link</a><br>';
			
			if ($flankseqflag){
				print OUTFILEFLANK5 ">Cluster_comp_flanking500:",($c_ctr-1),"_Start:",$j->[0]-500,"_End:",$j->[0],"\n",&comp(substr($chr_seq,($j->[0]-501),500)),"\n";
				print OUTFILEFLANK3 ">Cluster_comp_flanking500:",($c_ctr-1),"_Start:",$j->[1],"_End:",$j->[1]+500,"\n",&comp(substr($chr_seq,$j->[1],500)),"\n";
			}

			if ($clusterflag){
				print OUTFILECLUSTER ">Cluster_comp:",($c_ctr-1),"_Start:",$j->[0],"_End:",$j->[1],"_Len:",(($j->[1]-$j->[0])+1),"\n",&comp(substr($chr_seq,($j->[0]-1),($j->[1] - $j->[0]))),"\n";
				print OUTFILECLUSTER1K ">Cluster_comp:",($c_ctr-1),"_Start:",$j->[0]-1000,"_End:",$j->[1]+1000,"_Len:",(($j->[1]-$j->[0])+2001),"\n",&comp(substr($chr_seq,(($j->[0]-1)-1000),(($j->[1] - $j->[0])+2000))),"\n";
				print OUTFILECLUSTER2K ">Cluster_comp:",($c_ctr-1),"_Start:",$j->[0]-2000,"_End:",$j->[1]+2000,"_Len:",(($j->[1]-$j->[0])+4001),"\n",&comp(substr($chr_seq,(($j->[0]-1)-2000),(($j->[1] - $j->[0])+4000))),"\n";
				print OUTFILECLUSTER4K ">Cluster_comp:",($c_ctr-1),"_Start:",$j->[0]-4000,"_End:",$j->[1]+4000,"_Len:",(($j->[1]-$j->[0])+8001),"\n",&comp(substr($chr_seq,(($j->[0]-1)-4000),(($j->[1] - $j->[0])+8000))),"\n";
				print OUTFILECLUSTER5K ">Cluster_comp:",($c_ctr-1),"_Start:",$j->[0]-5000,"_End:",$j->[1]+5000,"_Len:",(($j->[1]-$j->[0])+10001),"\n",&comp(substr($chr_seq,(($j->[0]-1)-5000),(($j->[1] - $j->[0])+10000))),"\n";
			}
# 			print OUTFILEDATA "\n";
			print OUTFILEHTML "\n";
		}
		print OUTFILEDATA "********************************************\n\n";
		
		if ($flankseqflag){
			close (OUTFILEFLANK5);
			close (OUTFILEFLANK3);
		}
		
		if ($clusterflag){
			close(OUTFILECLUSTER);
			close(OUTFILECLUSTER1K);
			close(OUTFILECLUSTER2K);
			close(OUTFILECLUSTER4K);
			close(OUTFILECLUSTER5K);
		}
		if ($flankseqflag){
			close (OUTFILEGVIZ);
		}
	}
	print STDERR '.';
}
print STDERR "\n";

print OUTFILEHTML '</body><br></html>';
close (OUTFILEHTML);

print OUTFILEDATA "\n**********************************************************************\n\n";

print OUTFILEDATA "\n\# Graph statistics\n";

print OUTFILEDATA "\n\# Number of ccomps: $ctr\n";
print OUTFILEDATA "\# Number of ccomps with only annotated families: $fullyannotated\n";
print OUTFILEDATA "\# Number of ccomps with some annotated families: $partiallyannotated\n";
print OUTFILEDATA "\# Number of ccomps with no annotated families: $unannotated\n\n";

# getting the number of connections at each vertex
# keeping track in hash %conn_counts
@vertices=$UDgraph->vertices();
# get all counts
foreach $i (@vertices){
	@temp=$UDgraph->edges_at($i);
	$i= scalar @temp;
	if (!exists $conn_counts{"$i"}) {
		$conn_counts{"$i"} = 1;
	}
	else{
		$conn_counts{"$i"}++;
	}
}

# print to file, sorted by conn values, ascending
# Note that values first and keys later
print OUTFILEDATA "\#For inter family relationships only\n";
print OUTFILEDATA "\#NofConns\tNofNodes\n";
foreach $i (sort {$a <=> $b} (keys (%conn_counts))){
	if($i%2==0){print OUTFILEDATA "$i\t$conn_counts{$i}\n";}
}

print OUTFILEDATA "\n\#For intra family relationships only\n";
print OUTFILEDATA "\#NofConns\tNofNodes\n";
foreach $i (sort {$a <=> $b} (keys (%conn_counts))){
	if($i%2!=0){print OUTFILEDATA "$i\t$conn_counts{$i}\n";}
}

# getting the size of counts of sizes of ccomps
foreach $i(@UD_conn_comp){
	$j= scalar @$i;
	if (!exists $size_counts{"$j"}) {
		$size_counts{"$j"} = 1;
	}
	else{
		$size_counts{"$j"}++;
	}
}

print OUTFILEDATA "\n\#For all connected components\n";
print OUTFILEDATA "\#Size\tNofCComps\n";
foreach $i (sort {$a <=> $b} (keys (%size_counts))){
	if($i>1){ # to skip ccomps of size 1, nodes with self edges
		print OUTFILEDATA "$i\t$size_counts{$i}\n";
	}
}

if ($summflag){
	# getting the size of counts of sizes of ccomps for unannotated ccomps
	%size_counts=();
	foreach $i(@unannotated_ccomps){
		$j= scalar @$i;
		if (!exists $size_counts{"$j"}) {
			$size_counts{"$j"} = 1;
		}
		else{
			$size_counts{"$j"}++;
		}
	}

	print OUTFILEDATA "\n\#For unannotated ccomps\n";
	print OUTFILEDATA "\#Size\tNofCComps\n";
	foreach $i (sort {$a <=> $b} (keys (%size_counts))){
		if($i>1){ # to skip ccomps of size 1, nodes with self edges
			print OUTFILEDATA "$i\t$size_counts{$i}\n";
		}
	}
	
	# getting the size of counts of sizes of ccomps for fully annotated ccomps
	%size_counts=();
	foreach $i(@fullyannotated_ccomps){
		$j= scalar @$i;
		if (!exists $size_counts{"$j"}) {
			$size_counts{"$j"} = 1;
		}
		else{
			$size_counts{"$j"}++;
		}
	}
	
	print OUTFILEDATA "\n\#For fully annotated ccomps\n";
	print OUTFILEDATA "\#Size\tNofCComps\n";
	foreach $i (sort {$a <=> $b} (keys (%size_counts))){
		if($i>1){ # to skip ccomps of size 1, nodes with self edges
			print OUTFILEDATA "$i\t$size_counts{$i}\n";
		}
	}
	
	# getting the size of counts of sizes of ccomps for partially annotated ccomps
	%size_counts=();
	foreach $i(@partiallyannotated_ccomps){
		$j= scalar @$i;
		if (!exists $size_counts{"$j"}) {
			$size_counts{"$j"} = 1;
		}
		else{
			$size_counts{"$j"}++;
		}
	}
	
	print OUTFILEDATA "\n\#For partially annotated ccomps\n";
	print OUTFILEDATA "\#Size\tNofCComps\n";
	foreach $i (sort {$a <=> $b} (keys (%size_counts))){
		if($i>1){ # to skip ccomps of size 1, nodes with self edges
			print OUTFILEDATA "$i\t$size_counts{$i}\n";
		}
	}
}

print OUTFILEDATA "\n\n**********************************************************************\n\n";
print OUTFILEDATA "\# Runtime details : \n";
print OUTFILEDATA "\# System time for process: $system_t\n";
print OUTFILEDATA "\# User time for process: $user_t\n";

close (OUTFILEDATA);


# calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print "\# Runtime details : \n";
print "\# System time for process: $system_t\n";
print "\# User time for process: $user_t\n";


exit;
