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

# TODO: ADD IN MORE ANNOTATIONS FOR FAMILIES, LENGTHS FOR NON-ANNOTATED FAMILIES

#fam1, Category, fam2, Occ, Occ conf, Avg Rand Coin, Avg Rand Coin conf, Diff, Diff conf, Avg dist, Std Dev, fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen
# R=501	u2	R=468	8	0.8	0	0	8	0.8	1000 406.45 R=501	10	166	R=468	262	149
# R=810	u3	R=394	8	0.72	0	0	8	0.72	3854 746.43 R=810	11	224	R=394	15	289
# R=748	u3	R=225	7	0.7	0	0	7	0.7	634 37.22 R=748	10	241	R=225	22	457
         
my $ver=4.0;

use strict;
use warnings;
use POSIX;
use Graph;
use Graph::Directed;


unless (@ARGV == 3){
	print "USAGE: $0 <input merged.f_itemsets file> <.annot file> <U & D conf threshold> \n";
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

my ($ifname,$rec,@temp,@m_f_itemsets_table,@rmrb_annots_table,$ctr,$flag, $i,$j,$k,$tot_recs,$tot_rmrb_annots,%rmrb_annots,$UDgraph,$UD_conf_lim,@UD_conn_comp,
@vertices,%conn_counts,$fam,$length,$rmrb_percent_annot,%rels,
$user_t,$system_t,$cuser_t,$csystem_t);

$ifname=$ARGV[0];
chomp $ifname;
unless(open(INFILEFITEMSETS,$ifname)){print "not able to open ".$ifname."\n\n";exit;}
$UD_conf_lim=$ARGV[2];
chomp $UD_conf_lim;
unless(open(OUTFILEDATA,">$ifname.$UD_conf_lim.conn-comps.stg1.tab")){print "not able to open $ifname.$UD_conf_lim.conn-comps.stg1.tab\n\n";exit;}

$ifname=$ARGV[1];
chomp $ifname;
unless(open(INFILEANNOT,$ifname)){print "not able to open ".$ifname."\n\n";exit;}

#fam1, Category, fam2, Occ, Occ conf, Avg Rand Coin, Avg Rand Coin conf, Diff, Diff conf, Avg dist, Std Dev, fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen
# R=501	U	R=468	8	0.8	0	0	8	0.8	1000 406.45 R=501	10	166	R=468	262	149
# R=810	D	R=394	8	0.72	0	0	8	0.72	3854 746.43 R=810	11	224	R=394	15	289

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
$ctr=0;
while($rec=<INFILEANNOT>){
	if($rec =~ /^#/){next;}
	if(length ($rec) < 10){next;}#for empty lines
	push @rmrb_annots_table,[split(' ',$rec)];
	$ctr++;
}
close (INFILEANNOT);
# record tot recs
$tot_rmrb_annots = $ctr;


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
			
			# edge in opp. direction
			# check if fam2 U fam1 or fam1 D fam2 exists
			if (exists $rels{"$i->[2] U $i->[0]"} && exists $rels{"$i->[0] D $i->[2]"}){
				$j=$rels{"$i->[2] U $i->[0]"};
				$k=$rels{"$i->[0] D $i->[2]"};
				if(@$j[0] > @$k[0]){
					$UDgraph->add_weighted_edge($i->[2],$i->[0],@$j[0]);
					$UDgraph->set_edge_attribute($i->[2],$i->[0],"AvgDist",@$j[1]);
					$UDgraph->set_edge_attribute($i->[2],$i->[0],"StDev",@$j[2]);
				}
				else{
					$UDgraph->add_weighted_edge($i->[2],$i->[0],@$k[0]);
					$UDgraph->set_edge_attribute($i->[2],$i->[0],"AvgDist",@$k[1]);
					$UDgraph->set_edge_attribute($i->[2],$i->[0],"StDev",@$k[2]);
				}
			}
			elsif (exists $rels{"$i->[2] U $i->[0]"}){
				$j=$rels{"$i->[2] U $i->[0]"};
				$UDgraph->add_weighted_edge($i->[2],$i->[0],@$j[0]);
				$UDgraph->set_edge_attribute($i->[2],$i->[0],"AvgDist",@$j[1]);
				$UDgraph->set_edge_attribute($i->[2],$i->[0],"StDev",@$j[2]);
			}
			elsif (exists $rels{"$i->[0] D $i->[2]"}){
				$k=$rels{"$i->[0] D $i->[2]"};
				$UDgraph->add_weighted_edge($i->[2],$i->[0],@$k[0]);
				$UDgraph->set_edge_attribute($i->[2],$i->[0],"AvgDist",@$k[1]);
				$UDgraph->set_edge_attribute($i->[2],$i->[0],"StDev",@$k[2]);
			}
			else{
				$UDgraph->add_weighted_edge($i->[2],$i->[0],$i->[8]);
				$UDgraph->set_edge_attribute($i->[2],$i->[0],"AvgDist",$i->[9]);
				$UDgraph->set_edge_attribute($i->[2],$i->[0],"StDev",$i->[10]);
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
	
			# edge in opp. direction
			# check if fam2 D fam1 or fam1 U fam2 exists
			if (exists $rels{"$i->[2] D $i->[0]"} && exists $rels{"$i->[0] U $i->[2]"}){
				$j=$rels{"$i->[2] D $i->[0]"};
				$k=$rels{"$i->[0] U $i->[2]"};
				if(@$j[0] > @$k[0]){
					$UDgraph->add_weighted_edge($i->[0],$i->[2],@$j[0]);
					$UDgraph->set_edge_attribute($i->[0],$i->[2],"AvgDist",@$j[1]);
					$UDgraph->set_edge_attribute($i->[0],$i->[2],"StDev",@$j[2]);
				}
				else{
					$UDgraph->add_weighted_edge($i->[0],$i->[2],@$k[0]);
					$UDgraph->set_edge_attribute($i->[0],$i->[2],"AvgDist",@$k[1]);
					$UDgraph->set_edge_attribute($i->[0],$i->[2],"StDev",@$k[2]);
				}
			}
			elsif (exists $rels{"$i->[2] U $i->[0]"}){
				$j=$rels{"$i->[2] D $i->[0]"};
				$UDgraph->add_weighted_edge($i->[0],$i->[2],@$j[0]);
				$UDgraph->set_edge_attribute($i->[0],$i->[2],"AvgDist",@$j[1]);
				$UDgraph->set_edge_attribute($i->[0],$i->[2],"StDev",@$j[2]);
			}
			elsif (exists $rels{"$i->[0] D $i->[2]"}){
				$k=$rels{"$i->[0] U $i->[2]"};
				$UDgraph->add_weighted_edge($i->[0],$i->[2],@$k[0]);
				$UDgraph->set_edge_attribute($i->[0],$i->[2],"AvgDist",@$k[1]);
				$UDgraph->set_edge_attribute($i->[0],$i->[2],"StDev",@$k[2]);
			}
			else{
				$UDgraph->add_weighted_edge($i->[0],$i->[2],$i->[8]);
				$UDgraph->set_edge_attribute($i->[0],$i->[2],"AvgDist",$i->[9]);
				$UDgraph->set_edge_attribute($i->[0],$i->[2],"StDev",$i->[10]);
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
$flag=0;
foreach $i (@rmrb_annots_table){
	if($flag==0){# for first record
		@temp=();
		# get fam info
		$fam=$i->[0];
		$length=$i->[1];
		$rmrb_percent_annot=$i->[3];
		$flag=1;# for next iteration
		$ctr=0;# for @temp subscript
	}
	elsif ($flag == 1){
		if ($fam ne $i->[0]){
			# if the family exists in the graph
			if($UDgraph->has_vertex($fam)){
				# assign the annots to the hash
				# every family will have at least 1 annotation
				$rmrb_annots{$fam}=[@temp];
				# assign attributes to the vertex
				$UDgraph->set_vertex_attribute($fam,"Length",$length);
				$UDgraph->set_vertex_attribute($fam,"Percent_annot",&round2($rmrb_percent_annot));
			}
			@temp=();
			#reset fam info and fam name
			$fam=$i->[0];
			$length=$i->[1];
			$rmrb_percent_annot=$i->[3];
			# reset counter
			$ctr=0;
		}
		elsif ($fam eq $i->[0]){
			# add record to the temp array
			$temp[$ctr++]=$i
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


print OUTFILEDATA "\n**********************************************************************\n\n";

print OUTFILEDATA "\# Strongly connected components for Upstream and Downstream relationships using directed graph\n";
print OUTFILEDATA "\# Annotation file used : $ifname\n";
$i=$UDgraph->vertices;
print OUTFILEDATA "\# Number of vertices (families) :$i\n";
$i=$UDgraph->edges;
print OUTFILEDATA "\# Number of edges (relationships) :$i\n";
# $i=$UDgraph->average_degree;
# print OUTFILEDATA "\# Average degree :$i\n";
print OUTFILEDATA "\# Upstream/Downstream Confidence threshold is: $UD_conf_lim\n\n";

$ctr=0;
foreach $i(@UD_conn_comp){
	if ($#$i > 0){
		print OUTFILEDATA "\n\# COMPONENT :",$ctr++,"\n";
		print OUTFILEDATA "\# Members :",scalar @$i,"\n";
		my ($rfam_count,$nonrfam_count); $rfam_count=$nonrfam_count=0;
		foreach $j (@$i){
			if ($j =~ /^R/){ $rfam_count++;}
			else { $nonrfam_count++;}
		}
		if ($nonrfam_count> 0){
			print OUTFILEDATA "\# Repeat families :",$rfam_count,"\n";
			print OUTFILEDATA "\# Other members :",$nonrfam_count,"\n";
		}
		print OUTFILEDATA "\# Families (Count, AvgLen, ConSeqLen, PercentAnnot):\n";
		foreach $j (@$i){
			print OUTFILEDATA "\t",$j,"\t",$UDgraph->get_vertex_attribute($j,"Count"),"\t", $UDgraph->get_vertex_attribute($j,"AvgLen");
			if($UDgraph->has_vertex_attribute($j,'Length')){
				print OUTFILEDATA "\t",$UDgraph->get_vertex_attribute($j,'Length'),"\t", $UDgraph->get_vertex_attribute($j,'Percent_annot'),"\%";
			}
			print OUTFILEDATA "\n";
			# print annotations, if any
			if (exists $rmrb_annots{$j}){
				$k=$rmrb_annots{$j};
				@temp=@$k;
				print OUTFILEDATA "\#\tRMRB_Annot(Per_Divg, St, End, Rpt_Name\#Class, Strand(def:\+))\n";
				# R=1 16543 3.02 0.20 0.15 R=1 1 1984 (8686) C RIRE3A_I#LTR/Gypsy (2261) 1985 1 5
				# R=1 6653 16.11 3.38 3.66 R=1 2010 3517 (7153) TRUNCATOR#LTR/Gypsy 27 1530 (1388) 5
				foreach $k (@temp){
					if (defined ($k->[14])){ print OUTFILEDATA "\t\t$k->[2]\t$k->[6]\t$k->[7]\t$k->[10]\t$k->[9]\n";}
					else{ print OUTFILEDATA "\t\t$k->[2]\t$k->[6]\t$k->[7]\t$k->[9]\n";}
				}
			}
		}
		
		print OUTFILEDATA "\# Relationships (AvgDist, StdDev):\n";
		# print edge weights 
		for($j=0;$j<$#$i;$j++){
			for($k=$j+1;$k<($#$i+1);$k++){#get wt for all vertices to the right
				if($UDgraph->has_edge($i->[$j],$i->[$k])){
					# NOTE: Presuming that all relations are Upstream
					#       which is correct with all elements sorted by the pos strand.
					# print the rel with higher conf
					if ($UDgraph->get_edge_weight($i->[$j],$i->[$k]) > $UDgraph->get_edge_weight($i->[$k],$i->[$j])) {
						print OUTFILEDATA "\t",$i->[$j]," U(",$UDgraph->get_edge_weight($i->[$j],$i->[$k]),") ",$i->[$k]," ",$UDgraph->get_edge_attribute($i->[$j],$i->[$k],"AvgDist")," ",$UDgraph->get_edge_attribute($i->[$j],$i->[$k],"StDev"),"\n";
					}
					else {
						print OUTFILEDATA "\t",$i->[$k]," U(",$UDgraph->get_edge_weight($i->[$k],$i->[$j]),") ",$i->[$j]," ",$UDgraph->get_edge_attribute($i->[$k],$i->[$j],"AvgDist")," ",$UDgraph->get_edge_attribute($i->[$k],$i->[$j],"StDev"),"\n";
					}
				}
			}
		}
		print OUTFILEDATA "********************************************\n";
	}
}

print OUTFILEDATA "\n**********************************************************************\n\n";

print OUTFILEDATA "\n\# Graph statistics\n";

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
print OUTFILEDATA "\#Nodes\tConn's\n";
foreach $i (sort {$a <=> $b} (keys (%conn_counts))){
	print OUTFILEDATA "$conn_counts{$i}\t$i\n";
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
