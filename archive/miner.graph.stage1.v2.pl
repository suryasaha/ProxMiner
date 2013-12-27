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

#fam1, Category, fam2, Occ, Occ conf, Avg Rand Coin, Avg Rand Coin conf, Diff, Diff conf, Avg dist, Std Dev, fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen
# R=501	u2	R=468	8	0.8	0	0	8	0.8	1000 406.45 R=501	10	166	R=468	262	149
# R=810	u3	R=394	8	0.72	0	0	8	0.72	3854 746.43 R=810	11	224	R=394	15	289
# R=748	u3	R=225	7	0.7	0	0	7	0.7	634 37.22 R=748	10	241	R=225	22	457
           
my $ver=2.0;

use strict;
use warnings;
use POSIX;
use Graph;
use Graph::Directed;


unless (@ARGV == 2){
	print "USAGE: $0 <input merged.f_itemsets file> <U & D conf threshold> \n";
	exit;
}

my ($ifname,$rec,@temp,@m_f_itemsets_table,$ctr,$i,$j,$k,$tot_recs,$UDgraph,
$UD_conf_lim,@UD_conn_comp,@vertices,%conn_counts,$user_t,$system_t,$cuser_t,$csystem_t);

$ifname=$ARGV[0];
chomp $ifname;
$UD_conf_lim=$ARGV[1];
chomp $UD_conf_lim;

unless(open(INFILEFITEMSETS,$ifname)){print "not able to open ".$ifname."\n\n";exit;}

unless(open(OUTFILEDATA,">$ifname.$UD_conf_lim.connected-comps.stage1.tab")){print "not able to open $ifname.$UD_conf_lim.connected-comps.tab\n\n";exit;}

# slurping in the whole freq itemsets file
$ctr=0;
while($rec=<INFILEFITEMSETS>){
	if($rec =~ /^#/){next;}
	if(length ($rec) < 10){next;}#for avoiding last line
	
	push @m_f_itemsets_table,[split(' ',$rec)];
	$ctr++;
}
close (INFILEFITEMSETS);

# record tot recs
$tot_recs = $ctr;

# 1 graphs: One for U and D relstionships 
# Edges will be weighted by the confidence
$UDgraph = Graph::Directed->new;


#fam1, Category, fam2, Occ, Occ conf, Avg Rand Coin, Avg Rand Coin conf, Diff, Diff conf, Avg dist, Std Dev, fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen
# R=501	U	R=468	8	0.8	0	0	8	0.8	1000 406.45 R=501	10	166	R=468	262	149
# R=810	D	R=394	8	0.72	0	0	8	0.72	3854 746.43 R=810	11	224	R=394	15	289
# Add counts as attributes for each vertex

# adding the edges to the graphs
foreach $i (@m_f_itemsets_table){
	if(($i->[8] > $UD_conf_lim) && ($i->[1] =~ /^U/)){ # for U relationships
		# if edge already exists, change the weight if the current one is higher
		if($UDgraph->has_edge($i->[0],$i->[2]) && ($UDgraph->get_edge_weight($i->[0],$i->[2]) < $i->[8])) {
			$UDgraph->set_edge_weight($i->[0],$i->[2],$i->[8]);
		}
		else{# else create a new edge
			$UDgraph->add_weighted_edge($i->[0],$i->[2],$i->[8]);
		}
		# edge in opp. direction
		if($UDgraph->has_edge($i->[2],$i->[0]) && ($UDgraph->get_edge_weight($i->[2],$i->[0]) < $i->[8])) {
			$UDgraph->set_edge_weight($i->[2],$i->[0],$i->[8]);
		}
		else{# else create a new edge
			$UDgraph->add_weighted_edge($i->[2],$i->[0],$i->[8]);
		}

		#add attributes to vertex
		$UDgraph->set_vertex_attribute($i->[0],"Count",$i->[12]);
		$UDgraph->set_vertex_attribute($i->[0],"AvgLen",$i->[13]);
		$UDgraph->set_vertex_attribute($i->[2],"Count",$i->[15]);
		$UDgraph->set_vertex_attribute($i->[2],"AvgLen",$i->[16]);
		#add attributes to edge for avg dist and std dev
		$UDgraph->set_edge_attribute($i->[0],$i->[2],"AvgDist",$i->[9]);
		$UDgraph->set_edge_attribute($i->[0],$i->[2],"StDev",$i->[10]);
		$UDgraph->set_edge_attribute($i->[2],$i->[0],"AvgDist",$i->[9]);
		$UDgraph->set_edge_attribute($i->[2],$i->[0],"StDev",$i->[10]);
	}
	elsif(($i->[8] > $UD_conf_lim) && ($i->[1] =~ /^D/)){ # for D relationships
		# if edge already exists, change the weight if the current one is higher
		if($UDgraph->has_edge($i->[2],$i->[0]) && ($UDgraph->get_edge_weight($i->[2],$i->[0]) < $i->[8])) {
			$UDgraph->set_edge_weight($i->[2],$i->[0],$i->[8]);
		}
		else{# else create a new edge
			$UDgraph->add_weighted_edge($i->[2],$i->[0],$i->[8]);
		}
		# edge in opp direction
		if($UDgraph->has_edge($i->[0],$i->[2]) && ($UDgraph->get_edge_weight($i->[0],$i->[2]) < $i->[8])) {
			$UDgraph->set_edge_weight($i->[0],$i->[2],$i->[8]);
		}
		else{# else create a new edge
			$UDgraph->add_weighted_edge($i->[0],$i->[2],$i->[8]);
		}

		#add attributes to vertex
		$UDgraph->set_vertex_attribute($i->[0],"Count",$i->[12]);
		$UDgraph->set_vertex_attribute($i->[0],"AvgLen",$i->[13]);
		$UDgraph->set_vertex_attribute($i->[2],"Count",$i->[15]);
		$UDgraph->set_vertex_attribute($i->[2],"AvgLen",$i->[16]);
		#add attributes to edge for avg dist and std dev
		$UDgraph->set_edge_attribute($i->[2],$i->[0],"AvgDist",$i->[9]);
		$UDgraph->set_edge_attribute($i->[2],$i->[0],"StDev",$i->[10]);
		$UDgraph->set_edge_attribute($i->[0],$i->[2],"AvgDist",$i->[9]);
		$UDgraph->set_edge_attribute($i->[0],$i->[2],"StDev",$i->[10]);
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
		print OUTFILEDATA "\n\# Component :",$ctr++,"\n";
		print OUTFILEDATA "\# Families (Count, AvgLen):\n";
		foreach $j (@$i){
			print OUTFILEDATA "\t",$j,"\t",$UDgraph->get_vertex_attribute($j,"Count"),"\t", $UDgraph->get_vertex_attribute($j,"AvgLen"),"\n";

		}
		print OUTFILEDATA "\# Relationships (AvgDist, StdDev):\n";
		# print edge weights 
		for($j=0;$j<$#$i;$j++){
			for($k=$j+1;$k<($#$i+1);$k++){#get wt for all vertices to the right
				if($UDgraph->has_edge($i->[$j],$i->[$k])){
					print OUTFILEDATA "\t",$i->[$j]," U(",$UDgraph->get_edge_weight($i->[$j],$i->[$k]),") ",$i->[$k]," ",$UDgraph->get_edge_attribute($i->[$j],$i->[$k],"AvgDist")," ",$UDgraph->get_edge_attribute($i->[$j],$i->[$k],"StDev"),"\n";
				}
			}
		}
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
