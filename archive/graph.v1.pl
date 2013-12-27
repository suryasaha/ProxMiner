#!/usr/bin/perl -w
# MGEL
# Surya Saha 04/30/07
# reading cmd line input merged f_itemsets file which is sorted on confidence 
# writing out a new f_itemsets file with new interestingness measures


# v1: 

# #fam1, Category, fam2, Occ, Occ conf, Avg Rand Coin, Avg Rand Coin conf, Diff, Diff conf, Strand, fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen
# R=501	u2	R=468	8	0.8	0	0	8	0.8	+	R=501	10	166	R=468	262	149
# R=810	u3	R=394	8	0.72	0	0	8	0.72	+	R=810	11	224	R=394	15	289
# R=748	u3	R=225	7	0.7	0	0	7	0.7	+	R=748	10	241	R=225	22	457
# R=851	d1	R=825	9	0.69	0	0	9	0.69	+	R=851	13	130	R=825	13	206


use strict;
use warnings;
use POSIX;
use Graph;
use Graph::Directed;
use Graph::Undirected;

unless (@ARGV == 1){
	print "USAGE: $0 <input merged.f_itemsets file>\n";
	exit;
}

my ($ifname,$rec,@temp,@m_f_itemsets_table,$ctr,$i,$j,$tot_recs,$UDgraph,$OICgraph,
@UD_conn_comp,@OIC_conn_comp,$user_t,$system_t,$cuser_t,$csystem_t);

$ifname=$ARGV[0];
chomp $ifname;
unless(open(INFILEFITEMSETS,$ifname)){print "not able to open ".$ifname."\n\n";exit;}

unless(open(OUTFILEDATA,">$ifname.connected-comps.tab")){print "not able to open $ifname.connected-comps.tab\n\n";exit;}

# slurping in the whole freq itemsets file
$ctr=0;
while($rec=<INFILEFITEMSETS>){
	if($rec =~ /^#/){next;}
	if(length ($rec) < 10){next;}#for avoiding last line
	
	push @m_f_itemsets_table,[split(' ',$rec)];
	$ctr++;
}
# record tot recs
$tot_recs = $ctr;

# 2 graphs: One for U and D relstionships and one for Ovlap, IN and Cont relationships
# Edges will be weighted by the confidence

$UDgraph = Graph::Undirected->new;
$OICgraph = Graph::Undirected->new;

# #fam1, Category, fam2, Occ, Occ conf, Avg Rand Coin, Avg Rand Coin conf, Diff, Diff conf, Strand, fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen
# R=501	u2	R=468	8	0.8	0	0	8	0.8	+	R=501	10	166	R=468	262	149

# $UDgraph->add_weighted_edge(1,2,3.0);
# $UDgraph->add_weighted_edge(2,3,1.0);
# $UDgraph->set_edge_weight(1,2,4.0);

# adding the edges to the graphs
foreach $i (@m_f_itemsets_table){
	if($i->[1] =~ /[ud]/){ # for U and D relationships
		# if edge already exists, change the weight if the current one is higher
		if($UDgraph->has_edge($i->[0],$i->[2]) && ($UDgraph->get_edge_weight($i->[0],$i->[2]) < $i->[8])) {
			$UDgraph->set_edge_weight($i->[0],$i->[2],$i->[8]);
		}
		else{# else create a new edge
			$UDgraph->add_weighted_edge($i->[0],$i->[2],$i->[8]);
		}
	}
	elsif($i->[1] =~ /[ICp]/){ # for IN, CONT and Ovlap relationships
		# if edge already exists, change the weight if the current one is higher
		if($OICgraph->has_edge($i->[0],$i->[2]) && ($OICgraph->get_edge_weight($i->[0],$i->[2]) < $i->[8])){
			$OICgraph->set_edge_weight($i->[0],$i->[2],$i->[8]);
		}
		else{# else create a new edge
			$OICgraph->add_weighted_edge($i->[0],$i->[2],$i->[8]);
		}
	}
	else{
		print STDERR "Whoops!!...";
	}
}


# Get the connected components out
@UD_conn_comp=$UDgraph->connected_components();
@OIC_conn_comp=$OICgraph->connected_components();

# print the connected components
$i=localtime();
# calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;

print OUTFILEDATA "\# Version: 1.0\n";
print OUTFILEDATA "\# Time: $i\n\n";
print OUTFILEDATA "\# Runtime details after preparing graphs and getting the conn components:\n";
print OUTFILEDATA "\# System time for process: $system_t\n";
print OUTFILEDATA "\# User time for process: $user_t\n\n";

print OUTFILEDATA "\# Connected components for Upstream and Downstream relationships\n\n";
$ctr=0;
foreach $i(@UD_conn_comp){
	if ($#$i > 0){
		print OUTFILEDATA "\# Component :",$ctr++,"\n";
		foreach $j (@$i){
			print OUTFILEDATA $j,' ';
		}
		print OUTFILEDATA "\n";
	}
}

print OUTFILEDATA "\n**********************************************************************\n\n";
print OUTFILEDATA "\# Connected components for Within, Overlap and Contains relationships\n\n";
$ctr=0;
foreach $i(@OIC_conn_comp){
	if ($#$i > 0){
		print OUTFILEDATA "\# Component :",$ctr++,"\n";
		foreach $j (@$i){
			print OUTFILEDATA $j,' ';
		}
		print OUTFILEDATA "\n";
	}
}
print OUTFILEDATA "\n**********************************************************************\n\n";

print OUTFILEDATA "\# Runtime details : \n";
print OUTFILEDATA "\# System time for process: $system_t\n";
print OUTFILEDATA "\# User time for process: $user_t\n";


# calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print "\# Runtime details : \n";
print "\# System time for process: $system_t\n";
# print OUTFILEDATA "\# System time for children: $csystem_t\n";
print "\# User time for process: $user_t\n";
# print OUTFILEDATA "\# User time for children: $cuser_t\n";

close (INFILEFITEMSETS);
close (OUTFILEDATA);



exit;
