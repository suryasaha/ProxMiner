#!/usr/bin/perl -w
# MGEL
# Surya Saha 12/24/08
# reading cmd line input merged f_itemsets file which is sorted on confidence 
# writing out number of rules with confidence above certain thresholds


#fam1, Category, fam2, Occ, Occ conf, Avg Rand Coin, Avg Rand Coin conf, Diff, Diff conf, Avg dist, Std Dev, fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen
# R=501	u2	R=468	8	0.8	0	0	8	0.8	1000 406.45 R=501	10	166	R=468	262	149
# R=810	u3	R=394	8	0.72	0	0	8	0.72	3854 746.43 R=810	11	224	R=394	15	289
# R=748	u3	R=225	7	0.7	0	0	7	0.7	634 37.22 R=748	10	241	R=225	22	457
         
my $ver=1.0;

use strict;
use warnings;
use POSIX;


unless (@ARGV == 1){
	print "USAGE: $0 <input merged.f_itemsets file>\n";
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

my ($ifname,$rec,@temp,@m_f_itemsets_table,%rel_ctrs,$ctr,$flag,$i,$j,$k,$key,$value,$tot_recs
,$user_t,$system_t,$cuser_t,$csystem_t);

$ifname=$ARGV[0];
chomp $ifname;
unless(open(INFILEFITEMSETS,$ifname)){print "not able to open ".$ifname."\n\n";exit;}
unless(open(OUTFILEDATA,">$ifname.fdr.tab")){print "not able to open $ifname.fdr.tab\n\n";exit;}


#fam1, Category, fam2, Occ, Occ conf, Avg Rand Coin, Avg Rand Coin conf, Diff, Diff conf, Avg dist, Std Dev, fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen
# R=501	U	R=468	8	0.8	0	0	8	0.8	1000 406.45 R=501	10	166	R=468	262	149
# R=810	D	R=394	8	0.72	0	0	8	0.72	3854 746.43 R=810	11	224	R=394	15	289

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

# initialize counters
for ($i=0.05;$i<=1.05;$i+=0.05){
	$rel_ctrs{"$i"}=0;
}

#fam1, Category, fam2, Occ, Occ conf, Avg Rand Coin, Avg Rand Coin conf, Diff, Diff conf, Avg dist, Std Dev, fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen
# R=501	U	R=468	8	0.8	0	0	8	0.8	1000 406.45 R=501	10	166	R=468	262	149
# R=810	D	R=394	8	0.72	0	0	8	0.72	3854 746.43 R=810	11	224	R=394	15	289

# updating the counters
foreach $i (@m_f_itemsets_table){
	while(($key, $value)= each %rel_ctrs){
# 		if (($i->[12]<=50) && ($i->[15]<=50)){ for families with elements within a range
		if($i->[4] >= $key){
			$rel_ctrs{$key}++;} #increment ctr if valid
	}
}

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
print OUTFILEDATA "\# Merged fitemsets file used : $ifname\n";

# getting the ctrs in a sorted array
print OUTFILEDATA "Confidence threshold\tNumber of rules\n";
foreach $key (sort keys %rel_ctrs){
	print OUTFILEDATA $key,"\t",$rel_ctrs{$key},"\n";
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
