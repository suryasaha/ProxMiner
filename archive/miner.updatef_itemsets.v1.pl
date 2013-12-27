#!/usr/bin/perl -w
# MGEL
# Surya Saha 04/30/07
# reading cmd line input f_itemsets file which is sorted on the count 
# writing out a new f_itemsets file with new interestingness measures

# v1: 

use strict;
use warnings;
use POSIX;

unless (@ARGV == 1){
	print "USAGE: $0 <input .f_itemsets.tab file>\n";
	exit;
}

my ($ifname1,$rec,@temp,@f_itemsets_table,$ctr,$i,$j,$int_msr,$tot_recs,
$user_t,$system_t,$cuser_t,$csystem_t);

$ifname1=$ARGV[0];
chomp $ifname1;
unless(open(INFILEFITEMSETS,$ifname1)){print "not able to open ".$ifname1."\n\n";exit;}

unless(open(OUTFILEDATA,">$ifname1.updated")){print "not able to open $ifname1.updated\n\n";exit;}

# slurping in the whole freq itemsets file
$ctr=0;
while($rec=<INFILEFITEMSETS>){
	if($rec =~ /^#/){next;}
	if(length ($rec) < 10){next;}#for avoiding last line
	
	push @f_itemsets_table,[split(' ',$rec)];
	$ctr++;
}
# record tot recs
$tot_recs = $ctr;

# # now sort it by f1, f2 and category
# # R=133	334	493	R=133	334	493	696	u3
# @temp = sort {($a->[0] cmp $b->[0]) or ($a->[3] cmp $b->[3]) or
# ($a->[7] cmp $b->[7]) } @f_itemsets_table;
# @f_itemsets_table=@temp;


# 0       1          2             3     4           5            6          7       8          9
# #fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, Occurence, Noise, Difference, Category
# R=569	449	290	R=569	449	290	1253	0	1253	d3
# R=569	449	290	R=569	449	290	1264	29	1235	u3
# R=569	449	290	R=569	449	290	648	0	648	Ovlap-30to70
# R=569	449	290	R=569	449	290	502	1	501	u1

print OUTFILEDATA "#fam1\tfam1-cnt\tf1-avglen\tfam2\tfam2-cnt\tf2-avglen\tOccur_cnt\tOccur_sig\tNoise\tNoise_sig\tDiff\t Diff_sig\tCat\n";
foreach $i (@f_itemsets_table){
	print OUTFILEDATA "$i->[0]\t$i->[1]\t$i->[2]\t$i->[3]\t$i->[4]\t$i->[5]\t$i->[6]\t";
	
	if($i->[1] > $i->[4]){
		$int_msr=$i->[6]/$i->[4];
	}
	elsif($i->[4] > $i->[1]){
		$int_msr=$i->[6]/$i->[1];
	}
	else{#both are equal, same family relation??
		$int_msr=$i->[6]/$i->[1];
	}
	
	#chop off to 2 decimal places
	$int_msr=$int_msr*100;
	$int_msr=int($int_msr);
	$int_msr=$int_msr/100;
	
	print OUTFILEDATA "$int_msr\t$i->[7]\t";

	if ($i->[7] !=0){
		if($i->[1] > $i->[4]){
			$int_msr=$i->[7]/$i->[4];
		}
		elsif($i->[4] > $i->[1]){
			$int_msr=$i->[7]/$i->[1];
		}
		else{#both are equal, same family relation??
			$int_msr=$i->[7]/$i->[1];
		}
		#chop off to 2 decimal places
		$int_msr=$int_msr*100;
		$int_msr=int($int_msr);
		$int_msr=$int_msr/100;
	}
	else{
		$int_msr=0;
	}
	
	print OUTFILEDATA "$int_msr\t$i->[8]\t";

	if($i->[1] > $i->[4]){
		$int_msr=$i->[8]/$i->[4];
	}
	elsif($i->[4] > $i->[1]){
		$int_msr=$i->[8]/$i->[1];
	}
	else{#both are equal, same family relation??
		$int_msr=$i->[8]/$i->[1];
	}
	
	#chop off to 2 decimal places
	$int_msr=$int_msr*100;
	$int_msr=int($int_msr);
	$int_msr=$int_msr/100;
	
	print OUTFILEDATA "$int_msr\t$i->[9]\n";
}

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
