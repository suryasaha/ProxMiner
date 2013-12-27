#!/usr/bin/perl -w
# MGEL
# Surya Saha 04/30/07
# reading cmd line input fasta file of families and .annot file and updated f_itemsets (with new interestinges 
# measures) file
# writing out all the f_itemsets with unannotated families 


# v1: 

use strict;
use warnings;
use POSIX;

unless (@ARGV == 3){
	print "USAGE: $0 <fasta file> <.annot file> <updated f_itemsets file>\n";
	exit;
}

my ($ifname1,$ifname2,$ifname3,$rec,@fam_names,@annot_fam_names,@unannot_fam_names,@unannot_f_itemsets_fam_names
,@f_itemsets,%temphash,@temp,$ctr,$i,$j,$flag,$user_t,$system_t,$cuser_t,$csystem_t);

$ifname1=$ARGV[0];
chomp $ifname1;
unless(open(INFILEFASTA,$ifname1)){print "not able to open $ifname1\n\n";exit;}
$ifname2=$ARGV[1];
chomp $ifname2;
unless(open(INFILEANNOT,$ifname2)){print "not able to open $ifname2\n\n";exit;}
$ifname3=$ARGV[2];
chomp $ifname3;
unless(open(INFILEFITEMSETS,$ifname3)){print "not able to open $ifname3\n\n";exit;}

unless(open(OUTFILE,">$ifname3.unannot-fams")){print "not able to open $ifname3.unannot\n\n";exit;}
unless(open(OUTFILEUNANNOT,">$ifname1.unannot-fams")){print "not able to open $ifname1.unannot-fams\n\n";exit;}

# >R=1 (RR=2.  TRF=0.005 NSEG=0.110)
# TGGGTTCTGATTAACACGAGGCCTTACGGGTACGGCGTTTGCTCCAGTGACCACAATACGATTGGAATGTCTGCTTCGTG
# TAAATTGGCGGAAATGTGGAGGTCGATTTTGCCTTGGTTGGGGACAGTTTTTGGAGTAATGCCCTTCTTCACTGCAGTTG

# getting the names of the sequences in the fasta file
$ctr=0;
while($rec=<INFILEFASTA>){
	if($rec =~ /^>/){
		$rec=~ s/>//;
		$rec=~ s/\(.*//g;
		$rec=~ s/\s*//g;
		push @fam_names,$rec;
		$ctr++;
	}
}
print "$ctr family names read\n";

# Family:R=1  Len:10670  Nof annotations:4 %annotated:99.1752577319588
# 17406 3.02 0.20 0.15 R=1 1 1984 (8686) C RIRE3A_I#LTR/Gypsy (2261) 1985 1 5
# 7065 16.11 3.38 3.66 R=1 2010 3517 (7153) TRUNCATOR#LTR/Gypsy 27 1530 (1388) 5
# 9364 10.77 0.71 1.58 R=1 3581 4985 (5685) TRUNCATOR#LTR/Gypsy 1523 2915 (3) 5
# 48930 4.39 0.18 0.16 R=1 4984 10670 (0) C RIRE8B_I#LTR/Gypsy (0) 5981 294 5


# getting the names of the families in the annot file
$ctr=0;
while($rec=<INFILEANNOT>){
	if($rec =~ /^Family/){
		$rec=~ s/^Family://;
		$rec=~ s/Len:.*//g;
		$rec=~ s/\s*//g;
		push @annot_fam_names,$rec;
		$ctr++;
	}
}
print "$ctr annotations read\n";


# #fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, Occurence, Category
# R=569	449	290	R=569	449	290	1264	u3
# R=569	449	290	R=569	449	290	1253	d3
# R=569	449	290	R=569	449	290	648	Ovlap-30to70
# R=569	449	290	R=569	449	290	502	u1
# R=569	449	290	R=6	122	1105	376	u3
# R=569	449	290	R=569	449	290	353	u2
# R=6	122	1105	R=569	449	290	329	u3

# getting the names of the families in the f_itemsets file
# print "f_itemset families are:";
$ctr=0;
while($rec=<INFILEFITEMSETS>){
	if($rec =~ /^#/){next;}
	if(length ($rec) < 10){next;}#for avoiding empty lines
	
	push @f_itemsets,[split(' ',$rec)];
	$ctr++;
}
print "$ctr frequent itemsets read\n";


# getting the unannotated families
$ctr=0;
foreach $i (@fam_names){
	$flag=0;
	foreach $j (@annot_fam_names){
		if ($i eq $j){
			$flag=1;
			last;
		}
	}	
	
	if($flag==0){
		push @unannot_fam_names,$i;
		print OUTFILEUNANNOT $i,"\n";
		$ctr++;
	}
}
print "$ctr families are unannotated. This list has been written to $ifname1.unannot-fams\n";

print OUTFILE "#fam1\tfam1-cnt\tf1-avglen\tfam2\tfam2-cnt\tf2-avglen\tOccur_cnt\tOccur_sig\tNoise\tNoise_sig\tDiff\tDiff_sig\tCat\n";
# get f_itemset records that have at least 1 unannotated family
$ctr=0;
foreach $i (@f_itemsets){
	chomp $i->[0];
	chomp $i->[3];
	foreach $j (@unannot_fam_names){
		chomp $j;
		if ($j eq $i->[0] || $j eq $i->[3]){
			print OUTFILE "$i->[0]\t$i->[1]\t$i->[2]\t$i->[3]\t$i->[4]\t$i->[5]\t$i->[6]\t$i->[7]\t$i->[8]\t$i->[9]\t$i->[10]\t$i->[11]\t$i->[12]\n";
			$ctr++;
			last;#exit as soon as even 1 family found to be unannotated
		}
	}
}

print "$ctr f_itemsets have unannotated families\n";


# calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print "\# Runtime details : \n";
print "\# System time for process: $system_t\n";
# print OUTFILEDATA "\# System time for children: $csystem_t\n";
print "\# User time for process: $user_t\n";
# print OUTFILEDATA "\# User time for children: $cuser_t\n";

close (INFILEFASTA);
close (INFILEANNOT);
close (INFILEFITEMSETS);
close (OUTFILE);
close (OUTFILEUNANNOT);

exit;
