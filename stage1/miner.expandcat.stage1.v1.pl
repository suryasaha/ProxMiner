#!/usr/bin/perl -w
# MGEL
# Surya Saha 01/23/09
# reading the CAT file produced by RM and the fasta file of the repeat library used
# expand the repeat name in the cat file


# CAT file
# 50692 1.12 0.16 0.65 R=102 1 5587 (2601) ORSiTERT00200154 972 6531 (2855) 5
# 417 25.20 1.61 0.80 R=104 1 124 (7) C ORSiTETNOOT00148 (2315) 1093 969 5
# 16130 5.76 0.20 0.81 R=105 1 1995 (949) C ORSiTERTOOT00274 (0) 4749 2767 5
# 9361 6.14 1.71 0.08 R=110 1 1172 (0) C ORSiCMCM00100010 (8437) 1191 1 5
# 2118 5.08 0.00 0.00 R=112 1 256 (0) C ORSgTEMT03400028 (0) 256 1 5
# 90835 7.23 0.20 0.40 R=117 1 11317 (0) ORSgTERT00200317 1 11295 (0) 5



# v1: 


use strict;
use warnings;
use POSIX;

unless (@ARGV == 2){
	print "USAGE: $0 < .CAT file> < .fasta library>\n";
	exit;
}

my ($ifname,$rec,@temp,@catnames,@fastanames,%name_index,$ctr,$i,$j,$k,
$tot_repeats,$tot_records,$user_t,$system_t,$cuser_t,$csystem_t);

$ifname=$ARGV[0];
chomp $ifname;
unless(open(CATFILE,$ifname)){print "not able to open ".$ifname."\n\n";exit;}

unless(open(OUTFILE,">$ifname.fullnames")){print "not able to open ".$ifname."fullnames\n\n";exit;}

$ifname=$ARGV[1];
chomp $ifname;
unless(open(FASTAFILE,$ifname)){print "not able to open ".$ifname."\n\n";exit;}

# reading in names from the fasta file
$ctr=0;
while($rec=<FASTAFILE>){
	if($rec =~ /^#/){next;}
	if(length ($rec) < 2){next;}#for avoiding last line
	
	if($rec =~ /^>/){push @fastanames,$rec; $ctr++;}
}
# record tot recs
$tot_repeats = $ctr;
close (FASTAFILE);


# build hash of abbrvs to full names
foreach $i (@fastanames){
	$i=~ s/^>//;
	$j=$i;
	$j=~ s/ \S*//g;
	$i=~ s/ /_/g;
	$i=~ s/\s//g;
	$j=~ s/\s//g;
	$name_index{$j}=$i;# add to hash
# 	print $j,"==>",$i,"\n";
# 	print $name_index{$j};
}


# open cat file and write out complete records to fullnames file
$ctr=0;
while($rec=<CATFILE>){
	if($rec =~ /^#/){next;}
	if($rec =~ /^[rR]/){next;}
	if(length ($rec) < 2){next;}#for avoiding last line
	
	@temp=split(' ',$rec);
	print OUTFILE "$temp[0] $temp[1] $temp[2] $temp[3] $temp[4] $temp[5] $temp[6] $temp[7] ";
	# to handle matches to comp strand
	if ($temp[8] eq 'C'){
		$temp[9]=~ s/\s//g;
		print OUTFILE "$temp[8] $name_index{$temp[9]} ";
		print OUTFILE "$temp[10] $temp[11] $temp[12] $temp[13]\n"; 
	}
	else{
		$temp[8]=~ s/\s//g;
		print OUTFILE "$name_index{$temp[8]} ";
		print OUTFILE "$temp[9] $temp[10] $temp[11] $temp[12]\n"; 
	}
	$ctr++;
}
# record tot recs
$tot_records = $ctr;
close (CATFILE);
close (OUTFILE);


# calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print "\# Runtime details : \n";
print "\# System time for process: $system_t\n";
# print OUTFILEDATA "\# System time for children: $csystem_t\n";
print "\# User time for process: $user_t\n";
# print OUTFILEDATA "\# User time for children: $cuser_t\n";



exit;
