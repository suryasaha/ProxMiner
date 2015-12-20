#!/usr/bin/perl -w
# MGEL
# Surya Saha 4/24/07
# reading cmd line input .out file which is sorted on the start position
# and translate to GFF file
# GFF v 2.0 fields are: 
# <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]

# v1: 

use strict;
use warnings;
use POSIX;

unless (@ARGV == 1){
	print "USAGE: $0 <input .out file> \n";
	exit;
}

my ($ifname,$rec,@temp,@table,$tot_recs,$i,$ctr,$user_t,$system_t,$cuser_t,$csystem_t);

$ifname=$ARGV[0];
chomp $ifname;
unless(open(INFILEDATA,$ifname)){print "not able to open ".$ifname."\n\n";exit;}
unless(open(OUTFILEDATA,">$ifname.gff")){print "not able to open $ifname.gff \n\n";exit;}

print OUTFILEDATA "\##gff-version 2\n";
print OUTFILEDATA "\# Assuming repeat finder is : RptSct-1.0.1\n";


#slurping in the whole report file
$ctr=0;
while($rec=<INFILEDATA>){
	if($rec =~ /#/){next;}
	if(length ($rec) < 10){next;}#for avoiding last line
	push @table, [split(' ',$rec)];
	$ctr++;
}
# record tot recs
$tot_recs = $ctr;
print OUTFILEDATA "\# Total records read: $tot_recs\n";
print OUTFILEDATA "\# <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]\n";
#for $i (0 .. $#table){
#	for $j (0 .. $#{$table[$i]}){
#		print "record $i $j is $table[$i][$j]\n";
#	}
#}

#@table
#1935 10.6  0.0  2.8 chr12 8936  9225 27748096 C  R=286 Unknown (0) 283   2
#0    1     2    3   4     5     6    7        8  9     10      11  12    13

#sorting @table on family name
@temp = sort {$a->[9] cmp $b->[9]} @table;
@table=@temp;

# GFF fields are: <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
# printing to GFF format
foreach $i (@table){
	print OUTFILEDATA "$i->[4]\tRptSct-1.0.1\trepeat_unit\t$i->[5]\t$i->[6]\t$i->[0]\t";
	if ($i->[8] eq 'C'){
		print OUTFILEDATA "\-\t.\tRepeat family $i->[9]\; Sub type $i->[10]\n";
	}
	else{
		print OUTFILEDATA "\+\t.\tRepeat family $i->[9]\; Sub type $i->[10]\n";
	}

}


#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print  OUTFILEDATA "\# Runtime details after printing: \n";
print  OUTFILEDATA "\# System time for process: $system_t\n";
#print OUTFILEDATA "\# System time for children: $csystem_t\n";
print  OUTFILEDATA "\# User time for process: $user_t\n";
#print OUTFILEDATA "\# User time for children: $cuser_t\n";


close (INFILEDATA);
close (OUTFILEDATA);

exit;
