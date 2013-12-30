#!/usr/bin/perl -w
# MGEL
# Surya Saha 4/22/08
# Sorts the records in a GFF file to prepare the file for use by gffef_itemsets script

# GFF v 2.0 fields are: 
# <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
# chr1	RptSct-1.0.1	repeat_unit	9599047	9599208	606	+	.	R=825
# chr1	RptSct-1.0.1	repeat_unit	20233020	20233420	2669	+	.	R=825
# chr1	RptSct-1.0.1	repeat_unit	20233390	20233502	352	+	.	R=825
# chr1	RptSct-1.0.1	repeat_unit	20233477	20233567	237	+	.	R=825

# v1: 


use strict;
use warnings;
use POSIX;

my $ver=1.0;

unless (@ARGV == 1){
	print "USAGE: $0 <input .gff>\n";
	exit;
}

my ($ifname,$rec,@temp,@table,$tot_recs,$i,$ctr,$user_t,$system_t,$cuser_t,$csystem_t);

$ifname=$ARGV[0];
chomp $ifname;
unless(open(INFILEDATA,$ifname)){print "not able to open ".$ifname."\n\n";exit;}
unless(open(OUTFILEDATA,">$ifname.sorted")){print "not able to open $ifname.gff \n\n";exit;}

print OUTFILEDATA "\##gff-version 2\n";
print OUTFILEDATA "\# ver:$ver\n";

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
print OUTFILEDATA "\# Total records: $tot_recs\n";
print OUTFILEDATA "\# <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]\n";

#@table
# <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
# chr1	RptSct-1.0.1	repeat_unit	9599047	9599208	606	+	.	R=825


#sort ascending @table on start position
@temp = sort {$a->[3] <=> $b->[3]} @table;
@table=@temp;

# GFF fields are: <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
# printing to GFF format
foreach $i (@table){
	print OUTFILEDATA "$i->[0]\t$i->[1]\t$i->[2]\t$i->[3]\t$i->[4]\t$i->[5]\t$i->[6]\t$i->[7]\t$i->[8]\n";
}


#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print  OUTFILEDATA "\# Runtime details after printing: \n";
print  OUTFILEDATA "\# System time for process: $system_t\n";
#print OUTFILEDATA "\# System time for children: $csystem_t\n";
print  OUTFILEDATA "\# User time for process: $user_t\n";
#print OUTFILEDATA "\# User time for children: $cuser_t\n";

print  STDERR "Runtime details after printing: \n";
print  STDERR "System time for process: $system_t\n";
print  STDERR "User time for process: $user_t\n";
print STDERR "Does not print the [comments] column, if present\n";


close (INFILEDATA);
close (OUTFILEDATA);

exit;
