#!/usr/bin/perl -w
# MGEL
# Surya Saha 3/31/07
# reading 5 input noise .out file,
# and writing a file with the average randomly generated start pos 
# with a corresponding stop pos
# file is sorted on the start position



use strict;
use warnings;
use POSIX;

unless (@ARGV == 6){
	print "USAGE: $0 <input noise file1> <input noise file2> <input noise file3> <input noise file4> <input noise file5> <master out file>\n";
	exit;
}


#my ($ifname,$len,$rec,@temp,%temphash,@counts,@table,@famnames,@outarr,$temp,$ctr,$i,$j);
my ($orig,$ifname,$rec,@table1,@table2,@table3,@table4,@table5,@avg_table,
	$avg,@temp,$ctr,$image_len);

#get data from file 1
$ifname=$ARGV[0];
chomp $ifname;
unless(open(INFILEDATA,$ifname)){print "not able to open ".$ifname."\n\n";exit;}
while($rec=<INFILEDATA>){
	if($rec =~ /#/){next;}
	if(length ($rec) < 2){next;}
	
	push @table1, [split(' ',$rec)];
}
close (INFILEDATA);

#get data from file 2
$ifname=$ARGV[1];
chomp $ifname;
unless(open(INFILEDATA,$ifname)){print "not able to open ".$ifname."\n\n";exit;}
while($rec=<INFILEDATA>){
	if($rec =~ /#/){next;}
	if(length ($rec) < 2){next;}
	
	push @table2, [split(' ',$rec)];
}
close (INFILEDATA);

#get data from file 3
$ifname=$ARGV[2];
chomp $ifname;
unless(open(INFILEDATA,$ifname)){print "not able to open ".$ifname."\n\n";exit;}
while($rec=<INFILEDATA>){
	if($rec =~ /#/){next;}
	if(length ($rec) < 2){next;}
	
	push @table3, [split(' ',$rec)];
}
close (INFILEDATA);

#get data from file 4
$ifname=$ARGV[3];
chomp $ifname;
unless(open(INFILEDATA,$ifname)){print "not able to open ".$ifname."\n\n";exit;}
while($rec=<INFILEDATA>){
	if($rec =~ /#/){next;}
	if(length ($rec) < 2){next;}
	
	push @table4, [split(' ',$rec)];
}
close (INFILEDATA);

#get data from file 5
$ifname=$ARGV[4];
chomp $ifname;
unless(open(INFILEDATA,$ifname)){print "not able to open ".$ifname."\n\n";exit;}
while($rec=<INFILEDATA>){
	if($rec =~ /#/){next;}
	if(length ($rec) < 2){next;}
	
	push @table5, [split(' ',$rec)];
}
close (INFILEDATA);

$orig=$ARGV[5];
chomp $orig;

#initialize the avg table with one table
@avg_table=@table1;


unless(open(OUTFILEDATA,">$orig.noise.avg.tab")){print "not able to open $orig.noise.avg.tab \n\n";exit;}

#@table
#1935 10.6  0.0  2.8 chr12 8936  9225 27748096 C  R=286 Unknown (0) 283   2
#0    1     2    3   4     5     6    7        8  9     10      11  12    13

#TO AVG ALL START AND END POS AND PUT IN AVG table

for($ctr=0;$ctr<=$#avg_table;$ctr++){
	$avg=floor(($table1[$ctr][5] + $table2[$ctr][5] + $table3[$ctr][5] + $table4[$ctr][5] + $table5[$ctr][5])/5);
	$avg_table[$ctr][5] =  $avg;
	$avg=floor(($table1[$ctr][6] + $table2[$ctr][6] + $table3[$ctr][6] + $table4[$ctr][6] + $table5[$ctr][6])/5);
	$avg_table[$ctr][6] =  $avg;
}

#sorting @table on start position
@temp = sort {$a->[5] <=> $b->[5]} @avg_table;
@avg_table=@temp;

foreach(@avg_table){
	print OUTFILEDATA
	"$_->[0]\t$_->[1]\t$_->[2]\t$_->[3]\t$_->[4]\t$_->[5]\t$_->[6]\t$_->[7]\t$_->[8]\t$_->[9]\t$_->[10]\t$_->[11]\t$_->[12]\t$_->[13]\n";
}



close (OUTFILEDATA);

exit;
