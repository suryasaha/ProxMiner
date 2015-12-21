#!/usr/bin/perl -w
# MGEL
# Surya Saha 3/30/07
# reading cmd line input RM .out file,
# and writing a dummy hit details to output file
# the ONLY change is the randomly generated start pos with a corresponding stop pos
# Output is sorted on the start position



use strict;
use warnings;
use POSIX;

unless (@ARGV == 3){
	print "USAGE: $0 <input .out file> <noise file index> <length of chromosome>\n";
	exit;
}


#my ($ifname,$len,$rec,@temp,%temphash,@counts,@table,@famnames,@outarr,$temp,$ctr,$i,$j);
my ($ifname,$chr_len,$rec,@table,$index,@temp,$ctr,$image_len);

$ifname=$ARGV[0];
chomp $ifname;
$index=$ARGV[1];
chomp $index;
$chr_len=$ARGV[2];
chomp $chr_len;


unless(open(INFILEDATA,$ifname)){print "not able to open ".$ifname."\n\n";exit;}
unless(open(OUTFILEDATA,">$ifname.noise.$index")){print "not able to open ".$ifname.".noise.$index \n\n";exit;}

#@table
#1935 10.6  0.0  2.8 chr12 8936  9225 27748096 C  R=286 Unknown (0) 283   2
#0    1     2    3   4     5     6    7        8  9     10      11  12    13

#slurping in the whole report file
#adding random start and end positions 
#and writing it out

srand();

#skipping header lines in new .OUT format
for (1..3) { $rec=<INFILEDATA> }
while($rec=<INFILEDATA>){
	if(length ($rec) < 2){next;}
	
	push @table, [split(' ',$rec)];
	$image_len=$table[$#table][6] - $table[$#table][5];
	#generating random start point
	$table[$#table][5] = int(rand($chr_len - $image_len)) + 1;
	#generating corresponding end point
	$table[$#table][6] = $table[$#table][5] + $image_len;
}

#sorting @table on start position
@temp = sort {$a->[5] <=> $b->[5]} @table;
@table=@temp;

foreach(@table){
	print OUTFILEDATA "$_->[0]\t$_->[1]\t$_->[2]\t$_->[3]\t$_->[4]\t$_->[5]\t$_->[6]\t$_->[7]\t$_->[8]\t$_->[9]\t$_->[10]\t$_->[11]\t$_->[12]\t$_->[13]\n";
}

close (INFILEDATA);
close (OUTFILEDATA);

exit;
