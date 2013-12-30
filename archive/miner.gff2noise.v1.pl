#!/usr/bin/perl -w
# MGEL
# Surya Saha 3/30/07
# reading cmd line input .gff file,
# and writing a dummy hit details to output file
# the ONLY change is the randomly generated start pos with a corresponding stop pos
# Output is sorted on the start position



use strict;
use warnings;
use POSIX;

unless (@ARGV == 3){
	print "USAGE: $0 <input .gff file> <noise file index> <length of chromosome>\n";
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
# chr1	UCSC-galGal3	NM_001080714	13751	20858	.	+	.	gene
#  0         1               2           3          4       5    6    7     8  

#slurping in the whole report file
#adding random start and end positions 
#and writing it out

srand();
while($rec=<INFILEDATA>){
	if($rec =~ /#/){next;}
	if(length ($rec) < 2){next;}
	
	push @table, [split(' ',$rec)];
	$image_len=$table[$#table][4] - $table[$#table][3];
	#generating random start point
	$table[$#table][3] = int(rand($chr_len - $image_len)) + 1;
	#generating corresponding end point
	$table[$#table][4] = $table[$#table][3] + $image_len;
}

#sorting @table on start position
@temp = sort {$a->[3] <=> $b->[3]} @table;
@table=@temp;

#@table
# chr1	UCSC-galGal3	NM_001080714	13751	20858	.	+	.	gene
#  0         1               2           3          4       5    6    7     8  

foreach(@table){
	print OUTFILEDATA "$_->[0]\t$_->[1]\t$_->[2]\t$_->[3]\t$_->[4]\t$_->[5]\t$_->[6]\t$_->[7]\t$_->[8]\n";
}



close (INFILEDATA);
close (OUTFILEDATA);

exit;
