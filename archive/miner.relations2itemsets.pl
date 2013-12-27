#!/usr/bin/perl -w
# MGEL
# Surya Saha 3/26/07
# reading tabbed relationship file sorted on fam1
# and creating the itemsets in the format
# fam1 fam2 Occurence Category



use strict;
use warnings;

unless (@ARGV == 1){
	print "USAGE: $0 <input file> \n";
	exit;
}

#print "WARNING: Make sure relations are sorted on fam1, fam2 and cat\n";

my ($ifname,@itemsets,@temp,$rec,@table,$ctr,$i,$j,
	$fam1,$fam2,$count_1,$count_2,$count_3);

$ifname=$ARGV[0];
chomp $ifname;
unless(open(INFILEDATA,$ifname)){print "not able to open ".$ifname."\n\n";exit;}

unless(open(OUTFILEDATA,">$ifname.itemsets.tab")){print "not able to open ".$ifname.".itemsets.tab \n\n";exit;}

#slurping in the whole tab file
while($rec=<INFILEDATA>){
	if($rec =~ /#/){next;}
	push @table, [split(' ',$rec)];
}

#creating the itemsets
# @itemsets : fam1 fam2 Occurence Category
# @table (was @relations): fam1 fam2 Category (sorted on fam1,fam2,category)
$ctr=0;
for $i (0 .. $#table){
	#for the first time
	if($ctr==0){
		$fam1=$table[$i][0];
		$fam2=$table[$i][1];
		$count_1=$count_2=$count_3=0;
		if ($table[$i][2] eq '1'){ $count_1++;}
		elsif ($table[$i][2] == 2){ $count_2++;}
		elsif ($table[$i][2] == 3){ $count_3++;}
	}

	#for the rest
	if ($fam1 eq $table[$i][0] && $fam2 eq $table[$i][1]){
		if ($table[$i][2] == 1){ $count_1++;}
		elsif ($table[$i][2] == 2){ $count_2++;}
		elsif ($table[$i][2] == 3){ $count_3++;}
	}
	else {
		#write itemsets for that family pair
		$itemsets[$ctr][0]=$fam1;
		$itemsets[$ctr][1]=$fam2;
		$itemsets[$ctr][2]=$count_1;
		$itemsets[$ctr++][3]=1;
		$itemsets[$ctr][0]=$fam1;
		$itemsets[$ctr][1]=$fam2;
		$itemsets[$ctr][2]=$count_2;
		$itemsets[$ctr++][3]=2;
		$itemsets[$ctr][0]=$fam1;
		$itemsets[$ctr][1]=$fam2;
		$itemsets[$ctr][2]=$count_3;
		$itemsets[$ctr++][3]=3;
		
		#initialize for next family pair
		$fam1=$table[$i][0];
		$fam2=$table[$i][1];
		$count_1=$count_2=$count_3=0;

		if ($table[$i][2] == 1){ $count_1++;}
		elsif ($table[$i][2] == 2){ $count_2++;}
		elsif ($table[$i][2] == 3){ $count_3++;}
	}
}
#write itemsets for last family pair
$itemsets[$ctr][0]=$fam1;
$itemsets[$ctr][1]=$fam2;
$itemsets[$ctr][2]=$count_1;
$itemsets[$ctr++][3]=1;
$itemsets[$ctr][0]=$fam1;
$itemsets[$ctr][1]=$fam2;
$itemsets[$ctr][2]=$count_2;
$itemsets[$ctr++][3]=2;
$itemsets[$ctr][0]=$fam1;
$itemsets[$ctr][1]=$fam2;
$itemsets[$ctr][2]=$count_3;
$itemsets[$ctr++][3]=3;

#sorting @itemsets on occurences in decreasing order
@temp = sort {$b->[2] <=> $a->[2]} @itemsets;
@itemsets=@temp;

#printing the itemsets in tabbed format
# @itemsets : fam1 fam2 Occurence Category
for $i (0 .. $#itemsets){
	if($itemsets[$i][2] > 0){
		print OUTFILEDATA $itemsets[$i][0],"\t",$itemsets[$i][1],"\t",
		$itemsets[$i][2],"\t",$itemsets[$i][3],"\n";
	}
}


close (INFILEDATA);
close (OUTFILEDATA);

exit;
