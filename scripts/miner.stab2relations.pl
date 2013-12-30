#!/usr/bin/perl -w
# MGEL
# Surya Saha 3/25/07
# reading tabbed .out file sorted on start position
# and find the relationship among families
# Relationship types: 1 (0-500 bases),2 (500-1000 bases),3 (1000-5000 bases)
#


use strict;
use warnings;

unless (@ARGV == 1){
	print "USAGE: $0 <input file> \n";
	exit;
}

#print "WARNING : Make sure the input is sorted on start pos!!\n";

my ($ifname,@relations,@temp,$rec,@table,$ctr,$i,$j);

$ifname=$ARGV[0];
chomp $ifname;
unless(open(INFILEDATA,$ifname)){print "not able to open ".$ifname."\n\n";exit;}

unless(open(OUTFILEDATA,">$ifname.relations.tab")){print "not able to open ".$ifname.".relations.tab \n\n";exit;}


#slurping in the whole tab file
while($rec=<INFILEDATA>){
	if($rec =~ /#/){next;}
	push @table, [split(' ',$rec)];
}


# finding the relationships
# @relation : fam1, fam2, category
$ctr=0;
for $i (0 .. $#table){
	my $seed=$table[$i][4];
	$j=$i+1;#look from next record
	#for category 1 (0-500 bases)
	#to avoid looking past the end of @table
	# pos - ONLY from start to start + 500
	while(($j < $#table) && ($table[$j][4] <= $seed+500) 
	&& ($table[$j][4] > $seed)){
		$relations[$ctr][0]=$table[$i][6];
		$relations[$ctr][1]=$table[$j][6];
		$relations[$ctr++][2]=1;
		$j++;
	}
	#for category 2 (500-1000 bases)
	#to avoid looking past the end
	# pos - ONLY from start + 500 to start +1000
	while(($j < $#table) && ($table[$j][4] <= $seed+1000) 
	&& ($table[$j][4] > $seed+500)){
		$relations[$ctr][0]=$table[$i][6];
		$relations[$ctr][1]=$table[$j][6];
		$relations[$ctr++][2]=2;
		$j++;
	}
	#for category 3 (1000-5000 bases)
	#to avoid looking past the end
	# pos - ONLY from start + 1000 to start + 5000
	while(($j < $#table) && ($table[$j][4] <= $seed+5000) 
	&& ($table[$j][4] > $seed+1000)){
		$relations[$ctr][0]=$table[$i][6];
		$relations[$ctr][1]=$table[$j][6];
		$relations[$ctr++][2]=3;
		$j++;
	}
}

#sorting @relations on fam1, fam2 and category
@temp = sort {($a->[0] cmp $b->[0]) or ($a->[1] cmp $b->[1]) or
($a->[2] <=> $b->[2]) } @relations;
@relations=@temp;


#printing the relations in tabbed format
for $i (0 .. $#relations){
	print OUTFILEDATA $relations[$i][0],"\t",
	$relations[$i][1],"\t",$relations[$i][2],"\n";
}


close (INFILEDATA);
close (OUTFILEDATA);

exit;
