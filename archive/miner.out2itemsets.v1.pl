#!/usr/bin/perl -w
# MGEL
# Surya Saha 3/15/07
# reading cmd line input RM .out file,
# and writing hit details to tabbed output file
# which is sorted on the start position
# and find the relationship among families
# Relationship types: 1 (0-500 bases),2 (500-1000 bases),3 (1000-5000 bases)
# and creating the itemsets in the format
# fam1 fam2 Occurence Category

#RUNtime: 35 mins

use strict;
use warnings;
use POSIX;

unless (@ARGV == 1){
	print "USAGE: $0 <input file> \n";
	exit;
}


my ($ifname,$rec,@temp,@table,@famnames,@relations,@itemsets,@counts,%temphash,
$ctr,$i,$j,$fam1,$fam2,$count_1,$count_2,$count_3);

$ifname=$ARGV[0];
chomp $ifname;
unless(open(INFILEDATA,$ifname)){print "not able to open ".$ifname."\n\n";exit;}

unless(open(OUTFILEDATA,">$ifname.itemsets.tab")){print "not able to open ".$ifname.".tab \n\n";exit;}

#to get the count of a family
#params: $fam 
sub get_count{
	my ($fam);
	$fam=$_[0];
	$fam=~ s/\s*//g;
	
	foreach (@counts){
		$_->[0] =~ s/\s*//g;
		if ($_->[0] eq $fam){ 
		return $_->[1];
		last;
		}
	}
}

#to get the avg length of a family
#params: $fam
sub get_avglen{
	my ($fam);
	$fam=$_[0];
	$fam=~ s/\s*//g;

	foreach (@counts){
		$_->[0] =~ s/\s*//g;
		if ($_->[0] eq $fam){ 
		return $_->[2];
		last;
		}
	}
}

#slurping in the whole report file
while($rec=<INFILEDATA>){
	if($rec =~ /#/){next;}
	if(length ($rec) < 10){next;}#for avoiding last line
	push @table, [split(' ',$rec)];
}


#for $i (0 .. $#table){
#	for $j (0 .. $#{$table[$i]}){
#		print "record $i $j is $table[$i][$j]\n";
#	}
#}

#@table
#1935 10.6  0.0  2.8 chr12 8936  9225 27748096 C  R=286 Unknown (0) 283   2
#0    1     2    3   4     5     6    7        8  9     10      11  12    13

#sorting @table on start position
@temp = sort {$a->[5] <=> $b->[5]} @table;
@table=@temp;


#find the number of occurences of each family
#get family names
$ctr=0;
foreach(@table){
	$famnames[$ctr++]=$_->[9];
}

#removing duplicates
#sorting
@temp=sort @famnames;
@famnames=@temp;
%temphash = map { $_, 1 } @famnames;
@famnames = keys %temphash;
#sorting agn
@temp=sort @famnames;
@famnames=@temp;


#initializing the @counts 2D array
# @count: fam occurences avg-len
$ctr=0;
foreach(@famnames){
	$counts[$ctr][0]=$_;
	#initializing all counters to 0
	$counts[$ctr][1]=0;#occurences
	$counts[$ctr++][2]=0;#avg length
}

#count the number of times a family is found and its avg length
foreach $i (@counts){
	foreach $j (@table){
		if($i->[0] eq $j->[9]){
			$i->[1]++;#occurences
			$i->[2] = $i->[2] + ($j->[6] - $j->[5]);#length
		}
	}
	$i->[2]=floor($i->[2] / $i->[1]);#avg length
}

#testing
# print "\nPrinting \@counts...............\n";
# foreach(@counts){ print $_->[0],' ',$_->[1],"\t",$_->[2],"\t";}
# print "\n\n For R=997:",&get_count("R=997")," and ",&get_avglen("R=997"),"\n";


#finding the relationships
#@relations : fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, category
$ctr=0;
for $i (0 .. $#table){
	my $seed=$table[$i][5];
	$j=$i+1;#look from next record
	#for category 1 (0-500 bases)
	#to avoid looking past the end of @table
	# pos - ONLY from start to start + 500
	while(($j < $#table) && ($table[$j][5] <= $seed+500) 
	&& ($table[$j][5] > $seed)){
		$relations[$ctr][0]=$table[$i][9];
		$relations[$ctr][1]=&get_count($table[$i][9],@counts);
		$relations[$ctr][2]=&get_avglen($table[$i][9],@counts);
		$relations[$ctr][3]=$table[$j][9];
		$relations[$ctr][4]=&get_count($table[$j][9],@counts);
		$relations[$ctr][5]=&get_avglen($table[$j][9],@counts);
		$relations[$ctr++][6]=1;
		$j++;
	}
	#for category 2 (500-1000 bases)
	#to avoid looking past the end
	# pos - ONLY from start + 500 to start +1000
	while(($j < $#table) && ($table[$j][5] <= $seed+1000) 
	&& ($table[$j][5] > $seed+500)){
		$relations[$ctr][0]=$table[$i][9];
		$relations[$ctr][1]=&get_count($table[$i][9],@counts);
		$relations[$ctr][2]=&get_avglen($table[$i][9],@counts);
		$relations[$ctr][3]=$table[$j][9];
		$relations[$ctr][4]=&get_count($table[$j][9],@counts);
		$relations[$ctr][5]=&get_avglen($table[$j][9],@counts);
		$relations[$ctr++][6]=2;
		$j++;
	}
	#for category 3 (1000-5000 bases)
	#to avoid looking past the end
	# pos - ONLY from start + 1000 to start + 5000
	while(($j < $#table) && ($table[$j][5] <= $seed+5000) 
	&& ($table[$j][5] > $seed+1000)){
		$relations[$ctr][0]=$table[$i][9];
		$relations[$ctr][1]=&get_count($table[$i][9],@counts);
		$relations[$ctr][2]=&get_avglen($table[$i][9],@counts);
		$relations[$ctr][3]=$table[$j][9];
		$relations[$ctr][4]=&get_count($table[$j][9],@counts);
		$relations[$ctr][5]=&get_avglen($table[$j][9],@counts);
		$relations[$ctr++][6]=3;
		$j++;
	}
}



#sorting @relations on fam1, fam2 and category
#@relations : fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, category
@temp = sort {($a->[0] cmp $b->[0]) or ($a->[3] cmp $b->[3]) or
($a->[6] <=> $b->[6]) } @relations;
@relations=@temp;

#testing
# print "\nPrinting \@relations...............\n";
# foreach(@relations){ print $_->[0],' ',$_->[1],"\t",$_->[2],"\t",$_->[3],"\t",$_->[4],"\t",
# $_->[5],"\t",$_->[6],"\n";}


#creating the itemsets
# @itemsets : fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, Occurence, Category
# TODO @itemsets : fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, Occurence, Noise1, Noise2, Noise3, Noise4, Noise5, Noise-avg Category

#@relations : fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, category
#sorted on fam1, fam2 and category
$ctr=0;
$j=0;#marker for first iteration
for $i (0 .. $#relations){
	#cleaning up
	$relations[$i][0]=~ s/\s*//g;
	$relations[$i][3]=~ s/\s*//g;
	#for the first time
	if($j==0){
		$fam1=$relations[$i][0];
		$fam2=$relations[$i][3];
		$count_1=$count_2=$count_3=0;
		if ($relations[$i][6] == 1){ $count_1++;}
		elsif ($relations[$i][6] == 2){ $count_2++;}
		elsif ($relations[$i][6] == 3){ $count_3++;}
		$j++;
	}

	#for the rest
	if ($fam1 eq $relations[$i][0] && $fam2 eq $relations[$i][3]){
		if ($relations[$i][6] == 1){ $count_1++;}
		elsif ($relations[$i][6] == 2){ $count_2++;}
		elsif ($relations[$i][6] == 3){ $count_3++;}
	}
	else {
		#write itemsets for that family pair
		$itemsets[$ctr][0]=$fam1;
		$itemsets[$ctr][1]=&get_count($fam1);
		$itemsets[$ctr][2]=&get_avglen($fam1);
		$itemsets[$ctr][3]=$fam2;
		$itemsets[$ctr][4]=&get_count($fam2);
		$itemsets[$ctr][5]=&get_avglen($fam2);
		$itemsets[$ctr][6]=$count_1;
		$itemsets[$ctr++][7]=1;
		$itemsets[$ctr][0]=$fam1;
		$itemsets[$ctr][1]=&get_count($fam1);
		$itemsets[$ctr][2]=&get_avglen($fam1);
		$itemsets[$ctr][3]=$fam2;
		$itemsets[$ctr][4]=&get_count($fam2);
		$itemsets[$ctr][5]=&get_avglen($fam2);
		$itemsets[$ctr][6]=$count_2;
		$itemsets[$ctr++][7]=2;
		$itemsets[$ctr][0]=$fam1;
		$itemsets[$ctr][1]=&get_count($fam1);
		$itemsets[$ctr][2]=&get_avglen($fam1);
		$itemsets[$ctr][3]=$fam2;
		$itemsets[$ctr][4]=&get_count($fam2);
		$itemsets[$ctr][5]=&get_avglen($fam2);
		$itemsets[$ctr][6]=$count_3;
		$itemsets[$ctr++][7]=3;
		
		#initialize for next family pair
		$fam1=$relations[$i][0];
		$fam2=$relations[$i][3];
		$count_1=$count_2=$count_3=0;
		if ($relations[$i][6] == 1){ $count_1++;}
		elsif ($relations[$i][6] == 2){ $count_2++;}
		elsif ($relations[$i][6] == 3){ $count_3++;}
	}
}
#write itemsets for last family pair
$itemsets[$ctr][0]=$fam1;
$itemsets[$ctr][1]=&get_count($fam1);
$itemsets[$ctr][2]=&get_avglen($fam1);
$itemsets[$ctr][3]=$fam2;
$itemsets[$ctr][4]=&get_count($fam2);
$itemsets[$ctr][5]=&get_avglen($fam2);
$itemsets[$ctr][6]=$count_1;
$itemsets[$ctr++][7]=1;
$itemsets[$ctr][0]=$fam1;
$itemsets[$ctr][1]=&get_count($fam1);
$itemsets[$ctr][2]=&get_avglen($fam1);
$itemsets[$ctr][3]=$fam2;
$itemsets[$ctr][4]=&get_count($fam2);
$itemsets[$ctr][5]=&get_avglen($fam2);
$itemsets[$ctr][6]=$count_2;
$itemsets[$ctr++][7]=2;
$itemsets[$ctr][0]=$fam1;
$itemsets[$ctr][1]=&get_count($fam1,@counts);
$itemsets[$ctr][2]=&get_avglen($fam1,@counts);
$itemsets[$ctr][3]=$fam2;
$itemsets[$ctr][4]=&get_count($fam2,@counts);
$itemsets[$ctr][5]=&get_avglen($fam2,@counts);
$itemsets[$ctr][6]=$count_3;
$itemsets[$ctr++][7]=3;

#sorting @itemsets on occurences in decreasing order
@temp = sort {$b->[6] <=> $a->[6]} @itemsets;
@itemsets=@temp;

#printing the itemsets in tabbed format
# @itemsets : fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, Occurence, Category
# TODO @itemsets : fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, Occurence, Noise1, Noise2, Noise3, Noise4, Noise5, Noise-avg Category
print OUTFILEDATA "\#Printing \@itemsets...............\n";
print OUTFILEDATA "\#fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, Occurence, Category\n";
for $i (0 .. $#itemsets){
	if($itemsets[$i][6] > 0){#print only if >0 occurences
		print OUTFILEDATA $itemsets[$i][0],"\t",$itemsets[$i][1],"\t",
		$itemsets[$i][2],"\t",$itemsets[$i][3],"\t",$itemsets[$i][4],"\t",
		$itemsets[$i][5],"\t",$itemsets[$i][6],"\t",$itemsets[$i][7],"\n";
	}
}


close (INFILEDATA);
close (OUTFILEDATA);

exit;
