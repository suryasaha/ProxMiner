#!/usr/bin/perl -w
# MGEL
# Surya Saha 3/15/07
# reading cmd line input .out file which is sorted on the start position
# and finds the relationship among families
# Relationship types: 
# Upstream: u1 (0-500 bases),u2 (500-1000 bases),u3 (1000-5000 bases)
# Downstream: d1 (0-500 bases),d2 (500-1000 bases),d3 (1000-5000 bases)
# In: Location of fam2 is entirely within fam1 (IN)
# Overlap: 
# single linkage algo so consider overlap if > 10% of either
# 10% to 30% (Ovlap-10to30)
# 30% to 70% (Ovlap-30to70)
# 70% + (Ovlap>70)
# Creating the frequent itemsets in the format
# fam1 fam2 Occurence Category

# v3: Removed all duplicate counting
# v3: Counts all relationships 
# v4: Optimized the code to avoid recording itemsets with 0 count
# v4: Check for function call with large parameters

use strict;
use warnings;
use POSIX;

unless (@ARGV == 1){
	print "USAGE: $0 <input .out file> \n";
	exit;
}


my ($ifname,$rec,@temp,@table,@famnames,@relations,@itemsets,@pruned_itemsets,@counts,%temphash,
$ctr,$i,$j,$fam1,$fam2,$count_u1,$count_u2,$count_u3,$count_d1,$count_d2,$count_d3,
$count_IN,$count_OL1030,$count_OL3070,$count_OL70plus,$tot_rels,$tot_isets,$tot_recs,
$tot_fams,$user_t,$system_t,$cuser_t,$csystem_t);

$ifname=$ARGV[0];
chomp $ifname;
unless(open(INFILEDATA,$ifname)){print "not able to open ".$ifname."\n\n";exit;}

unless(open(OUTFILEDATA,">$ifname.f_itemsets.tab")){print "not able to open ".$ifname.".tab \n\n";exit;}

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
$ctr=0;
while($rec=<INFILEDATA>){
	if($rec =~ /#/){next;}
	if(length ($rec) < 10){next;}#for avoiding last line
	push @table, [split(' ',$rec)];
	$ctr++;
}
# record tot recs
$tot_recs = $ctr;

#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\# Runtime details after reading in the file: \n";
print OUTFILEDATA "\# System time for process: $system_t\n";
#print OUTFILEDATA "\# System time for children: $csystem_t\n";
print OUTFILEDATA "\# User time for process: $user_t\n";
#print OUTFILEDATA "\# User time for children: $cuser_t\n";


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
$tot_fams=$ctr;

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

#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\# Runtime details after preparing the \@counts and \@famnames array: \n";
print OUTFILEDATA "\# System time for process: $system_t\n";
#print OUTFILEDATA "\# System time for children: $csystem_t\n";
print OUTFILEDATA "\# User time for process: $user_t\n";
#print OUTFILEDATA "\# User time for children: $cuser_t\n";

#testing
# print "\nPrinting \@counts...............\n";
# foreach(@counts){ print $_->[0],' ',$_->[1],"\t",$_->[2],"\t";}
# print "\n\n For R=997:",&get_count("R=997")," and ",&get_avglen("R=997"),"\n";

#@table
#1935 10.6  0.0  2.8 chr12 8936  9225 27748096 C  R=286 Unknown (0) 283   2
#0    1     2    3   4     5     6    7        8  9     10      11  12    13

#finding the relationships
#@relations : fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, category
$ctr=0;
for $i (0 .. $#table){
	my $ref_start=$table[$i][5];
	my $ref_end=$table[$i][6];
	my $ref_fam=$table[$i][9];
	
	# look for all relationships for $table[$i][] 
	foreach $j (0 .. $#table){
		my $sub_start=$table[$j][5];
		my $sub_end=$table[$j][6];
		my $sub_fam=$table[$j][9];
		
		#skip if $table[$i][] record found
		if ( $ref_fam eq $sub_fam && $ref_start == $sub_start && $ref_end == $sub_end){ next;}
		
		# Note: since all relationship are exclusive, I have used elsif
		# Upstream: u1 (0-500 bases),u2 (500-1000 bases),u3 (1000-5000 bases)
		if (($sub_end <= $ref_start) && ($sub_end > $ref_start - 500)){
			$relations[$ctr][0]=$ref_fam;
			$relations[$ctr][1]=&get_count($ref_fam,@counts);
			$relations[$ctr][2]=&get_avglen($ref_fam,@counts);
			$relations[$ctr][3]=$sub_fam;
			$relations[$ctr][4]=&get_count($sub_fam,@counts);
			$relations[$ctr][5]=&get_avglen($sub_fam,@counts);
			$relations[$ctr++][6]="u1";
		}
		elsif (($sub_end <= $ref_start - 500) && ($sub_end > $ref_start - 1000)){
			$relations[$ctr][0]=$ref_fam;
			$relations[$ctr][1]=&get_count($ref_fam,@counts);
			$relations[$ctr][2]=&get_avglen($ref_fam,@counts);
			$relations[$ctr][3]=$sub_fam;
			$relations[$ctr][4]=&get_count($sub_fam,@counts);
			$relations[$ctr][5]=&get_avglen($sub_fam,@counts);
			$relations[$ctr++][6]="u2";
		}
		elsif (($sub_end <= $ref_start - 1000) && ($sub_end > $ref_start - 5000)){
			$relations[$ctr][0]=$ref_fam;
			$relations[$ctr][1]=&get_count($ref_fam,@counts);
			$relations[$ctr][2]=&get_avglen($ref_fam,@counts);
			$relations[$ctr][3]=$sub_fam;
			$relations[$ctr][4]=&get_count($sub_fam,@counts);
			$relations[$ctr][5]=&get_avglen($sub_fam,@counts);
			$relations[$ctr++][6]="u3";
		}
		# Downstream: d1 (0-500 bases),d2 (500-1000 bases),d3 (1000-5000 bases)
		elsif (($sub_start >= $ref_end) && ($sub_start < $ref_end + 500)){
			$relations[$ctr][0]=$ref_fam;
			$relations[$ctr][1]=&get_count($ref_fam,@counts);
			$relations[$ctr][2]=&get_avglen($ref_fam,@counts);
			$relations[$ctr][3]=$sub_fam;
			$relations[$ctr][4]=&get_count($sub_fam,@counts);
			$relations[$ctr][5]=&get_avglen($sub_fam,@counts);
			$relations[$ctr++][6]="d1";
		}
		elsif (($sub_start >= $ref_end + 500) && ($sub_start < $ref_end + 1000)){
			$relations[$ctr][0]=$ref_fam;
			$relations[$ctr][1]=&get_count($ref_fam,@counts);
			$relations[$ctr][2]=&get_avglen($ref_fam,@counts);
			$relations[$ctr][3]=$sub_fam;
			$relations[$ctr][4]=&get_count($sub_fam,@counts);
			$relations[$ctr][5]=&get_avglen($sub_fam,@counts);
			$relations[$ctr++][6]="d2";
		}
		elsif (($sub_start >= $ref_end + 1000) && ($sub_start < $ref_start + 5000)){
			$relations[$ctr][0]=$ref_fam;
			$relations[$ctr][1]=&get_count($ref_fam,@counts);
			$relations[$ctr][2]=&get_avglen($ref_fam,@counts);
			$relations[$ctr][3]=$sub_fam;
			$relations[$ctr][4]=&get_count($sub_fam,@counts);
			$relations[$ctr][5]=&get_avglen($sub_fam,@counts);
			$relations[$ctr++][6]="d3";
		}
		# In: Location of fam2 is entirely within fam1 (IN)
		elsif (($sub_start >= $ref_start) && ($ref_end >= $sub_end)){
			$relations[$ctr][0]=$ref_fam;
			$relations[$ctr][1]=&get_count($ref_fam,@counts);
			$relations[$ctr][2]=&get_avglen($ref_fam,@counts);
			$relations[$ctr][3]=$sub_fam;
			$relations[$ctr][4]=&get_count($sub_fam,@counts);
			$relations[$ctr][5]=&get_avglen($sub_fam,@counts);
			$relations[$ctr++][6]="IN";
		}
		# Overlap: If overlap is more than 10% of length of either family (Ovlap)
		elsif ((($sub_start < $ref_end) &&  ($sub_start > $ref_start)) 
		|| (($ref_end > $sub_end) && ( $ref_start < $sub_end))){
			my ($ovlap, $ref_ovlap, $sub_ovlap);
			
			#if subject fam starts within the reference fam
			if (($sub_start < $ref_end) &&  ($sub_start > $ref_start)){
				$ovlap = $ref_end - $sub_start;
			}
			#if subject fam ends within the reference fam
			elsif (($ref_end > $sub_end) && ( $ref_start < $sub_end)){
				$ovlap = $sub_end - $ref_start;
			}
			
			$ref_ovlap = ($ovlap / ($ref_end - $ref_start)) * 100;
			$sub_ovlap = ($ovlap / ($sub_end - $sub_start)) * 100;			
			
			# single linkage algo so consider overlap if > 10% of either
			# 10% to 30% (Ovlap-10to30)
			if ((($ref_ovlap > 10.00) && ($ref_ovlap <= 30.00)) || 
			(($sub_ovlap > 10.00) && ($sub_ovlap <= 30.00))) {
				$relations[$ctr][0]=$ref_fam;
				$relations[$ctr][1]=&get_count($ref_fam,@counts);
				$relations[$ctr][2]=&get_avglen($ref_fam,@counts);
				$relations[$ctr][3]=$sub_fam;
				$relations[$ctr][4]=&get_count($sub_fam,@counts);
				$relations[$ctr][5]=&get_avglen($sub_fam,@counts);
				$relations[$ctr++][6]="Ovlap-10to30";
			}
			# 30% to 70% (Ovlap-30to70)
			elsif ((($ref_ovlap > 30.00) && ($ref_ovlap <= 70.00)) || 
			(($sub_ovlap > 30.00) && ($sub_ovlap <= 70.00))) {
				$relations[$ctr][0]=$ref_fam;
				$relations[$ctr][1]=&get_count($ref_fam,@counts);
				$relations[$ctr][2]=&get_avglen($ref_fam,@counts);
				$relations[$ctr][3]=$sub_fam;
				$relations[$ctr][4]=&get_count($sub_fam,@counts);
				$relations[$ctr][5]=&get_avglen($sub_fam,@counts);
				$relations[$ctr++][6]="Ovlap-30to70";
			}
			# 70% + (Ovlap>70)
			elsif (($ref_ovlap > 70.00) || ($sub_ovlap > 70.00)) {
				$relations[$ctr][0]=$ref_fam;
				$relations[$ctr][1]=&get_count($ref_fam,@counts);
				$relations[$ctr][2]=&get_avglen($ref_fam,@counts);
				$relations[$ctr][3]=$sub_fam;
				$relations[$ctr][4]=&get_count($sub_fam,@counts);
				$relations[$ctr][5]=&get_avglen($sub_fam,@counts);
				$relations[$ctr++][6]="Ovlap-70plus";
			}
		}
	}
}

# record tot relations
$tot_rels= $ctr;


#sorting @relations on fam1, fam2 and category
#@relations : fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, category
@temp = sort {($a->[0] cmp $b->[0]) or ($a->[3] cmp $b->[3]) or
($a->[6] cmp $b->[6]) } @relations;
@relations=@temp;

#testing
# print "\nPrinting \@relations...............\n";
# foreach(@relations){ print $_->[0],' ',$_->[1],"\t",$_->[2],"\t",$_->[3],"\t",$_->[4],"\t",
# $_->[5],"\t",$_->[6],"\n";}

#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\# Runtime details after finding the relationships: \n";
print OUTFILEDATA "\# System time for process: $system_t\n";
#print OUTFILEDATA "\# System time for children: $csystem_t\n";
print OUTFILEDATA "\# User time for process: $user_t\n";
#print OUTFILEDATA "\# User time for children: $cuser_t\n";



#creating the itemsets
# @itemsets : fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, Occurence, Category

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
		#Note: have used fam1 and fam2 instead of ref_fam and sub_fam
		$fam1=$relations[$i][0];
		$fam2=$relations[$i][3];
		$count_u1=$count_u2=$count_u3=0;
		$count_d1=$count_d2=$count_d3=0;
		$count_IN=$count_OL1030=$count_OL3070=$count_OL70plus=0;
		
		#proces the record
		if ($relations[$i][6] eq "d1"){ $count_d1++;}
		elsif ($relations[$i][6] eq "d2"){ $count_d2++;}
		elsif ($relations[$i][6] eq "d3"){ $count_d3++;}
		elsif ($relations[$i][6] eq "u1"){ $count_u1++;}
		elsif ($relations[$i][6] eq "u2"){ $count_u2++;}
		elsif ($relations[$i][6] eq "u3"){ $count_u3++;}
		elsif ($relations[$i][6] eq "IN"){ $count_IN++;}
		elsif ($relations[$i][6] eq "Ovlap-10to30"){ $count_OL1030++;}
		elsif ($relations[$i][6] eq "Ovlap-30to70"){ $count_OL3070++;}
		elsif ($relations[$i][6] eq "Ovlap-70plus"){ $count_OL70plus++;}
		$j++;
	}

	#for the rest
	if ($fam1 eq $relations[$i][0] && $fam2 eq $relations[$i][3]){
		if ($relations[$i][6] eq "d1"){ $count_d1++;}
		elsif ($relations[$i][6] eq "d2"){ $count_d2++;}
		elsif ($relations[$i][6] eq "d3"){ $count_d3++;}
		elsif ($relations[$i][6] eq "u1"){ $count_u1++;}
		elsif ($relations[$i][6] eq "u2"){ $count_u2++;}
		elsif ($relations[$i][6] eq "u3"){ $count_u3++;}
		elsif ($relations[$i][6] eq "IN"){ $count_IN++;}
		elsif ($relations[$i][6] eq "Ovlap-10to30"){ $count_OL1030++;}
		elsif ($relations[$i][6] eq "Ovlap-30to70"){ $count_OL3070++;}
		elsif ($relations[$i][6] eq "Ovlap-70plus"){ $count_OL70plus++;}
	}
	else {
		#write itemsets for that family pair
		if ($count_u1 > 0){
			$itemsets[$ctr][0]=$fam1;
			$itemsets[$ctr][1]=&get_count($fam1);
			$itemsets[$ctr][2]=&get_avglen($fam1);
			$itemsets[$ctr][3]=$fam2;
			$itemsets[$ctr][4]=&get_count($fam2);
			$itemsets[$ctr][5]=&get_avglen($fam2);
			$itemsets[$ctr][6]=$count_u1;
			$itemsets[$ctr++][7]="u1";
		}
		
		if ($count_u2 > 0){
			$itemsets[$ctr][0]=$fam1;
			$itemsets[$ctr][1]=&get_count($fam1);
			$itemsets[$ctr][2]=&get_avglen($fam1);
			$itemsets[$ctr][3]=$fam2;
			$itemsets[$ctr][4]=&get_count($fam2);
			$itemsets[$ctr][5]=&get_avglen($fam2);
			$itemsets[$ctr][6]=$count_u2;
			$itemsets[$ctr++][7]="u2";
		}
		
		if ($count_u3 > 0){
			$itemsets[$ctr][0]=$fam1;
			$itemsets[$ctr][1]=&get_count($fam1);
			$itemsets[$ctr][2]=&get_avglen($fam1);
			$itemsets[$ctr][3]=$fam2;
			$itemsets[$ctr][4]=&get_count($fam2);
			$itemsets[$ctr][5]=&get_avglen($fam2);
			$itemsets[$ctr][6]=$count_u3;
			$itemsets[$ctr++][7]="u3";
		}
		
		if ($count_d1 > 0){
			$itemsets[$ctr][0]=$fam1;
			$itemsets[$ctr][1]=&get_count($fam1);
			$itemsets[$ctr][2]=&get_avglen($fam1);
			$itemsets[$ctr][3]=$fam2;
			$itemsets[$ctr][4]=&get_count($fam2);
			$itemsets[$ctr][5]=&get_avglen($fam2);
			$itemsets[$ctr][6]=$count_d1;
			$itemsets[$ctr++][7]="d1";
		}
		
		if ($count_d2 > 0){
			$itemsets[$ctr][0]=$fam1;
			$itemsets[$ctr][1]=&get_count($fam1);
			$itemsets[$ctr][2]=&get_avglen($fam1);
			$itemsets[$ctr][3]=$fam2;
			$itemsets[$ctr][4]=&get_count($fam2);
			$itemsets[$ctr][5]=&get_avglen($fam2);
			$itemsets[$ctr][6]=$count_d2;
			$itemsets[$ctr++][7]="d2";
		}
		
		if ($count_d3 > 0){
			$itemsets[$ctr][0]=$fam1;
			$itemsets[$ctr][1]=&get_count($fam1);
			$itemsets[$ctr][2]=&get_avglen($fam1);
			$itemsets[$ctr][3]=$fam2;
			$itemsets[$ctr][4]=&get_count($fam2);
			$itemsets[$ctr][5]=&get_avglen($fam2);
			$itemsets[$ctr][6]=$count_d3;
			$itemsets[$ctr++][7]="d3";
		}
		
		if ($count_IN > 0){
			$itemsets[$ctr][0]=$fam1;
			$itemsets[$ctr][1]=&get_count($fam1);
			$itemsets[$ctr][2]=&get_avglen($fam1);
			$itemsets[$ctr][3]=$fam2;
			$itemsets[$ctr][4]=&get_count($fam2);
			$itemsets[$ctr][5]=&get_avglen($fam2);
			$itemsets[$ctr][6]=$count_IN;
			$itemsets[$ctr++][7]="IN";
		}
		
		if ($count_OL1030 > 0){
			$itemsets[$ctr][0]=$fam1;
			$itemsets[$ctr][1]=&get_count($fam1);
			$itemsets[$ctr][2]=&get_avglen($fam1);
			$itemsets[$ctr][3]=$fam2;
			$itemsets[$ctr][4]=&get_count($fam2);
			$itemsets[$ctr][5]=&get_avglen($fam2);
			$itemsets[$ctr][6]=$count_OL1030;
			$itemsets[$ctr++][7]="Ovlap-10to30";
		}
		
		if ($count_OL3070 > 0){
			$itemsets[$ctr][0]=$fam1;
			$itemsets[$ctr][1]=&get_count($fam1);
			$itemsets[$ctr][2]=&get_avglen($fam1);
			$itemsets[$ctr][3]=$fam2;
			$itemsets[$ctr][4]=&get_count($fam2);
			$itemsets[$ctr][5]=&get_avglen($fam2);
			$itemsets[$ctr][6]=$count_OL3070;
			$itemsets[$ctr++][7]="Ovlap-30to70";
		}
		
		if ($count_OL70plus > 0){
			$itemsets[$ctr][0]=$fam1;
			$itemsets[$ctr][1]=&get_count($fam1);
			$itemsets[$ctr][2]=&get_avglen($fam1);
			$itemsets[$ctr][3]=$fam2;
			$itemsets[$ctr][4]=&get_count($fam2);
			$itemsets[$ctr][5]=&get_avglen($fam2);
			$itemsets[$ctr][6]=$count_OL70plus;
			$itemsets[$ctr++][7]="Ovlap-70plus";
		}

		#initialize for next family pair
		$fam1=$relations[$i][0];
		$fam2=$relations[$i][3];
		$count_u1=$count_u2=$count_u3=0;
		$count_d1=$count_d2=$count_d3=0;
		$count_IN=$count_OL1030=$count_OL3070=$count_OL70plus=0;
		
		#process the record
		if ($relations[$i][6] eq "d1"){ $count_d1++;}
		elsif ($relations[$i][6] eq "d2"){ $count_d2++;}
		elsif ($relations[$i][6] eq "d3"){ $count_d3++;}
		elsif ($relations[$i][6] eq "u1"){ $count_u1++;}
		elsif ($relations[$i][6] eq "u2"){ $count_u2++;}
		elsif ($relations[$i][6] eq "u3"){ $count_u3++;}
		elsif ($relations[$i][6] eq "IN"){ $count_IN++;}
		elsif ($relations[$i][6] eq "Ovlap-10to30"){ $count_OL1030++;}
		elsif ($relations[$i][6] eq "Ovlap-30to70"){ $count_OL3070++;}
		elsif ($relations[$i][6] eq "Ovlap-70plus"){ $count_OL70plus++;}
	}
}
#write itemsets for last family pair
if ($count_u1 > 0){
	$itemsets[$ctr][0]=$fam1;
	$itemsets[$ctr][1]=&get_count($fam1);
	$itemsets[$ctr][2]=&get_avglen($fam1);
	$itemsets[$ctr][3]=$fam2;
	$itemsets[$ctr][4]=&get_count($fam2);
	$itemsets[$ctr][5]=&get_avglen($fam2);
	$itemsets[$ctr][6]=$count_u1;
	$itemsets[$ctr++][7]="u1";
}

if ($count_u2 > 0){
	$itemsets[$ctr][0]=$fam1;
	$itemsets[$ctr][1]=&get_count($fam1);
	$itemsets[$ctr][2]=&get_avglen($fam1);
	$itemsets[$ctr][3]=$fam2;
	$itemsets[$ctr][4]=&get_count($fam2);
	$itemsets[$ctr][5]=&get_avglen($fam2);
	$itemsets[$ctr][6]=$count_u2;
	$itemsets[$ctr++][7]="u2";
}

if ($count_u3 > 0){
	$itemsets[$ctr][0]=$fam1;
	$itemsets[$ctr][1]=&get_count($fam1);
	$itemsets[$ctr][2]=&get_avglen($fam1);
	$itemsets[$ctr][3]=$fam2;
	$itemsets[$ctr][4]=&get_count($fam2);
	$itemsets[$ctr][5]=&get_avglen($fam2);
	$itemsets[$ctr][6]=$count_u3;
	$itemsets[$ctr++][7]="u3";
}

if ($count_d1 > 0){
	$itemsets[$ctr][0]=$fam1;
	$itemsets[$ctr][1]=&get_count($fam1);
	$itemsets[$ctr][2]=&get_avglen($fam1);
	$itemsets[$ctr][3]=$fam2;
	$itemsets[$ctr][4]=&get_count($fam2);
	$itemsets[$ctr][5]=&get_avglen($fam2);
	$itemsets[$ctr][6]=$count_d1;
	$itemsets[$ctr++][7]="d1";
}

if ($count_d2 > 0){
	$itemsets[$ctr][0]=$fam1;
	$itemsets[$ctr][1]=&get_count($fam1);
	$itemsets[$ctr][2]=&get_avglen($fam1);
	$itemsets[$ctr][3]=$fam2;
	$itemsets[$ctr][4]=&get_count($fam2);
	$itemsets[$ctr][5]=&get_avglen($fam2);
	$itemsets[$ctr][6]=$count_d2;
	$itemsets[$ctr++][7]="d2";
}

if ($count_d3 > 0){
	$itemsets[$ctr][0]=$fam1;
	$itemsets[$ctr][1]=&get_count($fam1);
	$itemsets[$ctr][2]=&get_avglen($fam1);
	$itemsets[$ctr][3]=$fam2;
	$itemsets[$ctr][4]=&get_count($fam2);
	$itemsets[$ctr][5]=&get_avglen($fam2);
	$itemsets[$ctr][6]=$count_d3;
	$itemsets[$ctr++][7]="d3";
}

if ($count_IN > 0){
	$itemsets[$ctr][0]=$fam1;
	$itemsets[$ctr][1]=&get_count($fam1);
	$itemsets[$ctr][2]=&get_avglen($fam1);
	$itemsets[$ctr][3]=$fam2;
	$itemsets[$ctr][4]=&get_count($fam2);
	$itemsets[$ctr][5]=&get_avglen($fam2);
	$itemsets[$ctr][6]=$count_IN;
	$itemsets[$ctr++][7]="IN";
}

if ($count_OL1030 > 0){
	$itemsets[$ctr][0]=$fam1;
	$itemsets[$ctr][1]=&get_count($fam1);
	$itemsets[$ctr][2]=&get_avglen($fam1);
	$itemsets[$ctr][3]=$fam2;
	$itemsets[$ctr][4]=&get_count($fam2);
	$itemsets[$ctr][5]=&get_avglen($fam2);
	$itemsets[$ctr][6]=$count_OL1030;
	$itemsets[$ctr++][7]="Ovlap-10to30";
}

if ($count_OL3070 > 0){
	$itemsets[$ctr][0]=$fam1;
	$itemsets[$ctr][1]=&get_count($fam1);
	$itemsets[$ctr][2]=&get_avglen($fam1);
	$itemsets[$ctr][3]=$fam2;
	$itemsets[$ctr][4]=&get_count($fam2);
	$itemsets[$ctr][5]=&get_avglen($fam2);
	$itemsets[$ctr][6]=$count_OL3070;
	$itemsets[$ctr++][7]="Ovlap-30to70";
}

if ($count_OL70plus > 0){
	$itemsets[$ctr][0]=$fam1;
	$itemsets[$ctr][1]=&get_count($fam1);
	$itemsets[$ctr][2]=&get_avglen($fam1);
	$itemsets[$ctr][3]=$fam2;
	$itemsets[$ctr][4]=&get_count($fam2);
	$itemsets[$ctr][5]=&get_avglen($fam2);
	$itemsets[$ctr][6]=$count_OL70plus;
	$itemsets[$ctr++][7]="Ovlap-70plus";
}

# record tot itemsets
$tot_isets= $ctr;

# @itemsets : fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, Occurence, Category
#sorting @itemsets on occurences in decreasing order
@temp = sort {$b->[6] <=> $a->[6]} @itemsets;
@itemsets=@temp;

#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\# Runtime details after finding the itemsets: \n";
print OUTFILEDATA "\# System time for process: $system_t\n";
#print OUTFILEDATA "\# System time for children: $csystem_t\n";
print OUTFILEDATA "\# User time for process: $user_t\n";
#print OUTFILEDATA "\# User time for children: $cuser_t\n";


#pruning the itemset of reciprocal relations
$ctr=0;

# initializing the array
$pruned_itemsets[$ctr++]=$itemsets[0];

foreach $i (1 .. $#itemsets){#skipping the 1st record
	my ($flag,$ref_rel, $ref_fam1,$ref_fam2,$ref_occ,$sub_fam1,
	$sub_fam2,$sub_occ,$sub_rel);
	
	$ref_fam1=$itemsets[$i][0];
	$ref_fam2=$itemsets[$i][3];
	$ref_occ=$itemsets[$i][6];
	$ref_rel=$itemsets[$i][7];

	#debugging
	#print "Ref fam1:$ref_fam1 fam2:$ref_fam2 occ:$ref_occ rel:$rel\n";
	
	# filtering out records with 0 occurences
	# SHOULD BE NO NEED FOR THIS
	if($ref_occ == 0) {
		print STDERR "\n Surprise from Ref fam1:$ref_fam1 fam2:$ref_fam2 occ:$ref_occ rel:$ref_rel\n";
		next;
	}
	
	# changing to corresponding rel 
	if ($ref_rel eq "u1") { $ref_rel="d1";}
	elsif ($ref_rel eq "u2") { $ref_rel="d2";}
	elsif ($ref_rel eq "u3") { $ref_rel="d3";}
	elsif ($ref_rel eq "d1") { $ref_rel="u1";}
	elsif ($ref_rel eq "d2") { $ref_rel="u2";}
	elsif ($ref_rel eq "d3") { $ref_rel="u3";}
	
	$flag=0;
	# filtering out double counted records
	# looking for reverse relations
	# if duplicate exists in pruned_itemsets, then flag it
	foreach $j (@pruned_itemsets){
		# go to next record if the class of counts has not been 
		# encountered yet
		if ($ref_occ < $j->[6]){
			next;
		}

		$sub_fam1=$j->[0];
		$sub_fam2=$j->[3];
		$sub_occ=$j->[6];
		$sub_rel=$j->[7];
		
		#debugging
		#print "Sub fam1:$sub_fam1 fam2:$sub_fam2 occ:$sub_occ rel:$rel\n";
		
		
		# actually detect the duplicates
		if (($ref_fam1 eq $sub_fam2) && ($ref_fam2 eq $sub_fam1) 
		&& ($ref_occ == $sub_occ) && ($ref_rel eq $sub_rel)){
			$flag=1;
			last;
		}
		# since the occurences are in decreasing order
		# so exit if no reverse found and next set of occurences are
		# encountered
		elsif ($ref_occ > $sub_occ){
			last;
		}
	}
	if($flag == 0){
		$pruned_itemsets[$ctr++]=$itemsets[$i];
	}
}


#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\# Runtime details after pruning the itemsets: \n";
print OUTFILEDATA "\# System time for process: $system_t\n";
#print OUTFILEDATA "\# System time for children: $csystem_t\n";
print OUTFILEDATA "\# User time for process: $user_t\n";
#print OUTFILEDATA "\# User time for children: $cuser_t\n";

#printing the itemsets in tabbed format
# @itemsets : fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, Occurence, Category
print OUTFILEDATA "\# Total records in out file: $tot_recs\n";
print OUTFILEDATA "\# Total families: $tot_fams\n";
print OUTFILEDATA "\# Total relations: $tot_rels\n";
print OUTFILEDATA "\# Total itemsets: $tot_isets\n";
print OUTFILEDATA "\# Total itemsets after pruning: $ctr\n";
# print OUTFILEDATA "\# Runtime details: \n";
# print OUTFILEDATA "\# System time for process: $system_t\n";
# print OUTFILEDATA "\# System time for children: $csystem_t\n";
# print OUTFILEDATA "\# User time for process: $user_t\n";
# print OUTFILEDATA "\# User time for children: $cuser_t\n";
print OUTFILEDATA "\# Printing pruned \@itemsets...............\n";
print OUTFILEDATA "\#fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, Occurence, Category\n";

for $i (0 .. $#pruned_itemsets){
	print OUTFILEDATA $pruned_itemsets[$i][0],"\t",$pruned_itemsets[$i][1],"\t",
	$pruned_itemsets[$i][2],"\t",$pruned_itemsets[$i][3],"\t",$pruned_itemsets[$i][4],"\t",
	$pruned_itemsets[$i][5],"\t",$pruned_itemsets[$i][6],"\t",$pruned_itemsets[$i][7],"\n";
}

#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\# Runtime details after printing the pruned itemsets: \n";
print OUTFILEDATA "\# System time for process: $system_t\n";
#print OUTFILEDATA "\# System time for children: $csystem_t\n";
print OUTFILEDATA "\# User time for process: $user_t\n";
#print OUTFILEDATA "\# User time for children: $cuser_t\n";


close (INFILEDATA);
close (OUTFILEDATA);

exit;
