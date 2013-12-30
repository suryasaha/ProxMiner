#!/usr/bin/perl -w
# MGEL
# Surya Saha 3/15/07
# reading cmd line input .out file which is sorted on the start position
# and finds the relationship among images and families
# Relationship types: 
# Upstream: u1 (0-500 bases),u2 (500-1000 bases),u3 (1000-5000 bases), u4 (5000-10000 bases), u5 (10000-15000 bases)
# Downstream: d1 (0-500 bases),d2 (500-1000 bases),d3 (1000-5000 bases), d4 (5000-10000 bases), d5 (10000-15000 bases)
# In: Location of fam2 is entirely within fam1 (IN)
# Overlap: 
# single linkage algo so consider overlap if > 10% of either
# 10% to 30% (Ovlap-10to30)
# 30% to 70% (Ovlap-30to70)
# 70% + (Ovlap>70)
# Creating the frequent itemsets in the format
# fam1, imageno, fam1-count, fam1-avglen, imagelen, fam2, imageno, fam2-count, fam2-avglen, imagelen, Occurence, Strand, Category



# v3: Removed all duplicate counting
# v3: Counts all relationships 
# v4: Optimized the code to avoid recording itemsets with 0 count
# v4: Check for function call with large parameters

# v5: count relations for images INSTEAD of families
# v5: Use the strand information to calculate the relationships (See rLog)

# v6: Optimize the code (remove duplicates)
# v6: Fixed the bug where false relations were being counted for 'B' strand because of missing ELSE 

# v7: Better progress messages

use strict;
use warnings;
use POSIX;

unless (@ARGV == 1){
	print "USAGE: $0 <input .out file> \n";
	exit;
}

print STDERR "This version hangs with chr12.con.out!!\n";

my  ($ifname,$rec,@temp,%temphash,$ctr,$i,$j,$k,$l,$flag);
my (@table,@famnames,@relations,@itemsets,@pruned_itemsets,@counts,@ref_relations,$ref_fam_index,$ref_fam_img,$ref_fam,$sub_fam_index,$sub_fam,$sub_fam_img,$count_u1,$count_u2,$count_u3,$count_u4,$count_u5,$count_u1_comp,$count_u2_comp,$count_u3_comp,$count_u4_comp,$count_u5_comp,$count_u1_both,$count_u2_both,$count_u3_both,$count_u4_both,$count_u5_both,$count_d1,$count_d2,$count_d3,$count_d4,$count_d5,$count_d1_comp,$count_d2_comp,$count_d3_comp,$count_d4_comp,$count_d5_comp,$count_d1_both,$count_d2_both,$count_d3_both,$count_d4_both,$count_d5_both,$count_IN,$count_IN_comp,$count_IN_both,$count_OL1030,$count_OL3070,$count_OL70plus,$count_OL1030_comp,$count_OL3070_comp,$count_OL70plus_comp,$count_OL1030_both,$count_OL3070_both,$count_OL70plus_both,$tot_rels,$tot_isets,$tot_recs,$tot_fams,$user_t,$system_t,$cuser_t,$csystem_t);

$ifname=$ARGV[0];
chomp $ifname;
unless(open(INFILEDATA,$ifname)){print "not able to open ".$ifname."\n\n";exit;}

unless(open(OUTFILEDATA,">$ifname.f_itemsets.tab")){print "not able to open ".$ifname.".tab \n\n";exit;}

#to get the count of a family using global @counts array
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

#to get the avg length of a family using global @counts array
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

#to get the index position of a family in the @counts array
#it might be faster to just get info from the @counts array 
#once we have the index pos
#params: $fam
sub get_index{
	my ($fam,$ctr);
	$fam=$_[0];
	$fam=~ s/\s*//g;
	$ctr=0;
	foreach (@counts){
		$_->[0] =~ s/\s*//g;
		if ($_->[0] eq $fam){ 
			return $ctr;
			last;
		}
		else{
			$ctr++;
		}
	}
}

# SLURPING IN THE WHOLE .OUT REPORT FILE
$ctr=0;
while($rec=<INFILEDATA>){
	if($rec =~ /#/){next;}
	if(length ($rec) < 10){next;}#for avoiding last line
	push @table, [split(' ',$rec)];
	$ctr++;
}
# record tot recs
$tot_recs = $ctr;

print OUTFILEDATA "\# Version: 7.0\n";

#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\# Runtime details after reading in the file: \n";
print OUTFILEDATA "\# System time for process: $system_t\n";
#print OUTFILEDATA "\# System time for children: $csystem_t\n";
print OUTFILEDATA "\# User time for process: $user_t\n";
#print OUTFILEDATA "\# User time for children: $cuser_t\n";

print STDERR "\# Runtime details after reading $tot_recs from file: \n";
print STDERR "\# System time for process: $system_t\n";
print STDERR "\# User time for process: $user_t\n";


#@table
#1935 10.6  0.0  2.8 chr12 8936  9225 27748096 C  R=286 Unknown (0) 283   2
#0    1     2    3   4     5     6    7        8  9     10      11  12    13

#sorting @table on start position (should already be sorted by RM)
@temp = sort {$a->[5] <=> $b->[5]} @table;
@table=@temp;


# FIND THE NUMBER OF OCCURENCES OF EACH FAMILY
# get family names
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


# INITIALIZING THE @COUNTS 2D ARRAY
# @count: fam occurences avg-len imagenum
$ctr=0;
foreach(@famnames){
	$counts[$ctr][0]=$_;
	#initializing all counters to 0
	$counts[$ctr][1]=0;#occurences
	$counts[$ctr][2]=0;#avg length
	$counts[$ctr++][3]=0;#number of images (mined till now)
}
$tot_fams=$ctr;

# count the number of times a family is found and its avg length
foreach $i (@counts){
	foreach $j (@table){
		if($i->[0] eq $j->[9]){
			$i->[1]++;#occurences
			$i->[2] = $i->[2] + ($j->[6] - $j->[5]);#total length till now
		}
	}
	$i->[2]=floor($i->[2] / $i->[1]);#avg length
}

#@table
#1935 10.6  0.0  2.8 chr12 8936  9225 27748096 C  R=286 Unknown (0) 283   2
#0    1     2    3   4     5     6    7        8  9     10      11  12    13

# Add a field to end of @table 
# where @table[][14]=image number
foreach (@table){
	# since @counts[][3] is initialized to 0
	$_->[14]=1+$counts[&get_index($_->[9])][3]++; 
}


#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\# Runtime details after preparing \@counts and \@famnames arrays and appending \@table array: \n";
print OUTFILEDATA "\# System time for process: $system_t\n";
#print OUTFILEDATA "\# System time for children: $csystem_t\n";
print OUTFILEDATA "\# User time for process: $user_t\n";
#print OUTFILEDATA "\# User time for children: $cuser_t\n";

# TODO (FOR OPTIMIZATION)
# look for all relationships for image at $table[$i][] 
# first look at images located before it but within range of 5000
# now look at images located within it
# now look at images located after it but within range of 5000
# look for relations where $ref_start - 5000 < $sub_start < $ref_start + 5000 




# FINDING ALL RELATIONS
# @table sorted on start position
# 1935 10.6  0.0  2.8 chr12 8936  9225 27748096 C  R=286 Unknown (0) 283   2
# 0    1     2    3   4     5     6    7        8  9     10      11  12    13
# @count: fam occurences avg-len imagenum
# finding the relationships
# @relations : fam1,fam1-img_no, fam1-count, fam1-avglen, fam1-img_no-len, fam2, fam2-img_no, fam2-count, fam2-avglen, fam2-img_no-len, strand (+,C,B), category

$ctr=0;
for $i (0 .. $#table){
	my $ref_start=$table[$i][5]; my $ref_end=$table[$i][6];
	my $ref_strand=$table[$i][8]; my $ref_fam=$table[$i][9];
	my $ref_img=$table[$i][14];
	my $ref_index=&get_index($ref_fam); #to use later
	my $ref_img_len=$ref_end - $ref_start + 1;
	
	print STDERR ':';
	
	foreach $j (0 .. $#table){
		my $sub_start=$table[$j][5]; my $sub_end=$table[$j][6];
		my $sub_strand=$table[$j][8]; my $sub_fam=$table[$j][9];
		my $sub_img=$table[$j][14];
		my $sub_index=&get_index($sub_fam);
		
		#skip if $table[$i][] record found
		if ( $ref_fam eq $sub_fam && $ref_start == $sub_start && $ref_end == $sub_end){ next;}
		
		# Note: since all relationship are exclusive, I have used elsif
		# Upstream: u1 (0-500 bases),u2 (500-1000 bases),u3 (1000-5000 bases), u4 (5000-10000 bases), u5 (10000-15000 bases)
		if (($sub_end <= $ref_start) && ($sub_end > $ref_start - 500)){
			$relations[$ctr][0]=$ref_fam;# ref fam
			$relations[$ctr][1]=$ref_img;# ref image no
			$relations[$ctr][2]=$counts[$ref_index][1];# ref occurences
			$relations[$ctr][3]=$counts[$ref_index][2];# ref avg len
			$relations[$ctr][4]=$ref_img_len;# ref image length
			$relations[$ctr][5]=$sub_fam;# sub fam
			$relations[$ctr][6]=$sub_img;# sub image no
			$relations[$ctr][7]=$counts[$sub_index][1];# sub occurences
			$relations[$ctr][8]=$counts[$sub_index][2];# sub avg len
			$relations[$ctr][9]=$sub_end - $sub_start + 1;# sub image length
			
			if ($ref_strand eq $sub_strand){# for + and C strand
				$relations[$ctr][10]=$ref_strand;# ref strand
				$relations[$ctr++][11]="u1";
			}
			else{
				# irrespective of strand
				$relations[$ctr][10]='B';# for both strands
				$relations[$ctr++][11]="u1";
			}
		}
		elsif (($sub_end <= $ref_start - 500) && ($sub_end > $ref_start - 1000)){
			$relations[$ctr][0]=$ref_fam;# ref fam
			$relations[$ctr][1]=$ref_img;# ref image no
			$relations[$ctr][2]=$counts[$ref_index][1];# ref occurences
			$relations[$ctr][3]=$counts[$ref_index][2];# ref avg len
			$relations[$ctr][4]=$ref_img_len;# ref image length
			$relations[$ctr][5]=$sub_fam;# sub fam
			$relations[$ctr][6]=$sub_img;# sub image no
			$relations[$ctr][7]=$counts[$sub_index][1];# sub occurences
			$relations[$ctr][8]=$counts[$sub_index][2];# sub avg len
			$relations[$ctr][9]=$sub_end - $sub_start + 1;# sub image length
			
			if ($ref_strand eq $sub_strand){
				$relations[$ctr][10]=$ref_strand;# ref strand
				$relations[$ctr++][11]="u2";
			}
			else{
				# irrespective of strand
				$relations[$ctr][10]='B';# for both strands
				$relations[$ctr++][11]="u2";
			}
		}
		elsif (($sub_end <= $ref_start - 1000) && ($sub_end > $ref_start - 5000)) {
			$relations[$ctr][0]=$ref_fam;# ref fam
			$relations[$ctr][1]=$ref_img;# ref image no
			$relations[$ctr][2]=$counts[$ref_index][1];# ref occurences
			$relations[$ctr][3]=$counts[$ref_index][2];# ref avg len
			$relations[$ctr][4]=$ref_img_len;# ref image length
			$relations[$ctr][5]=$sub_fam;# sub fam
			$relations[$ctr][6]=$sub_img;# sub image no
			$relations[$ctr][7]=$counts[$sub_index][1];# sub occurences
			$relations[$ctr][8]=$counts[$sub_index][2];# sub avg len
			$relations[$ctr][9]=$sub_end - $sub_start + 1;# sub image length
			
			if ($ref_strand eq $sub_strand){
				$relations[$ctr][10]=$ref_strand;# ref strand
				$relations[$ctr++][11]="u3";
			}
			else{
				# irrespective of strand
				$relations[$ctr][10]='B';# for both strands
				$relations[$ctr++][11]="u3";
			}
		}
		elsif (($sub_end <= $ref_start - 5000) && ($sub_end > $ref_start - 10000)){
			$relations[$ctr][0]=$ref_fam;# ref fam
			$relations[$ctr][1]=$ref_img;# ref image no
			$relations[$ctr][2]=$counts[$ref_index][1];# ref occurences
			$relations[$ctr][3]=$counts[$ref_index][2];# ref avg len
			$relations[$ctr][4]=$ref_img_len;# ref image length
			$relations[$ctr][5]=$sub_fam;# sub fam
			$relations[$ctr][6]=$sub_img;# sub image no
			$relations[$ctr][7]=$counts[$sub_index][1];# sub occurences
			$relations[$ctr][8]=$counts[$sub_index][2];# sub avg len
			$relations[$ctr][9]=$sub_end - $sub_start + 1;# sub image length
			
			if ($ref_strand eq $sub_strand){
				$relations[$ctr][10]=$ref_strand;# ref strand
				$relations[$ctr++][11]="u4";
			}
			else{
				# irrespective of strand
				$relations[$ctr][10]='B';# for both strands
				$relations[$ctr++][11]="u4";
			}
		}
		elsif (($sub_end <= $ref_start - 10000) && ($sub_end > $ref_start - 15000)){
			$relations[$ctr][0]=$ref_fam;# ref fam
			$relations[$ctr][1]=$ref_img;# ref image no
			$relations[$ctr][2]=$counts[$ref_index][1];# ref occurences
			$relations[$ctr][3]=$counts[$ref_index][2];# ref avg len
			$relations[$ctr][4]=$ref_img_len;# ref image length
			$relations[$ctr][5]=$sub_fam;# sub fam
			$relations[$ctr][6]=$sub_img;# sub image no
			$relations[$ctr][7]=$counts[$sub_index][1];# sub occurences
			$relations[$ctr][8]=$counts[$sub_index][2];# sub avg len
			$relations[$ctr][9]=$sub_end - $sub_start + 1;# sub image length
			
			if ($ref_strand eq $sub_strand){
				$relations[$ctr][10]=$ref_strand;# ref strand
				$relations[$ctr++][11]="u5";
			}
			else{
				# irrespective of strand
				$relations[$ctr][10]='B';# for both strands
				$relations[$ctr++][11]="u5";
			}
		}
		# Downstream: d1 (0-500 bases),d2 (500-1000 bases),d3 (1000-5000 bases), d4 (5000-10000 bases)
		elsif (($sub_start >= $ref_end) && ($sub_start < $ref_end + 500)){
			$relations[$ctr][0]=$ref_fam;# ref fam
			$relations[$ctr][1]=$ref_img;# ref image no
			$relations[$ctr][2]=$counts[$ref_index][1];# ref occurences
			$relations[$ctr][3]=$counts[$ref_index][2];# ref avg len
			$relations[$ctr][4]=$ref_img_len;# ref image length
			$relations[$ctr][5]=$sub_fam;# sub fam
			$relations[$ctr][6]=$sub_img;# sub image no
			$relations[$ctr][7]=$counts[$sub_index][1];# sub occurences
			$relations[$ctr][8]=$counts[$sub_index][2];# sub avg len
			$relations[$ctr][9]=$sub_end - $sub_start + 1;# sub image length
			
			if ($ref_strand eq $sub_strand){
				$relations[$ctr][10]=$ref_strand;# ref strand
				$relations[$ctr++][11]="d1";
			}
			else{
				# irrespective of strand
				$relations[$ctr][10]='B';# for both strands
				$relations[$ctr++][11]="d1";
			}
		}
		elsif (($sub_start >= $ref_end + 500) && ($sub_start < $ref_end + 1000)){
			$relations[$ctr][0]=$ref_fam;# ref fam
			$relations[$ctr][1]=$ref_img;# ref image no
			$relations[$ctr][2]=$counts[$ref_index][1];# ref occurences
			$relations[$ctr][3]=$counts[$ref_index][2];# ref avg len
			$relations[$ctr][4]=$ref_img_len;# ref image length
			$relations[$ctr][5]=$sub_fam;# sub fam
			$relations[$ctr][6]=$sub_img;# sub image no
			$relations[$ctr][7]=$counts[$sub_index][1];# sub occurences
			$relations[$ctr][8]=$counts[$sub_index][2];# sub avg len
			$relations[$ctr][9]=$sub_end - $sub_start + 1;# sub image length
			
			if ($ref_strand eq $sub_strand){
				$relations[$ctr][10]=$ref_strand;# ref strand
				$relations[$ctr++][11]="d2";
			}
			else{
				# irrespective of strand
				$relations[$ctr][10]='B';# for both strands
				$relations[$ctr++][11]="d2";
			}
		}
		elsif (($sub_start >= $ref_end + 1000) && ($sub_start < $ref_start + 5000)){
			$relations[$ctr][0]=$ref_fam;# ref fam
			$relations[$ctr][1]=$ref_img;# ref image no
			$relations[$ctr][2]=$counts[$ref_index][1];# ref occurences
			$relations[$ctr][3]=$counts[$ref_index][2];# ref avg len
			$relations[$ctr][4]=$ref_img_len;# ref image length
			$relations[$ctr][5]=$sub_fam;# sub fam
			$relations[$ctr][6]=$sub_img;# sub image no
			$relations[$ctr][7]=$counts[$sub_index][1];# sub occurences
			$relations[$ctr][8]=$counts[$sub_index][2];# sub avg len
			$relations[$ctr][9]=$sub_end - $sub_start + 1;# sub image length
			
			if ($ref_strand eq $sub_strand){
				$relations[$ctr][10]=$ref_strand;# ref strand
				$relations[$ctr++][11]="d3";
			}
			else{
				# irrespective of strand
				$relations[$ctr][10]='B';# for both strands
				$relations[$ctr++][11]="d3";
			}
		}
		elsif (($sub_start >= $ref_end + 5000) && ($sub_start < $ref_start + 10000)){
			$relations[$ctr][0]=$ref_fam;# ref fam
			$relations[$ctr][1]=$ref_img;# ref image no
			$relations[$ctr][2]=$counts[$ref_index][1];# ref occurences
			$relations[$ctr][3]=$counts[$ref_index][2];# ref avg len
			$relations[$ctr][4]=$ref_img_len;# ref image length
			$relations[$ctr][5]=$sub_fam;# sub fam
			$relations[$ctr][6]=$sub_img;# sub image no
			$relations[$ctr][7]=$counts[$sub_index][1];# sub occurences
			$relations[$ctr][8]=$counts[$sub_index][2];# sub avg len
			$relations[$ctr][9]=$sub_end - $sub_start + 1;# sub image length
			
			if ($ref_strand eq $sub_strand){
				$relations[$ctr][10]=$ref_strand;# ref strand
				$relations[$ctr++][11]="d4";
			}
			else{
				# irrespective of strand
				$relations[$ctr][10]='B';# for both strands
				$relations[$ctr++][11]="d4";
			}
		}
		elsif (($sub_start >= $ref_end + 10000) && ($sub_start < $ref_start + 15000)){
			$relations[$ctr][0]=$ref_fam;# ref fam
			$relations[$ctr][1]=$ref_img;# ref image no
			$relations[$ctr][2]=$counts[$ref_index][1];# ref occurences
			$relations[$ctr][3]=$counts[$ref_index][2];# ref avg len
			$relations[$ctr][4]=$ref_img_len;# ref image length
			$relations[$ctr][5]=$sub_fam;# sub fam
			$relations[$ctr][6]=$sub_img;# sub image no
			$relations[$ctr][7]=$counts[$sub_index][1];# sub occurences
			$relations[$ctr][8]=$counts[$sub_index][2];# sub avg len
			$relations[$ctr][9]=$sub_end - $sub_start + 1;# sub image length
			
			if ($ref_strand eq $sub_strand){
				$relations[$ctr][10]=$ref_strand;# ref strand
				$relations[$ctr++][11]="d5";
			}
			else{
				# irrespective of strand
				$relations[$ctr][10]='B';# for both strands
				$relations[$ctr++][11]="d5";
			}
		}
		# In: Location of fam2 is entirely within fam1 (IN)
		elsif (($sub_start >= $ref_start) && ($ref_end >= $sub_end)){
			$relations[$ctr][0]=$ref_fam;# ref fam
			$relations[$ctr][1]=$ref_img;# ref image no
			$relations[$ctr][2]=$counts[$ref_index][1];# ref occurences
			$relations[$ctr][3]=$counts[$ref_index][2];# ref avg len
			$relations[$ctr][4]=$ref_img_len;# ref image length
			$relations[$ctr][5]=$sub_fam;# sub fam
			$relations[$ctr][6]=$sub_img;# sub image no
			$relations[$ctr][7]=$counts[$sub_index][1];# sub occurences
			$relations[$ctr][8]=$counts[$sub_index][2];# sub avg len
			$relations[$ctr][9]=$sub_end - $sub_start + 1;# sub image length
			
			if ($ref_strand eq $sub_strand){
				$relations[$ctr][10]=$ref_strand;# ref strand
				$relations[$ctr++][11]="IN";
			}
			else{
				# irrespective of strand
				$relations[$ctr][10]='B';# for both strands
				$relations[$ctr++][11]="IN";
			}
		}
		# Overlap: If overlap is more than 10% of length of either family (Ovlap)
		elsif ((($sub_start < $ref_end) &&  ($sub_start > $ref_start)) 
		|| (($ref_end > $sub_end) && ( $ref_start < $sub_end))) {
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
			
			if  ($ref_strand eq $sub_strand){
				# single linkage algo so consider overlap if > 10% of either
				# 10% to 30% (Ovlap-10to30)
				$relations[$ctr][0]=$ref_fam;# ref fam
				$relations[$ctr][1]=$ref_img;# ref image no
				$relations[$ctr][2]=$counts[$ref_index][1];# ref occurences
				$relations[$ctr][3]=$counts[$ref_index][2];# ref avg len
				$relations[$ctr][4]=$ref_img_len;# ref image length
				$relations[$ctr][5]=$sub_fam;# sub fam
				$relations[$ctr][6]=$sub_img;# sub image no
				$relations[$ctr][7]=$counts[$sub_index][1];# sub occurences
				$relations[$ctr][8]=$counts[$sub_index][2];# sub avg len
				$relations[$ctr][9]=$sub_end - $sub_start + 1;# sub image length
				$relations[$ctr][10]=$ref_strand;# ref strand
				
				if ((($ref_ovlap > 10.00) && ($ref_ovlap <= 30.00)) || 
				(($sub_ovlap > 10.00) && ($sub_ovlap <= 30.00))) {
					$relations[$ctr++][11]="Ovlap-10to30";
				}
				# 30% to 70% (Ovlap-30to70)
				elsif ((($ref_ovlap > 30.00) && ($ref_ovlap <= 70.00)) || 
				(($sub_ovlap > 30.00) && ($sub_ovlap <= 70.00))) {
					$relations[$ctr++][11]="Ovlap-30to70";
				}
				# 70% + (Ovlap>70)
				elsif (($ref_ovlap > 70.00) || ($sub_ovlap > 70.00)) {
					$relations[$ctr++][11]="Ovlap-70plus";
				}
			}
			else{
				# irrespective of strand
				$relations[$ctr][0]=$ref_fam;# ref fam
				$relations[$ctr][1]=$ref_img;# ref image no
				$relations[$ctr][2]=$counts[$ref_index][1];# ref occurences
				$relations[$ctr][3]=$counts[$ref_index][2];# ref avg len
				$relations[$ctr][4]=$ref_img_len;# ref image length
				$relations[$ctr][5]=$sub_fam;# sub fam
				$relations[$ctr][6]=$sub_img;# sub image no
				$relations[$ctr][7]=$counts[$sub_index][1];# sub occurences
				$relations[$ctr][8]=$counts[$sub_index][2];# sub avg len
				$relations[$ctr][9]=$sub_end - $sub_start + 1;# sub image length
				$relations[$ctr][10]='B';# for both strands
				
				if ((($ref_ovlap > 10.00) && ($ref_ovlap <= 30.00)) || 
				(($sub_ovlap > 10.00) && ($sub_ovlap <= 30.00))) {
					$relations[$ctr++][11]="Ovlap-10to30";
				}
				# 30% to 70% (Ovlap-30to70)
				elsif ((($ref_ovlap > 30.00) && ($ref_ovlap <= 70.00)) || 
				(($sub_ovlap > 30.00) && ($sub_ovlap <= 70.00))) {
					$relations[$ctr++][11]="Ovlap-30to70";
				}
				# 70% + (Ovlap>70)
				elsif (($ref_ovlap > 70.00) || ($sub_ovlap > 70.00)) {
					$relations[$ctr++][11]="Ovlap-70plus";
				}
			}
		}
		print STDERR '.';
	}
}

# record tot relations
$tot_rels= $ctr;


#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\# Runtime details after finding the relationships: \n";
print OUTFILEDATA "\# System time for process: $system_t\n";
#print OUTFILEDATA "\# System time for children: $csystem_t\n";
print OUTFILEDATA "\# User time for process: $user_t\n";
#print OUTFILEDATA "\# User time for children: $cuser_t\n";


# sorting @relations on fam1 and fam2 
@temp = sort {($a->[0] cmp $b->[0]) or ($a->[5] cmp $b->[5])} @relations;
@relations=@temp;

#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\# Runtime details after sorting the relationships: \n";
print OUTFILEDATA "\# System time for process: $system_t\n";
#print OUTFILEDATA "\# System time for children: $csystem_t\n";
print OUTFILEDATA "\# User time for process: $user_t\n";
#print OUTFILEDATA "\# User time for children: $cuser_t\n";

print STDERR "\n\# Runtime details after finding and sorting the $tot_rels relationships: \n";
print STDERR "\# System time for process: $system_t\n";
print STDERR "\# User time for process: $user_t\n";


# CREATING THE ITEMSETS
# No image info as we are correcting our counting routine
# @itemsets : fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, Occurence, Strand, Category
# @relations : fam1,fam1-img_no, fam1-count, fam1-avglen, fam1-img_no-len, fam2, fam2-img_no, fam2-count, fam2-avglen, fam2-img_no-len, strand (+,C,B), category
# sorted on fam1 and fam2

#marker for first iteration
$ctr=0;$j=0;

for $i (0 .. $#relations){
	#cleaning up
	$relations[$i][1]=~ s/\s*//g;
	$relations[$i][6]=~ s/\s*//g;
	
	#for the first time
	if($j==0){
		# initializing all variables
		$ref_fam=$relations[$i][0];
		$ref_fam_img=$relations[$i][1];
		#$ref_fam_index=&get_index($ref_fam);
		$sub_fam=$relations[$i][5];
		#$sub_fam_img=$relations[$i][6]; # not being used
		#$sub_fam_index=&get_index($sub_fam);
		$count_u1=$count_u2=$count_u3=$count_u4=$count_u5=0;
		$count_u1_comp=$count_u2_comp=$count_u3_comp=$count_u4_comp=$count_u5_comp=0;
		$count_u1_both=$count_u2_both=$count_u3_both=$count_u4_both=$count_u5_both=0;
		$count_d1=$count_d2=$count_d3=$count_d4=$count_d5=0;
		$count_d1_comp=$count_d2_comp=$count_d3_comp=$count_d4_comp=$count_d5_comp=0;
		$count_d1_both=$count_d2_both=$count_d3_both=$count_d4_both=$count_d5_both=0;
		$count_IN=$count_IN_comp=$count_IN_both=0;
		$count_OL1030=$count_OL3070=$count_OL70plus=0;
		$count_OL1030_comp=$count_OL3070_comp=$count_OL70plus_comp=0;
		$count_OL1030_both=$count_OL3070_both=$count_OL70plus_both=0;
		
		# @ref_relations: ref_fam_img, sub_fam_img, strand, category
		# remember which ref image , sub image, strand and category were involved
		$k=0;
		$ref_relations[$k][0]=$relations[$i][1]; $ref_relations[$k][1]=$relations[$i][6];
		$ref_relations[$k][2]=$relations[$i][10]; $ref_relations[$k++][3]=$relations[$i][11];
		
		#process the record
		if ($relations[$i][11] eq "d1" && $relations[$i][10] eq '+'){ $count_d1++;}
		elsif ($relations[$i][11] eq "d1" && $relations[$i][10] eq 'C'){ $count_d1_comp++;}
		elsif ($relations[$i][11] eq "d1" && $relations[$i][10] eq 'B'){ $count_d1_both++;}
		elsif ($relations[$i][11] eq "d2" && $relations[$i][10] eq '+'){ $count_d2++;}
		elsif ($relations[$i][11] eq "d2" && $relations[$i][10] eq 'C'){ $count_d2_comp++;}
		elsif ($relations[$i][11] eq "d2" && $relations[$i][10] eq 'B'){ $count_d2_both++;}
		elsif ($relations[$i][11] eq "d3" && $relations[$i][10] eq '+'){ $count_d3++;}
		elsif ($relations[$i][11] eq "d3" && $relations[$i][10] eq 'C'){ $count_d3_comp++;}
		elsif ($relations[$i][11] eq "d3" && $relations[$i][10] eq 'B'){ $count_d3_both++;}
		elsif ($relations[$i][11] eq "d4" && $relations[$i][10] eq '+'){ $count_d4++;}
		elsif ($relations[$i][11] eq "d4" && $relations[$i][10] eq 'C'){ $count_d4_comp++;}
		elsif ($relations[$i][11] eq "d4" && $relations[$i][10] eq 'B'){ $count_d4_both++;}
		elsif ($relations[$i][11] eq "d5" && $relations[$i][10] eq '+'){ $count_d5++;}
		elsif ($relations[$i][11] eq "d5" && $relations[$i][10] eq 'C'){ $count_d5_comp++;}
		elsif ($relations[$i][11] eq "d5" && $relations[$i][10] eq 'B'){ $count_d5_both++;}
		elsif ($relations[$i][11] eq "u1" && $relations[$i][10] eq 'C'){ $count_u1_comp++;}
		elsif ($relations[$i][11] eq "u1" && $relations[$i][10] eq 'B'){ $count_u1_both++;}
		elsif ($relations[$i][11] eq "u2" && $relations[$i][10] eq '+'){ $count_u2++;}
		elsif ($relations[$i][11] eq "u2" && $relations[$i][10] eq 'C'){ $count_u2_comp++;}
		elsif ($relations[$i][11] eq "u2" && $relations[$i][10] eq 'B'){ $count_u2_both++;}
		elsif ($relations[$i][11] eq "u3" && $relations[$i][10] eq '+'){ $count_u3++;}
		elsif ($relations[$i][11] eq "u3" && $relations[$i][10] eq 'C'){ $count_u3_comp++;}
		elsif ($relations[$i][11] eq "u3" && $relations[$i][10] eq 'B'){ $count_u3_both++;}
		elsif ($relations[$i][11] eq "u4" && $relations[$i][10] eq '+'){ $count_u4++;}
		elsif ($relations[$i][11] eq "u4" && $relations[$i][10] eq 'C'){ $count_u4_comp++;}
		elsif ($relations[$i][11] eq "u4" && $relations[$i][10] eq 'B'){ $count_u4_both++;}
		elsif ($relations[$i][11] eq "u5" && $relations[$i][10] eq '+'){ $count_u5++;}
		elsif ($relations[$i][11] eq "u5" && $relations[$i][10] eq 'C'){ $count_u5_comp++;}
		elsif ($relations[$i][11] eq "u5" && $relations[$i][10] eq 'B'){ $count_u5_both++;}
		elsif ($relations[$i][11] eq "IN" && $relations[$i][10] eq '+'){ $count_IN++;}
		elsif ($relations[$i][11] eq "IN" && $relations[$i][10] eq 'C'){ $count_IN_comp++;}
		elsif ($relations[$i][11] eq "IN" && $relations[$i][10] eq 'B'){ $count_IN_both++;}
		elsif ($relations[$i][11] eq "Ovlap-10to30" && $relations[$i][10] eq '+'){ $count_OL1030++;}
		elsif ($relations[$i][11] eq "Ovlap-10to30" && $relations[$i][10] eq 'C'){ $count_OL1030_comp++;}
		elsif ($relations[$i][11] eq "Ovlap-10to30" && $relations[$i][10] eq 'B'){ $count_OL1030_both++;}
		elsif ($relations[$i][11] eq "Ovlap-30to70" && $relations[$i][10] eq '+'){ $count_OL3070++;}
		elsif ($relations[$i][11] eq "Ovlap-30to70" && $relations[$i][10] eq 'C'){ $count_OL3070_comp++;}
		elsif ($relations[$i][11] eq "Ovlap-30to70" && $relations[$i][10] eq 'B'){ $count_OL3070_both++;}
		elsif ($relations[$i][11] eq "Ovlap-70plus" && $relations[$i][10] eq '+'){ $count_OL70plus++;}
		elsif ($relations[$i][11] eq "Ovlap-70plus" && $relations[$i][10] eq 'C'){ $count_OL70plus_comp++;}
		elsif ($relations[$i][11] eq "Ovlap-70plus" && $relations[$i][10] eq 'B'){ $count_OL70plus_both++;}
		
		$j++;# increment ctr
	}
	elsif (($ref_fam eq $relations[$i][0]) && ($sub_fam eq $relations[$i][5])){ 
		# for other relations between the same families
		
		# check if this $ref_fam_img already had this relation with $sub_fam for this strand
		# or if this $sub_fam_img already had this relation with $ref_fam for this strand
		$flag=0;
		for $l (0 .. $#ref_relations){
			#debugging
# 			if (!exists $ref_relations[$l][0]){ print "rrel 0 missing!\n";}
			
			# @ref_relations: ref_fam_img, sub_fam_img, strand, category
			if (($ref_relations[$l][0] eq $relations[$i][1]) &&
			($ref_relations[$l][2] eq $relations[$i][10]) &&
			($ref_relations[$l][3] eq $relations[$i][11])){
				# means that the ref image already has this relation 
				# for this sub family on this strand
				$flag=1; last;
			}
			elsif (($ref_relations[$l][1] eq $relations[$i][6]) &&
			($ref_relations[$l][2] eq $relations[$i][10]) &&
			($ref_relations[$l][3] eq $relations[$i][11])){
				# means that the sub image already has this relation 
				# for this ref family on this strand
				$flag=1; last;
			}
		}
		
		# if this relation is new
		if ($flag == 0){
			#add to @ref_relations
			$ref_relations[$k][0]=$relations[$i][1];
			$ref_relations[$k][1]=$relations[$i][6];
			$ref_relations[$k][2]=$relations[$i][10];
			
			if ($relations[$i][11] eq "d1"){
				if ($relations[$i][10] eq '+'){ $count_d1++;}
				elsif ($relations[$i][10] eq 'C'){ $count_d1_comp++;}
				elsif ($relations[$i][10] eq 'B'){ $count_d1_both++;}
				
				#add to @ref_relations
				$ref_relations[$k++][3]=$relations[$i][11];
			}
			elsif ($relations[$i][11] eq "d2"){ 
				if ($relations[$i][10] eq '+'){ $count_d2++;}
				elsif ($relations[$i][10] eq 'C'){ $count_d2_comp++;}
				elsif ($relations[$i][10] eq 'B'){ $count_d2_both++;}

				#add to @ref_relations
				$ref_relations[$k++][3]=$relations[$i][11];
			}
			elsif ($relations[$i][11] eq "d3"){ 
				if ($relations[$i][10] eq '+'){ $count_d3++;}
				elsif ($relations[$i][10] eq 'C'){ $count_d3_comp++;}
				elsif ($relations[$i][10] eq 'B'){ $count_d3_both++;}

				#add to @ref_relations
				$ref_relations[$k++][3]=$relations[$i][11];
			}
			elsif ($relations[$i][11] eq "d4"){ 
				if ($relations[$i][10] eq '+'){ $count_d4++;}
				elsif ($relations[$i][10] eq 'C'){ $count_d4_comp++;}
				elsif ($relations[$i][10] eq 'B'){ $count_d4_both++;}

				#add to @ref_relations
				$ref_relations[$k++][3]=$relations[$i][11];
			}
			elsif ($relations[$i][11] eq "d5"){ 
				if ($relations[$i][10] eq '+'){ $count_d5++;}
				elsif ($relations[$i][10] eq 'C'){ $count_d5_comp++;}
				elsif ($relations[$i][10] eq 'B'){ $count_d5_both++;}

				#add to @ref_relations
				$ref_relations[$k++][3]=$relations[$i][11];
			}
			elsif ($relations[$i][11] eq "u1"){ 
				if ($relations[$i][10] eq '+'){ $count_u1++;}
				elsif ($relations[$i][10] eq 'C'){ $count_u1_comp++;}
				elsif ($relations[$i][10] eq 'B'){ $count_u1_both++;}

				#add to @ref_relations
				$ref_relations[$k++][3]=$relations[$i][11];
			}
			elsif ($relations[$i][11] eq "u2"){ 
				if ($relations[$i][10] eq '+'){ $count_u2++;}
				elsif ($relations[$i][10] eq 'C'){ $count_u2_comp++;}
				elsif ($relations[$i][10] eq 'B'){ $count_u2_both++;}

				#add to @ref_relations
				$ref_relations[$k++][3]=$relations[$i][11];
			}
			elsif ($relations[$i][11] eq "u3"){ 
				if ($relations[$i][10] eq '+'){ $count_u3++;}
				elsif ($relations[$i][10] eq 'C'){ $count_u3_comp++;}
				elsif ($relations[$i][10] eq 'B'){ $count_u3_both++;}

				#add to @ref_relations
				$ref_relations[$k++][3]=$relations[$i][11];
			}
			elsif ($relations[$i][11] eq "u4"){ 
				if ($relations[$i][10] eq '+'){ $count_u4++;}
				elsif ($relations[$i][10] eq 'C'){ $count_u4_comp++;}
				elsif ($relations[$i][10] eq 'B'){ $count_u4_both++;}

				#add to @ref_relations
				$ref_relations[$k++][3]=$relations[$i][11];
			}
			elsif ($relations[$i][11] eq "u5"){ 
				if ($relations[$i][10] eq '+'){ $count_u5++;}
				elsif ($relations[$i][10] eq 'C'){ $count_u5_comp++;}
				elsif ($relations[$i][10] eq 'B'){ $count_u5_both++;}

				#add to @ref_relations
				$ref_relations[$k++][3]=$relations[$i][11];
			}
			elsif ($relations[$i][11] eq "IN"){ 
				if ($relations[$i][10] eq '+'){ $count_IN++;}
				elsif ($relations[$i][10] eq 'C'){ $count_IN_comp++;}
				elsif ($relations[$i][10] eq 'B'){ $count_IN_both++;}

				#add to @ref_relations
				$ref_relations[$k++][3]=$relations[$i][11];
			}
			elsif ($relations[$i][11] eq "Ovlap-10to30"){ 
				if ($relations[$i][10] eq '+'){ $count_OL1030++;}
				elsif ($relations[$i][10] eq 'C'){ $count_OL1030_comp++;}
				elsif ($relations[$i][10] eq 'B'){ $count_OL1030_both++;}

				#add to @ref_relations
				$ref_relations[$k++][3]=$relations[$i][11];
				}
			elsif ($relations[$i][11] eq "Ovlap-30to70"){ 
				if ($relations[$i][10] eq '+'){ $count_OL3070++;}
				elsif ($relations[$i][10] eq 'C'){ $count_OL3070_comp++;}
				elsif ($relations[$i][10] eq 'B'){ $count_OL3070_both++;}

				#add to @ref_relations
				$ref_relations[$k++][3]=$relations[$i][11];
				}
			elsif ($relations[$i][11] eq "Ovlap-70plus"){ 
				if ($relations[$i][10] eq '+'){ $count_OL70plus++;}
				elsif ($relations[$i][10] eq 'C'){ $count_OL70plus_comp++;}
				elsif ($relations[$i][10] eq 'B'){ $count_OL70plus_both++;}

				#add to @ref_relations
				$ref_relations[$k++][3]=$relations[$i][11];
			}
		}
	}
	elsif ($sub_fam ne $relations[$i][5]){ 
		# for the next family pair
		
		#the current @relations record is for the next family pair
		#write itemsets for previous family pair
		if ($count_u1 > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_u1;
			$itemsets[$ctr][7]="+";
			$itemsets[$ctr++][8]="u1";
		}
		if ($count_u1_comp > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_u1_comp;
			$itemsets[$ctr][7]="C";
			$itemsets[$ctr++][8]="u1";
		}
		if ($count_u1_both > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_u1_both;
			$itemsets[$ctr][7]="B";
			$itemsets[$ctr++][8]="u1";
		}
		if ($count_u2 > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_u2;
			$itemsets[$ctr][7]="+";
			$itemsets[$ctr++][8]="u2";
		}
		if ($count_u2_comp > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_u2_comp;
			$itemsets[$ctr][7]="C";
			$itemsets[$ctr++][8]="u2";
		}
		if ($count_u2_both > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_u2_both;
			$itemsets[$ctr][7]="B";
			$itemsets[$ctr++][8]="u2";
		}
		if ($count_u3 > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_u3;
			$itemsets[$ctr][7]="+";
			$itemsets[$ctr++][8]="u3";
		}
		if ($count_u3_comp > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_u3_comp;
			$itemsets[$ctr][7]="C";
			$itemsets[$ctr++][8]="u3";

		}
		if ($count_u3_both > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_u3_both;
			$itemsets[$ctr][7]="B";
			$itemsets[$ctr++][8]="u3";
		}
		if ($count_u4 > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_u4;
			$itemsets[$ctr][7]="+";
			$itemsets[$ctr++][8]="u4";
		}
		if ($count_u4_comp > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_u4_comp;
			$itemsets[$ctr][7]="C";
			$itemsets[$ctr++][8]="u4";
		}
		if ($count_u4_both > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_u4_both;
			$itemsets[$ctr][7]="B";
			$itemsets[$ctr++][8]="u4";
		}
		if ($count_u5 > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_u5;
			$itemsets[$ctr][7]="+";
			$itemsets[$ctr++][8]="u5";
		}
		if ($count_u5_comp > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_u5_comp;
			$itemsets[$ctr][7]="C";
			$itemsets[$ctr++][8]="u5";
		}
		if ($count_u5_both > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_u5_both;
			$itemsets[$ctr][7]="B";
			$itemsets[$ctr++][8]="u5";
		}
		if ($count_d1 > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_d1;
			$itemsets[$ctr][7]="+";
			$itemsets[$ctr++][8]="d1";
		}
		if ($count_d1_comp > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_d1_comp;
			$itemsets[$ctr][7]="C";
			$itemsets[$ctr++][8]="d1";
		}
		if ($count_d1_both > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_d1_both;
			$itemsets[$ctr][7]="B";
			$itemsets[$ctr++][8]="d1";
		}
		if ($count_d2 > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_d2;
			$itemsets[$ctr][7]="+";
			$itemsets[$ctr++][8]="d2";
		}
		if ($count_d2_comp > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_d2_comp;
			$itemsets[$ctr][7]="C";
			$itemsets[$ctr++][8]="d2";
		}
		if ($count_d2_both > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_d2_both;
			$itemsets[$ctr][7]="B";
			$itemsets[$ctr++][8]="d2";
		}
		if ($count_d3 > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_d3;
			$itemsets[$ctr][7]="+";
			$itemsets[$ctr++][8]="d3";
		}
		if ($count_d3_comp > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_d3_comp;
			$itemsets[$ctr][7]="C";
			$itemsets[$ctr++][8]="d3";
		}
		if ($count_d3_both > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_d3_both;
			$itemsets[$ctr][7]="B";
			$itemsets[$ctr++][8]="d3";
		}
		if ($count_d4 > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_d4;
			$itemsets[$ctr][7]="+";
			$itemsets[$ctr++][8]="d4";
		}
		if ($count_d4_comp > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_d4_comp;
			$itemsets[$ctr][7]="C";
			$itemsets[$ctr++][8]="d4";
		}
		if ($count_d4_both > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_d4_both;
			$itemsets[$ctr][7]="B";
			$itemsets[$ctr++][8]="d4";
		}
		if ($count_d5 > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_d5;
			$itemsets[$ctr][7]="+";
			$itemsets[$ctr++][8]="d5";
		}
		if ($count_d5_comp > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_d5_comp;
			$itemsets[$ctr][7]="C";
			$itemsets[$ctr++][8]="d5";
		}
		if ($count_d5_both > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_d5_both;
			$itemsets[$ctr][7]="B";
			$itemsets[$ctr++][8]="d5";
		}
		if ($count_IN > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_IN;
			$itemsets[$ctr][7]="+";
			$itemsets[$ctr++][8]="IN";
		}
		if ($count_IN_comp > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_IN_comp;
			$itemsets[$ctr][7]="C";
			$itemsets[$ctr++][8]="IN";
		}
		if ($count_IN_both > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_IN_both;
			$itemsets[$ctr][7]="B";
			$itemsets[$ctr++][8]="IN";
		}
		if ($count_OL1030 > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_OL1030;
			$itemsets[$ctr][7]="+";
			$itemsets[$ctr++][8]="Ovlap-10to30";
		}
		if ($count_OL1030_comp > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_OL1030_comp;
			$itemsets[$ctr][7]="C";
			$itemsets[$ctr++][8]="Ovlap-10to30";
		}
		if ($count_OL1030_both > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_OL1030_both;
			$itemsets[$ctr][7]="B";
			$itemsets[$ctr++][8]="Ovlap-10to30";
		}
		if ($count_OL3070 > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_OL3070;
			$itemsets[$ctr][7]="+";
			$itemsets[$ctr++][8]="Ovlap-30to70";
		}
		if ($count_OL3070_comp > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_OL3070_comp;
			$itemsets[$ctr][7]="C";
			$itemsets[$ctr++][8]="Ovlap-30to70";
		}
		if ($count_OL3070_both > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_OL3070_both;
			$itemsets[$ctr][7]="B";
			$itemsets[$ctr++][8]="Ovlap-30to70";
		}
		if ($count_OL70plus > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_OL70plus;
			$itemsets[$ctr][7]="+";
			$itemsets[$ctr++][8]="Ovlap-70plus";
		}
		if ($count_OL70plus_comp > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_OL70plus_comp;
			$itemsets[$ctr][7]="C";
			$itemsets[$ctr++][8]="Ovlap-70plus";
		}
		if ($count_OL70plus_both > 0){
			$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
			$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
			$itemsets[$ctr][6]=$count_OL70plus_both;
			$itemsets[$ctr][7]="B";
			$itemsets[$ctr++][8]="Ovlap-70plus";
		}

		#initialize for next family pair
		$ref_fam=$relations[$i][0];
		$sub_fam=$relations[$i][5];
		$count_u1=$count_u2=$count_u3=$count_u4=$count_u5=0;
		$count_u1_comp=$count_u2_comp=$count_u3_comp=$count_u4_comp=$count_u5_comp=0;
		$count_u1_both=$count_u2_both=$count_u3_both=$count_u4_both=$count_u5_both=0;
		$count_d1=$count_d2=$count_d3=$count_d4=$count_d5=0;
		$count_d1_comp=$count_d2_comp=$count_d3_comp=$count_d4_comp=$count_d5_comp=0;
		$count_d1_both=$count_d2_both=$count_d3_both=$count_d4_both=$count_d5_both=0;
		$count_IN=$count_IN_comp=$count_IN_both=0;
		$count_OL1030=$count_OL3070=$count_OL70plus=0;
		$count_OL1030_comp=$count_OL3070_comp=$count_OL70plus_comp=0;
		$count_OL1030_both=$count_OL3070_both=$count_OL70plus_both=0;
		
		# reset @ref_relations
		undef(@ref_relations);
		# @ref_relations: ref_fam_img, sub_fam_img, strand, category
		# remember which ref image , sub image, strand and category were involved
		$k=0;
		$ref_relations[$k][0]=$relations[$i][1]; $ref_relations[$k][1]=$relations[$i][6];
		$ref_relations[$k][2]=$relations[$i][10]; $ref_relations[$k++][3]=$relations[$i][11];
		
		#process the record
		if ($relations[$i][11] eq "d1" && $relations[$i][10] eq '+'){ $count_d1++;}
		elsif ($relations[$i][11] eq "d1" && $relations[$i][10] eq 'C'){ $count_d1_comp++;}
		elsif ($relations[$i][11] eq "d1" && $relations[$i][10] eq 'B'){ $count_d1_both++;}
		elsif ($relations[$i][11] eq "d2" && $relations[$i][10] eq '+'){ $count_d2++;}
		elsif ($relations[$i][11] eq "d2" && $relations[$i][10] eq 'C'){ $count_d2_comp++;}
		elsif ($relations[$i][11] eq "d2" && $relations[$i][10] eq 'B'){ $count_d2_both++;}
		elsif ($relations[$i][11] eq "d3" && $relations[$i][10] eq '+'){ $count_d3++;}
		elsif ($relations[$i][11] eq "d3" && $relations[$i][10] eq 'C'){ $count_d3_comp++;}
		elsif ($relations[$i][11] eq "d3" && $relations[$i][10] eq 'B'){ $count_d3_both++;}
		elsif ($relations[$i][11] eq "d4" && $relations[$i][10] eq '+'){ $count_d4++;}
		elsif ($relations[$i][11] eq "d4" && $relations[$i][10] eq 'C'){ $count_d4_comp++;}
		elsif ($relations[$i][11] eq "d4" && $relations[$i][10] eq 'B'){ $count_d4_both++;}
		elsif ($relations[$i][11] eq "d5" && $relations[$i][10] eq '+'){ $count_d5++;}
		elsif ($relations[$i][11] eq "d5" && $relations[$i][10] eq 'C'){ $count_d5_comp++;}
		elsif ($relations[$i][11] eq "d5" && $relations[$i][10] eq 'B'){ $count_d5_both++;}
		elsif ($relations[$i][11] eq "u1" && $relations[$i][10] eq 'C'){ $count_u1_comp++;}
		elsif ($relations[$i][11] eq "u1" && $relations[$i][10] eq 'B'){ $count_u1_both++;}
		elsif ($relations[$i][11] eq "u2" && $relations[$i][10] eq '+'){ $count_u2++;}
		elsif ($relations[$i][11] eq "u2" && $relations[$i][10] eq 'C'){ $count_u2_comp++;}
		elsif ($relations[$i][11] eq "u2" && $relations[$i][10] eq 'B'){ $count_u2_both++;}
		elsif ($relations[$i][11] eq "u3" && $relations[$i][10] eq '+'){ $count_u3++;}
		elsif ($relations[$i][11] eq "u3" && $relations[$i][10] eq 'C'){ $count_u3_comp++;}
		elsif ($relations[$i][11] eq "u3" && $relations[$i][10] eq 'B'){ $count_u3_both++;}
		elsif ($relations[$i][11] eq "u4" && $relations[$i][10] eq '+'){ $count_u4++;}
		elsif ($relations[$i][11] eq "u4" && $relations[$i][10] eq 'C'){ $count_u4_comp++;}
		elsif ($relations[$i][11] eq "u4" && $relations[$i][10] eq 'B'){ $count_u4_both++;}
		elsif ($relations[$i][11] eq "u5" && $relations[$i][10] eq '+'){ $count_u5++;}
		elsif ($relations[$i][11] eq "u5" && $relations[$i][10] eq 'C'){ $count_u5_comp++;}
		elsif ($relations[$i][11] eq "u5" && $relations[$i][10] eq 'B'){ $count_u5_both++;}
		elsif ($relations[$i][11] eq "IN" && $relations[$i][10] eq '+'){ $count_IN++;}
		elsif ($relations[$i][11] eq "IN" && $relations[$i][10] eq 'C'){ $count_IN_comp++;}
		elsif ($relations[$i][11] eq "IN" && $relations[$i][10] eq 'B'){ $count_IN_both++;}
		elsif ($relations[$i][11] eq "Ovlap-10to30" && $relations[$i][10] eq '+'){ $count_OL1030++;}
		elsif ($relations[$i][11] eq "Ovlap-10to30" && $relations[$i][10] eq 'C'){ $count_OL1030_comp++;}
		elsif ($relations[$i][11] eq "Ovlap-10to30" && $relations[$i][10] eq 'B'){ $count_OL1030_both++;}
		elsif ($relations[$i][11] eq "Ovlap-30to70" && $relations[$i][10] eq '+'){ $count_OL3070++;}
		elsif ($relations[$i][11] eq "Ovlap-30to70" && $relations[$i][10] eq 'C'){ $count_OL3070_comp++;}
		elsif ($relations[$i][11] eq "Ovlap-30to70" && $relations[$i][10] eq 'B'){ $count_OL3070_both++;}
		elsif ($relations[$i][11] eq "Ovlap-70plus" && $relations[$i][10] eq '+'){ $count_OL70plus++;}
		elsif ($relations[$i][11] eq "Ovlap-70plus" && $relations[$i][10] eq 'C'){ $count_OL70plus_comp++;}
		elsif ($relations[$i][11] eq "Ovlap-70plus" && $relations[$i][10] eq 'B'){ $count_OL70plus_both++;}
	}
	print STDERR '.';
}
#write itemsets for last family pair
if ($count_u1 > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_u1; $itemsets[$ctr][7]="+"; $itemsets[$ctr++][8]="u1";
}
if ($count_u1_comp > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_u1_comp; $itemsets[$ctr][7]="C"; $itemsets[$ctr++][8]="u1"; 
}
if ($count_u1_both > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_u1_both; $itemsets[$ctr][7]="B"; $itemsets[$ctr++][8]="u1";
}
if ($count_u2 > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_u2; $itemsets[$ctr][7]="+"; $itemsets[$ctr++][8]="u2";
}
if ($count_u2_comp > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_u2_comp; $itemsets[$ctr][7]="C"; $itemsets[$ctr++][8]="u2";
}
if ($count_u2_both > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_u2_both; $itemsets[$ctr][7]="B"; $itemsets[$ctr++][8]="u2";
}
if ($count_u3 > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_u3; $itemsets[$ctr][7]="+"; $itemsets[$ctr++][8]="u3";
}
if ($count_u3_comp > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_u3_comp; $itemsets[$ctr][7]="C"; $itemsets[$ctr++][8]="u3";

}
if ($count_u3_both > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_u3_both; $itemsets[$ctr][7]="B"; $itemsets[$ctr++][8]="u3";
}
if ($count_u4 > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_u4; $itemsets[$ctr][7]="+"; $itemsets[$ctr++][8]="u4";
}
if ($count_u4_comp > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_u4_comp; $itemsets[$ctr][7]="C"; $itemsets[$ctr++][8]="u4";
}
if ($count_u4_both > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_u4_both; $itemsets[$ctr][7]="B"; $itemsets[$ctr++][8]="u4";
}
if ($count_u5 > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_u5; $itemsets[$ctr][7]="+"; $itemsets[$ctr++][8]="u5";
}
if ($count_u5_comp > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_u5_comp; $itemsets[$ctr][7]="C"; $itemsets[$ctr++][8]="u5";
}
if ($count_u5_both > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_u5_both; $itemsets[$ctr][7]="B"; $itemsets[$ctr++][8]="u5";
}
if ($count_d1 > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_d1; $itemsets[$ctr][7]="+"; $itemsets[$ctr++][8]="d1";
}
if ($count_d1_comp > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_d1_comp; $itemsets[$ctr][7]="C"; $itemsets[$ctr++][8]="d1";
}
if ($count_d1_both > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_d1_both; $itemsets[$ctr][7]="B"; $itemsets[$ctr++][8]="d1";
}
if ($count_d2 > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_d2; $itemsets[$ctr][7]="+"; $itemsets[$ctr++][8]="d2";
}
if ($count_d2_comp > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_d2_comp; $itemsets[$ctr][7]="C"; $itemsets[$ctr++][8]="d2";
}
if ($count_d2_both > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_d2_both; $itemsets[$ctr][7]="B"; $itemsets[$ctr++][8]="d2";
}
if ($count_d3 > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_d3; $itemsets[$ctr][7]="+"; $itemsets[$ctr++][8]="d3";
}
if ($count_d3_comp > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_d3_comp; $itemsets[$ctr][7]="C"; $itemsets[$ctr++][8]="d3";
}
if ($count_d3_both > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_d3_both; $itemsets[$ctr][7]="B"; $itemsets[$ctr++][8]="d3";
}
if ($count_d4 > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_d4; $itemsets[$ctr][7]="+"; $itemsets[$ctr++][8]="d4";
}
if ($count_d4_comp > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_d4_comp; $itemsets[$ctr][7]="C"; $itemsets[$ctr++][8]="d4";
}
if ($count_d4_both > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_d4_both; $itemsets[$ctr][7]="B"; $itemsets[$ctr++][8]="d4";
}
if ($count_d5 > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_d5; $itemsets[$ctr][7]="+"; $itemsets[$ctr++][8]="d5";
}
if ($count_d5_comp > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_d5_comp; $itemsets[$ctr][7]="C"; $itemsets[$ctr++][8]="d5";
}
if ($count_d5_both > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_d5_both; $itemsets[$ctr][7]="B"; $itemsets[$ctr++][8]="d5";
}
if ($count_IN > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_IN; $itemsets[$ctr][7]="+"; $itemsets[$ctr++][8]="IN";
}
if ($count_IN_comp > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_IN_comp; $itemsets[$ctr][7]="C"; $itemsets[$ctr++][8]="IN";
}
if ($count_IN_both > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_IN_both; $itemsets[$ctr][7]="B"; $itemsets[$ctr++][8]="IN";
}
if ($count_OL1030 > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_OL1030; $itemsets[$ctr][7]="+"; $itemsets[$ctr++][8]="Ovlap-10to30";
}
if ($count_OL1030_comp > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_OL1030_comp; $itemsets[$ctr][7]="C"; $itemsets[$ctr++][8]="Ovlap-10to30";
}
if ($count_OL1030_both > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_OL1030_both; $itemsets[$ctr][7]="B"; $itemsets[$ctr++][8]="Ovlap-10to30";
}
if ($count_OL3070 > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_OL3070; $itemsets[$ctr][7]="+"; $itemsets[$ctr++][8]="Ovlap-30to70";
}
if ($count_OL3070_comp > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_OL3070_comp; $itemsets[$ctr][7]="C"; $itemsets[$ctr++][8]="Ovlap-30to70";
}
if ($count_OL3070_both > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_OL3070_both; $itemsets[$ctr][7]="B"; $itemsets[$ctr++][8]="Ovlap-30to70";
}
if ($count_OL70plus > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_OL70plus; $itemsets[$ctr][7]="+"; $itemsets[$ctr++][8]="Ovlap-70plus";
}
if ($count_OL70plus_comp > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_OL70plus_comp; $itemsets[$ctr][7]="C"; $itemsets[$ctr++][8]="Ovlap-70plus";
}
if ($count_OL70plus_both > 0){
	$itemsets[$ctr][0]=$ref_fam; $itemsets[$ctr][1]=&get_count($ref_fam); $itemsets[$ctr][2]=&get_avglen($ref_fam);
	$itemsets[$ctr][3]=$sub_fam; $itemsets[$ctr][4]=&get_count($sub_fam); $itemsets[$ctr][5]=&get_avglen($sub_fam);
	$itemsets[$ctr][6]=$count_OL70plus_both; $itemsets[$ctr][7]="B"; $itemsets[$ctr++][8]="Ovlap-70plus";
}

print STDERR '.';
# record tot itemsets
$tot_isets= $ctr;

# @itemsets : fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, Occurence, Strand, Category
#sorting @itemsets on occurences in decreasing order
@temp = sort {($a->[7] cmp $b->[7]) or ($b->[6] <=> $a->[6])} @itemsets;
@itemsets=@temp;

#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\# Runtime details after finding the itemsets: \n";
print OUTFILEDATA "\# System time for process: $system_t\n";
#print OUTFILEDATA "\# System time for children: $csystem_t\n";
print OUTFILEDATA "\# User time for process: $user_t\n";
#print OUTFILEDATA "\# User time for children: $cuser_t\n";


print STDERR "\n\# Runtime details after finding the $tot_isets itemsets: \n";
print STDERR "\# System time for process: $system_t\n";
print STDERR "\# User time for process: $user_t\n";


#PRUNING THE ITEMSET OF RECIPROCAL RELATIONS
$ctr=0;
# @itemsets : fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, Occurence, Strand, Category
# initializing the array
$pruned_itemsets[$ctr++]=$itemsets[0];

foreach $i (1 .. $#itemsets){#skipping the 1st record
	my ($flag,$ref_rel, $ref_fam1,$ref_fam2,$ref_occ,$sub_fam1,
	$sub_fam2,$sub_occ,$sub_rel);
	
	$ref_fam1=$itemsets[$i][0];
	$ref_fam2=$itemsets[$i][3];
	$ref_occ=$itemsets[$i][6];
	$ref_rel=$itemsets[$i][8];

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
	elsif ($ref_rel eq "u4") { $ref_rel="d4";}
	elsif ($ref_rel eq "u5") { $ref_rel="d5";}
	elsif ($ref_rel eq "d1") { $ref_rel="u1";}
	elsif ($ref_rel eq "d2") { $ref_rel="u2";}
	elsif ($ref_rel eq "d3") { $ref_rel="u3";}
	elsif ($ref_rel eq "d4") { $ref_rel="u4";}
	elsif ($ref_rel eq "d5") { $ref_rel="u5";}
	
	$flag=0;
	# filtering out overlap relations among images of the same family
	if (($ref_fam1 eq $ref_fam2) && ($ref_rel=~/^Ovlap/)){ 
		$flag=1;
	}
	elsif(($ref_rel=~/^u/) || ($ref_rel=~/^d/)){# to avoid processing Ovlap and IN itemsets
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
			$sub_rel=$j->[8];
			
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
	}
	
	if($flag == 0){
		$pruned_itemsets[$ctr++]=$itemsets[$i];
	}
	print STDERR '.';
}


#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\# Runtime details after pruning the itemsets: \n";
print OUTFILEDATA "\# System time for process: $system_t\n";
#print OUTFILEDATA "\# System time for children: $csystem_t\n";
print OUTFILEDATA "\# User time for process: $user_t\n";
#print OUTFILEDATA "\# User time for children: $cuser_t\n";

print STDERR "\n\# Runtime details after pruning the itemsets to $ctr: \n";
print STDERR "\# System time for process: $system_t\n";
print STDERR "\# User time for process: $user_t\n";

#printing the itemsets in tabbed format
# @itemsets : fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, Occurence, Category
print OUTFILEDATA "\n\# Total records in out file: $tot_recs\n";
print OUTFILEDATA "\# Total families: $tot_fams\n";
print OUTFILEDATA "\# Total relations: $tot_rels\n";
print OUTFILEDATA "\# Total itemsets: $tot_isets\n";
print OUTFILEDATA "\# Total itemsets after pruning: $ctr\n\n";
# print OUTFILEDATA "\# Runtime details: \n";
# print OUTFILEDATA "\# System time for process: $system_t\n";
# print OUTFILEDATA "\# System time for children: $csystem_t\n";
# print OUTFILEDATA "\# User time for process: $user_t\n";
# print OUTFILEDATA "\# User time for children: $cuser_t\n";
print OUTFILEDATA "\# Printing pruned \@itemsets...............\n";
print OUTFILEDATA "\# - Removed exact reciprocal relationships\n";
print OUTFILEDATA "\# - Removed same family overlap relationships\n";

print OUTFILEDATA "\#fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, Occurence, Strand, Category\n";

for $i (0 .. $#pruned_itemsets){
	print OUTFILEDATA $pruned_itemsets[$i][0],"\t",$pruned_itemsets[$i][1],"\t",
	$pruned_itemsets[$i][2],"\t",$pruned_itemsets[$i][3],"\t",$pruned_itemsets[$i][4],"\t",
	$pruned_itemsets[$i][5],"\t",$pruned_itemsets[$i][6],"\t",$pruned_itemsets[$i][7],
	"\t",$pruned_itemsets[$i][8],"\n";
}

#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\n\# Runtime details after printing the pruned itemsets: \n";
print OUTFILEDATA "\# System time for process: $system_t\n";
#print OUTFILEDATA "\# System time for children: $csystem_t\n";
print OUTFILEDATA "\# User time for process: $user_t\n";
#print OUTFILEDATA "\# User time for children: $cuser_t\n";

$i=localtime();
print OUTFILEDATA "\n\#Time: $i\n";

print STDERR "\n\# Runtime details after printing the pruned itemsets: \n";
print STDERR "\# System time for process: $system_t\n";
print STDERR "\# User time for process: $user_t\n";



close (INFILEDATA);
close (OUTFILEDATA);

exit;
