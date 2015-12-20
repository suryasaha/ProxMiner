#!/usr/bin/perl -w
# MGEL
# Surya Saha 3/15/07
# reading cmd line input .out file which is sorted on the start position
# and finds the relationship among images and families
# Relationship types: 
# Upstream: u1 (0-500 bases),u2 (500-1000 bases),u3 (1000-5000 bases), u4 (5000-10000 bases), u5 (10000-15000 bases)
# Downstream: d1 (0-500 bases),d2 (500-1000 bases),d3 (1000-5000 bases), d4 (5000-10000 bases), d5 (10000-15000 bases)
# In: Location of fam2 is entirely within fam1 (IN)
# Contains: Location of fam2 is entirely within fam1 (Cont)
# Overlap: 
# single linkage algo so consider overlap if > 10% of either
# 10% to 30% (Ovlap-10to30)
# 30% to 70% (Ovlap-30to70)
# 70% + (Ovlap>70)

# Creating the frequent itemsets in the format
# fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, Occurence, Strand, Category



# v3: Removed all duplicate counting
# v3: Counts all relationships 
# v4: Optimized the code to avoid recording itemsets with 0 count
# v4: Check for function call with large parameters

# v5: count relations for images INSTEAD of families
# v5: Use the strand information to calculate the relationships (See rLog)

# v6: Optimize the code (remove duplicates)
# v6: Fixed the bug where false relations were being counted for 'B' strand because of missing ELSE 

# v7: Better progress messages
# v7: Hangs with chr12.con.out

# v8 : 07/01/07
# v8: Reducing the number of loops
# v8: No pruning.
# v8 : F1 O F2 is equal to F2 O F1 if F1==F2
# v8 : Huge improvement in complexity (8+ hours to 36 mins for chr12.con.out)
# v8 : Both the sub_fam/img and ref_fam/img will not take part in relationships with ref_fam and 
#      sub_fam respec. in the future

# v9 : Added a reciprocal relationship for IN called CONTAINS to handle the new confidence calulation

# v10: Writing out the copy information for each relationship
# v10: Writing out the copy information for each relationship separately for each strand
# v10.1: Introduced a flag to prevent writing out copies file (for noise files)

use strict;
use warnings;
use POSIX;

unless (@ARGV == 2){
	print "USAGE: $0 <input .out file> <write copies??(0/1)>\n";
	exit;
}


my  ($ifname,$rec,@temp,%temphash,$ctr,$i,$j,$copy_file_flag);

my (@table,@famnames,@counts,$ups_ctr,$dns_ctr, $pos_copy_ctr, 
$comp_copy_ctr, $both_copy_ctr,$ref_img,$ref_fam,$ref_start,
$ref_end,$ref_strand,$sub_fam,$sub_img,$sub_start,$sub_end,$sub_strand,
%pos_relationships, %comp_relationships, %both_relationships, 
%pos_rel_history, %comp_rel_history, %both_rel_history, @pos_copies, 
@comp_copies, @both_copies,$tot_fams,$tot_recs,
$user_t,$system_t, $cuser_t,$csystem_t);

$ifname=$ARGV[0];
chomp $ifname;
$copy_file_flag=$ARGV[1];
chomp $copy_file_flag;
if(!($copy_file_flag == 0 || $copy_file_flag == 1) ){ print STDERR "flag can be only 0 or 1\nExiting..\n"; exit;}

unless(open(INFILEDATA,$ifname)){print "not able to open ".$ifname."\n\n";exit;}

unless(open(OUTFILEDATA,">$ifname.f_itemsets.tab")){print "not able to open ".$ifname."f_itemsets.tab\n\n";exit;}

if($copy_file_flag){
	unless(open(OUTFILECOPIESPOS,">$ifname.pos.copies.tab")){print "not able to open ".$ifname."pos.copies.tab \n\n";exit;}
	unless(open(OUTFILECOPIESCOMP,">$ifname.comp.copies.tab")){print "not able to open ".$ifname."comp.copies.tab \n\n";exit;}
	unless(open(OUTFILECOPIESBOTH,">$ifname.both.copies.tab")){print "not able to open ".$ifname."both.copies.tab \n\n";exit;}
}

# debugging
# unless(open(ERRFILE,">ERRFILE")){print "not able to open ERRFILE \n\n";exit;}

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

print OUTFILEDATA "\# Version: 10.1\n";
if($copy_file_flag){
	print OUTFILECOPIESPOS "\# Version: 10.1\n";
	print OUTFILECOPIESCOMP "\# Version: 10.1\n";
	print OUTFILECOPIESBOTH "\# Version: 10.1\n";
}
$i=localtime();
print OUTFILEDATA "\# Time: $i\n";
if($copy_file_flag){
	print OUTFILECOPIESPOS "\# Time: $i\n";
	print OUTFILECOPIESCOMP "\# Time: $i\n";
	print OUTFILECOPIESBOTH "\# Time: $i\n";
}

#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\# Runtime details after reading in the file: \n";
print OUTFILEDATA "\# System time for process: $system_t\n";
print OUTFILEDATA "\# User time for process: $user_t\n";
print OUTFILEDATA "\n";


print STDERR "\# Runtime details after reading $tot_recs from file: \n";
print STDERR "\# System time for process: $system_t\n";
print STDERR "\# User time for process: $user_t\n\n";


#@table
#1935 10.6  0.0  2.8 chr12 8936  9225 27748096 C  R=286 Unknown (0) 283   2
#0    1     2    3   4     5     6    7        8  9     10      11  12    13


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

# Add a field to end of @table 
# where @table[][14]=image number
foreach (@table){
	# since @counts[][3] is initialized to 0
	$_->[14]=1+$counts[&get_index($_->[9])][3]++; 
}


#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\# Runtime details after preparing \@counts and appending \@table: \n";
print OUTFILEDATA "\# System time for process: $system_t\n";
print OUTFILEDATA "\# User time for process: $user_t\n";
print OUTFILEDATA "\n";

print STDERR "\# Runtime details after preparing \@counts and appending \@table: \n";
print STDERR "\# System time for process: $system_t\n";
print STDERR "\# User time for process: $user_t\n";
print STDERR "\n";

# FINDING ALL RELATIONS
# @table sorted on start position
# 1935 10.6  0.0  2.8 chr12 8936  9225 27748096 C  R=286 Unknown (0) 283   2   3
# 0    1     2    3   4     5     6    7        8  9     10      11  12    13 14
# @count: fam occurences avg-len imagenum
# finding the relationships
# %pos_relationships : [fam1 fam2 category] = count
# %comp_relationships : [fam1 fam2 category] = count
# %both_relationships : [fam1 fam2 category] = count
# %pos_rel_history : [fam1 fam1-img category] = fam2
# %comp_rel_history : [fam1 fam1-img category] = fam2
# %both_rel_history : [fam1 fam1-img category] = fam2

$ups_ctr=$dns_ctr=$pos_copy_ctr=$comp_copy_ctr=$both_copy_ctr=0;
for $i (0 .. $#table){
	$ref_start=$table[$i][5]; $ref_end=$table[$i][6];
	$ref_strand=$table[$i][8]; $ref_fam=$table[$i][9];
	$ref_img=$table[$i][14];
	
	# cleaning up
	$ref_start=~ s/\s//g; $ref_end=~ s/\s//g;
	$ref_strand=~ s/\s//g; $ref_fam=~ s/\s//g;
	$ref_img=~ s/\s//g;

	print STDERR '.';
	$j=$i;
	
	# only look for relationships with images located before it and
	# ending within 15k bases before ref_start or anytime after it
	while(($j!=0) && ($table[$j-1][6] > $ref_start-15000)) {
		$ups_ctr++;
		$j--;
		$sub_start=$table[$j][5]; $sub_end=$table[$j][6];
		$sub_strand=$table[$j][8]; $sub_fam=$table[$j][9];
		$sub_img=$table[$j][14];
		
		
		# cleaning up
		$sub_start=~ s/\s//g; $sub_end=~ s/\s//g;
		$sub_strand=~ s/\s//g; $sub_fam=~ s/\s//g;
		$sub_img=~ s/\s//g;
		
		# Note: since all relationship are exclusive, I have used elsif
		
		# In: Location of ref fam is entirely within sub fam (IN,CONT)
		# IN should be first bcos if sub start is near the ref start, it will
		# be listed right before the ref record in the list
		if(($sub_start <= $ref_start) && ($sub_end >= $ref_end)){
			if ($ref_strand eq $sub_strand){
				if($ref_strand eq '+'){# pos strand
					# check if the ref image has this relation with this family already
					if ((!(exists $pos_rel_history{"$ref_fam $ref_img $sub_fam IN"}) || ($pos_rel_history{"$ref_fam $ref_img $sub_fam IN"} ne $sub_fam)) && (!(exists $pos_rel_history{"$sub_fam $sub_img $ref_fam IN"}) || ($pos_rel_history{"$sub_fam $sub_img $ref_fam IN"} ne $ref_fam))) {
						# create history entry
						$pos_rel_history{"$ref_fam $ref_img $sub_fam IN"} =$sub_fam;
						$pos_rel_history{"$sub_fam $sub_img $ref_fam IN"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $pos_relationships{"$ref_fam $sub_fam IN"}) {
							$pos_relationships{"$ref_fam $sub_fam IN"} = 1;
						}
						else{
							$pos_relationships{"$ref_fam $sub_fam IN"}++;
						}
						
						# add record for IN relationship
						$pos_copies[$pos_copy_ctr][0]=$ref_fam; $pos_copies[$pos_copy_ctr][1]="IN";
						$pos_copies[$pos_copy_ctr][2]=$sub_fam; $pos_copies[$pos_copy_ctr][3]='+';
						$pos_copies[$pos_copy_ctr][4]=$ref_start; $pos_copies[$pos_copy_ctr][5]=$ref_end;
						$pos_copies[$pos_copy_ctr][6]=$sub_start; $pos_copies[$pos_copy_ctr++][7]=$sub_end;

						# increment reciprocal relationship count or create relationship entry
						if (!exists $pos_relationships{"$sub_fam $ref_fam CONT"}) {
							$pos_relationships{"$sub_fam $ref_fam CONT"} = 1;
						}
						else{
							$pos_relationships{"$sub_fam $ref_fam CONT"}++;
						}
						
						# add record for CONT relationship
						$pos_copies[$pos_copy_ctr][0]=$sub_fam; $pos_copies[$pos_copy_ctr][1]="CONT";
						$pos_copies[$pos_copy_ctr][2]=$ref_fam; $pos_copies[$pos_copy_ctr][3]='+';
						$pos_copies[$pos_copy_ctr][4]=$sub_start; $pos_copies[$pos_copy_ctr][5]=$sub_end;
						$pos_copies[$pos_copy_ctr][6]=$ref_start; $pos_copies[$pos_copy_ctr++][7]=$ref_end;
					}
				}
				elsif($ref_strand eq 'C'){# comp strand
					# check if the ref image has this relation with this family already
					if ((!(exists $comp_rel_history{"$ref_fam $ref_img $sub_fam IN"}) || ($comp_rel_history{"$ref_fam $ref_img $sub_fam IN"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_img $ref_fam IN"}) || ($comp_rel_history{"$sub_fam $sub_img $ref_fam IN"} ne $ref_fam))) {
						# create history entry
						$comp_rel_history{"$ref_fam $ref_img $sub_fam IN"} = $sub_fam;
						$comp_rel_history{"$sub_fam $sub_img $ref_fam IN"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $comp_relationships{"$ref_fam $sub_fam IN"}) {
							$comp_relationships{"$ref_fam $sub_fam IN"} = 1;
						}
						else{
							$comp_relationships{"$ref_fam $sub_fam IN"}++;
						}
						
						# add record for IN relationship
						$comp_copies[$comp_copy_ctr][0]=$ref_fam; $comp_copies[$comp_copy_ctr][1]="IN";
						$comp_copies[$comp_copy_ctr][2]=$sub_fam; $comp_copies[$comp_copy_ctr][3]='C';
						$comp_copies[$comp_copy_ctr][4]=$ref_start; $comp_copies[$comp_copy_ctr][5]=$ref_end;
						$comp_copies[$comp_copy_ctr][6]=$sub_start; $comp_copies[$comp_copy_ctr++][7]=$sub_end;
						
						# increment reciprocal relationship count or create relationship entry
						if (!exists $comp_relationships{"$sub_fam $ref_fam CONT"}) {
							$comp_relationships{"$sub_fam $ref_fam CONT"} = 1;
						}
						else{
							$comp_relationships{"$sub_fam $ref_fam CONT"}++;
						}
						
						# add record for CONT relationship
						$comp_copies[$comp_copy_ctr][0]=$sub_fam; $comp_copies[$comp_copy_ctr][1]="CONT";
						$comp_copies[$comp_copy_ctr][2]=$ref_fam; $comp_copies[$comp_copy_ctr][3]='C';
						$comp_copies[$comp_copy_ctr][4]=$sub_start; $comp_copies[$comp_copy_ctr][5]=$sub_end;
						$comp_copies[$comp_copy_ctr][6]=$ref_start; $comp_copies[$comp_copy_ctr++][7]=$ref_end;
					}
				}
			}
			# irrespective of strand
			# check if the ref image has this relation with this family already
			if ((!(exists $both_rel_history{"$ref_fam $ref_img $sub_fam IN"}) || ($both_rel_history{"$ref_fam $ref_img $sub_fam IN"} ne $sub_fam)) && (!(exists $both_rel_history{"$sub_fam $sub_img $ref_fam IN"}) || ($both_rel_history{"$sub_fam $sub_img $ref_fam IN"} ne $ref_fam))) {
				# create history entry
				$both_rel_history{"$ref_fam $ref_img $sub_fam IN"} = $sub_fam;
				$both_rel_history{"$sub_fam $sub_img $ref_fam IN"} = $ref_fam;
				
				# increment relationship count or create relationship entry
				if (!exists $both_relationships{"$ref_fam $sub_fam IN"}) {
					$both_relationships{"$ref_fam $sub_fam IN"} = 1;
				}
				else{
					$both_relationships{"$ref_fam $sub_fam IN"}++;
				}
				
				# add record for IN relationship
				$both_copies[$both_copy_ctr][0]=$ref_fam; $both_copies[$both_copy_ctr][1]="IN";
				$both_copies[$both_copy_ctr][2]=$sub_fam; $both_copies[$both_copy_ctr][3]='B';
				$both_copies[$both_copy_ctr][4]=$ref_start; $both_copies[$both_copy_ctr][5]=$ref_end;
				$both_copies[$both_copy_ctr][6]=$sub_start; $both_copies[$both_copy_ctr++][7]=$sub_end;

				# increment reciprocal relationship count or create relationship entry
				if (!exists $both_relationships{"$sub_fam $ref_fam CONT"}) {
					$both_relationships{"$sub_fam $ref_fam CONT"} = 1;
				}
				else{
					$both_relationships{"$sub_fam $ref_fam CONT"}++;
				}
				
				# add record for CONT relationship
				$both_copies[$both_copy_ctr][0]=$sub_fam; $both_copies[$both_copy_ctr][1]="CONT";
				$both_copies[$both_copy_ctr][2]=$ref_fam; $both_copies[$both_copy_ctr][3]='B';
				$both_copies[$both_copy_ctr][4]=$sub_start; $both_copies[$both_copy_ctr][5]=$sub_end;
				$both_copies[$both_copy_ctr][6]=$ref_start; $both_copies[$both_copy_ctr++][7]=$ref_end;
			}
		}# IN end
		
		# Overlap: If overlap is more than 10% of length of either family (Ovlap)
		# now if subject fam ends within the reference fam
		elsif (($sub_end > $ref_start) &&  ($sub_end < $ref_end)) {
			my ($ovlap, $ref_ovlap, $sub_ovlap);
			$ovlap = $sub_end - $ref_start;
			
			$ref_ovlap = ($ovlap / ($ref_end - $ref_start)) * 100;
			$sub_ovlap = ($ovlap / ($sub_end - $sub_start)) * 100;
			
			# Overlap :10% to 30% (Ovlap-10to30)
			if ((($ref_ovlap > 10.00) && ($ref_ovlap <= 30.00)) || 
			(($sub_ovlap > 10.00) && ($sub_ovlap <= 30.00))) {
				if ($ref_strand eq $sub_strand){
					if($ref_strand eq '+'){# pos strand
						# check if the ref image has this relation with this family already
						if ((!(exists $pos_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-10to30"}) || ($pos_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-10to30"} ne $sub_fam)) && (!(exists $pos_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-10to30"}) || ($pos_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-10to30"} ne $ref_fam))) {
							# create history entry
							$pos_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-10to30"} = $sub_fam;
							$pos_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-10to30"} = $ref_fam;
							
							# increment relationship count or create relationship entry
							if (!exists $pos_relationships{"$sub_fam $ref_fam Ovlap-10to30"}) {
								$pos_relationships{"$sub_fam $ref_fam Ovlap-10to30"} = 1;
							}
							else{
								$pos_relationships{"$sub_fam $ref_fam Ovlap-10to30"}++;
							}
							
							# add record for Ovlap-10to30 relationship
							$pos_copies[$pos_copy_ctr][0]=$sub_fam; $pos_copies[$pos_copy_ctr][1]="Ovlap-10to30";
							$pos_copies[$pos_copy_ctr][2]=$ref_fam; $pos_copies[$pos_copy_ctr][3]='+';
							$pos_copies[$pos_copy_ctr][4]=$sub_start; $pos_copies[$pos_copy_ctr][5]=$sub_end;
							$pos_copies[$pos_copy_ctr][6]=$ref_start; $pos_copies[$pos_copy_ctr++][7]=$ref_end;
						}
					}
					elsif($ref_strand eq 'C'){# comp strand
						# check if the ref image has this relation with this family already
						if ((!(exists $comp_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-10to30"}) || ($comp_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-10to30"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-10to30"}) || ($comp_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-10to30"} ne $ref_fam))) {
							# create history entry
							$comp_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-10to30"} = $sub_fam;
							$comp_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-10to30"} = $ref_fam;
							
							# increment relationship count or create relationship entry
							if (!exists $comp_relationships{"$sub_fam $ref_fam Ovlap-10to30"}) {
								$comp_relationships{"$sub_fam $ref_fam Ovlap-10to30"} = 1;
							}
							else{
								$comp_relationships{"$sub_fam $ref_fam Ovlap-10to30"}++;
							}
							
							# add record for Ovlap-10to30 relationship
							$comp_copies[$comp_copy_ctr][0]=$sub_fam; $comp_copies[$comp_copy_ctr][1]="Ovlap-10to30";
							$comp_copies[$comp_copy_ctr][2]=$ref_fam; $comp_copies[$comp_copy_ctr][3]='C';
							$comp_copies[$comp_copy_ctr][4]=$sub_start; $comp_copies[$comp_copy_ctr][5]=$sub_end;
							$comp_copies[$comp_copy_ctr][6]=$ref_start; $comp_copies[$comp_copy_ctr++][7]=$ref_end;
						}
					}
				}
				# irrespective of strand
				# check if the ref image has this relation with this family already
				if ((!(exists $both_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-10to30"}) || ($both_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-10to30"} ne $sub_fam)) && (!(exists $both_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-10to30"}) || ($both_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-10to30"} ne $ref_fam))) {
					# create history entry
					$both_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-10to30"} = $sub_fam;
					$both_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-10to30"} = $ref_fam;

					
					# increment relationship count or create relationship entry
					if (!exists $both_relationships{"$sub_fam $ref_fam Ovlap-10to30"}) {
						$both_relationships{"$sub_fam $ref_fam Ovlap-10to30"} = 1;
					}
					else{
						$both_relationships{"$sub_fam $ref_fam Ovlap-10to30"}++;
					}
					
					# add record for Ovlap-10to30 relationship
					$both_copies[$both_copy_ctr][0]=$sub_fam; $both_copies[$both_copy_ctr][1]="Ovlap-10to30";
					$both_copies[$both_copy_ctr][2]=$ref_fam; $both_copies[$both_copy_ctr][3]='B';
					$both_copies[$both_copy_ctr][4]=$sub_start; $both_copies[$both_copy_ctr][5]=$sub_end;
					$both_copies[$both_copy_ctr][6]=$ref_start; $both_copies[$both_copy_ctr++][7]=$ref_end;
				}
			}
			# Overlap :30% to 70% (Ovlap-30to70)
			elsif ((($ref_ovlap > 30.00) && ($ref_ovlap <= 70.00)) || 
			(($sub_ovlap > 30.00) && ($sub_ovlap <= 70.00))) {
				if ($ref_strand eq $sub_strand){
					if($ref_strand eq '+'){# pos strand
						# check if the ref image has this relation with this family already
						if ((!(exists $pos_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-30to70"}) || ($pos_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-30to70"} ne $sub_fam)) && (!(exists $pos_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-30to70"}) || ($pos_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-30to70"} ne $ref_fam))) {
							# create history entry
							$pos_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-30to70"} = $sub_fam;
							$pos_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-30to70"} = $ref_fam;
							
							# increment relationship count or create relationship entry
							if (!exists $pos_relationships{"$sub_fam $ref_fam Ovlap-30to70"}) {
								$pos_relationships{"$sub_fam $ref_fam Ovlap-30to70"} = 1;
							}
							else{
								$pos_relationships{"$sub_fam $ref_fam Ovlap-30to70"}++;
							}
							
							# add record for Ovlap-30to70 relationship
							$pos_copies[$pos_copy_ctr][0]=$sub_fam; $pos_copies[$pos_copy_ctr][1]="Ovlap-30to70";
							$pos_copies[$pos_copy_ctr][2]=$ref_fam; $pos_copies[$pos_copy_ctr][3]='+';
							$pos_copies[$pos_copy_ctr][4]=$sub_start; $pos_copies[$pos_copy_ctr][5]=$sub_end;
							$pos_copies[$pos_copy_ctr][6]=$ref_start; $pos_copies[$pos_copy_ctr++][7]=$ref_end;
						}
					}
					elsif($ref_strand eq 'C'){# comp strand
						# check if the ref image has this relation with this family already
						if ((!(exists $comp_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-30to70"}) || ($comp_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-30to70"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-30to70"}) || ($comp_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-30to70"} ne $ref_fam))) {
							# create history entry
							$comp_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-30to70"} = $sub_fam;
							$comp_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-30to70"} = $ref_fam;
							
							# increment relationship count or create relationship entry
							if (!exists $comp_relationships{"$sub_fam $ref_fam Ovlap-30to70"}) {
								$comp_relationships{"$sub_fam $ref_fam Ovlap-30to70"} = 1;
							}
							else{
								$comp_relationships{"$sub_fam $ref_fam Ovlap-30to70"}++;
							}
							
							# add record for Ovlap-30to70 relationship
							$comp_copies[$comp_copy_ctr][0]=$sub_fam; $comp_copies[$comp_copy_ctr][1]="Ovlap-30to70";
							$comp_copies[$comp_copy_ctr][2]=$ref_fam; $comp_copies[$comp_copy_ctr][3]='C';
							$comp_copies[$comp_copy_ctr][4]=$sub_start; $comp_copies[$comp_copy_ctr][5]=$sub_end;
							$comp_copies[$comp_copy_ctr][6]=$ref_start; $comp_copies[$comp_copy_ctr++][7]=$ref_end;
						}
					}
				}
				# irrespective of strand
				# check if the ref image has this relation with this family already
				if ((!(exists $both_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-30to70"}) || ($both_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-30to70"} ne $sub_fam)) && (!(exists $both_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-30to70"}) || ($both_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-30to70"} ne $ref_fam))) {
					# create history entry
					$both_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-30to70"} = $sub_fam;
					$both_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-30to70"} = $ref_fam;
					
					# increment relationship count or create relationship entry
					if (!exists $both_relationships{"$sub_fam $ref_fam Ovlap-30to70"}) {
						$both_relationships{"$sub_fam $ref_fam Ovlap-30to70"} = 1;
					}
					else{
						$both_relationships{"$sub_fam $ref_fam Ovlap-30to70"}++;
					}
					
					# add record for Ovlap-30to70 relationship
					$both_copies[$both_copy_ctr][0]=$sub_fam; $both_copies[$both_copy_ctr][1]="Ovlap-30to70";
					$both_copies[$both_copy_ctr][2]=$ref_fam; $both_copies[$both_copy_ctr][3]='B';
					$both_copies[$both_copy_ctr][4]=$sub_start; $both_copies[$both_copy_ctr][5]=$sub_end;
					$both_copies[$both_copy_ctr][6]=$ref_start; $both_copies[$both_copy_ctr++][7]=$ref_end;
				}
			}
			# Overlap : >70% (Ovlap-70plus)
			elsif (($ref_ovlap > 70.00) || ($sub_ovlap > 70.00)) {
				if ($ref_strand eq $sub_strand){
					if($ref_strand eq '+'){# pos strand
						# check if the ref image has this relation with this family already
						if ((!(exists $pos_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-70plus"}) || ($pos_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-70plus"} ne $sub_fam)) && (!(exists $pos_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-70plus"}) || ($pos_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-70plus"} ne $ref_fam))) {
							# create history entry
							$pos_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-70plus"} = $sub_fam;
							$pos_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-70plus"} = $ref_fam;
							
							# increment relationship count or create relationship entry
							if (!exists $pos_relationships{"$sub_fam $ref_fam Ovlap-70plus"}) {
								$pos_relationships{"$sub_fam $ref_fam Ovlap-70plus"} = 1;
							}
							else{
								$pos_relationships{"$sub_fam $ref_fam Ovlap-70plus"}++;
							}
							
							# add record for Ovlap-70plus relationship
							$pos_copies[$pos_copy_ctr][0]=$sub_fam; $pos_copies[$pos_copy_ctr][1]="Ovlap-70plus";
							$pos_copies[$pos_copy_ctr][2]=$ref_fam; $pos_copies[$pos_copy_ctr][3]='+';
							$pos_copies[$pos_copy_ctr][4]=$sub_start; $pos_copies[$pos_copy_ctr][5]=$sub_end;
							$pos_copies[$pos_copy_ctr][6]=$ref_start; $pos_copies[$pos_copy_ctr++][7]=$ref_end;
						}
					}
					elsif($ref_strand eq 'C'){# comp strand
						# check if the ref image has this relation with this family already
						if ((!(exists $comp_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-70plus"}) || ($comp_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-70plus"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-70plus"}) || ($comp_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-70plus"} ne $ref_fam))) {
							# create history entry
							$comp_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-70plus"} = $sub_fam;
							$comp_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-70plus"} = $ref_fam;
							
							# increment relationship count or create relationship entry
							if (!exists $comp_relationships{"$sub_fam $ref_fam Ovlap-70plus"}) {
								$comp_relationships{"$sub_fam $ref_fam Ovlap-70plus"} = 1;
							}
							else{
								$comp_relationships{"$sub_fam $ref_fam Ovlap-70plus"}++;
							}
							
							# add record for Ovlap-70plus relationship
							$comp_copies[$comp_copy_ctr][0]=$sub_fam; $comp_copies[$comp_copy_ctr][1]="Ovlap-70plus";
							$comp_copies[$comp_copy_ctr][2]=$ref_fam; $comp_copies[$comp_copy_ctr][3]='C';
							$comp_copies[$comp_copy_ctr][4]=$sub_start; $comp_copies[$comp_copy_ctr][5]=$sub_end;
							$comp_copies[$comp_copy_ctr][6]=$ref_start; $comp_copies[$comp_copy_ctr++][7]=$ref_end;
						}
					}
				}
				# irrespective of strand
				# check if the ref image has this relation with this family already
				if ((!(exists $both_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-70plus"}) || ($both_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-70plus"} ne $sub_fam)) && (!(exists $both_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-70plus"}) || ($both_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-70plus"} ne $ref_fam))) {
					# create history entry
					$both_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-70plus"} = $sub_fam;
					$both_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-70plus"} = $ref_fam;

					
					# increment relationship count or create relationship entry
					if (!exists $both_relationships{"$sub_fam $ref_fam Ovlap-70plus"}) {
						$both_relationships{"$sub_fam $ref_fam Ovlap-70plus"} = 1;
					}
					else{
						$both_relationships{"$sub_fam $ref_fam Ovlap-70plus"}++;
					}
					
					# add record for Ovlap-70plus relationship
					$both_copies[$both_copy_ctr][0]=$sub_fam; $both_copies[$both_copy_ctr][1]="Ovlap-70plus";
					$both_copies[$both_copy_ctr][2]=$ref_fam; $both_copies[$both_copy_ctr][3]='B';
					$both_copies[$both_copy_ctr][4]=$sub_start; $both_copies[$both_copy_ctr][5]=$sub_end;
					$both_copies[$both_copy_ctr][6]=$ref_start; $both_copies[$both_copy_ctr++][7]=$ref_end;
				}
			}
		}# overlap end
		# Upstream: u1 (0-500 bases)
		elsif(($sub_end <= $ref_start) && ($sub_end > $ref_start-500)){
			if ($ref_strand eq $sub_strand){
				if($ref_strand eq '+'){# pos strand
					# check if the ref image has this relation with this family already
					if ((!(exists $pos_rel_history{"$ref_fam $ref_img $sub_fam u1"}) || ($pos_rel_history{"$ref_fam $ref_img $sub_fam u1"} ne $sub_fam)) && (!(exists $pos_rel_history{"$sub_fam $sub_img $ref_fam u1"}) || ($pos_rel_history{"$sub_fam $sub_img $ref_fam u1"} ne $ref_fam))) {
						# create history entry
						$pos_rel_history{"$ref_fam $ref_img $sub_fam u1"} = $sub_fam;
						$pos_rel_history{"$sub_fam $sub_img $ref_fam u1"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $pos_relationships{"$sub_fam $ref_fam u1"}) {
							$pos_relationships{"$sub_fam $ref_fam u1"} = 1;
						}
						else{
							$pos_relationships{"$sub_fam $ref_fam u1"}++;
						}
						# add record for u1 relationship
						$pos_copies[$pos_copy_ctr][0]=$sub_fam; $pos_copies[$pos_copy_ctr][1]="u1";
						$pos_copies[$pos_copy_ctr][2]=$ref_fam; $pos_copies[$pos_copy_ctr][3]='+';
						$pos_copies[$pos_copy_ctr][4]=$sub_start; $pos_copies[$pos_copy_ctr][5]=$sub_end;
						$pos_copies[$pos_copy_ctr][6]=$ref_start; $pos_copies[$pos_copy_ctr++][7]=$ref_end;
					}
				}
				elsif($ref_strand eq 'C'){# comp strand
					# check if the ref image has this relation with this family already
					if ((!(exists $comp_rel_history{"$ref_fam $ref_img $sub_fam u1"}) || ($comp_rel_history{"$ref_fam $ref_img $sub_fam u1"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_img $ref_fam u1"}) || ($comp_rel_history{"$sub_fam $sub_img $ref_fam u1"} ne $ref_fam))) {
						# create history entry
						$comp_rel_history{"$ref_fam $ref_img $sub_fam u1"} = $sub_fam;
						$comp_rel_history{"$sub_fam $sub_img $ref_fam u1"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $comp_relationships{"$ref_fam $sub_fam u1"}) {
							# now ref fam is upstream of sub fam as we are 
							# counting from right
							$comp_relationships{"$ref_fam $sub_fam u1"} = 1;
						}
						else{
							$comp_relationships{"$ref_fam $sub_fam u1"}++;
						}
						
						# add record for u1 relationship
						$comp_copies[$comp_copy_ctr][0]=$ref_fam; $comp_copies[$comp_copy_ctr][1]="u1";
						$comp_copies[$comp_copy_ctr][2]=$sub_fam; $comp_copies[$comp_copy_ctr][3]='C';
						$comp_copies[$comp_copy_ctr][4]=$ref_start; $comp_copies[$comp_copy_ctr][5]=$ref_end;
						$comp_copies[$comp_copy_ctr][6]=$sub_start; $comp_copies[$comp_copy_ctr++][7]=$sub_end;
					}
				}
			}
			# irrespective of strand
			# check if the ref image has this relation with this family already
			if ((!(exists $both_rel_history{"$ref_fam $ref_img $sub_fam u1"}) || ($both_rel_history{"$ref_fam $ref_img $sub_fam u1"} ne $sub_fam)) && (!(exists $both_rel_history{"$sub_fam $sub_img $ref_fam u1"}) || ($both_rel_history{"$sub_fam $sub_img $ref_fam u1"} ne $ref_fam))) {
				# create history entry
				$both_rel_history{"$ref_fam $ref_img $sub_fam u1"} = $sub_fam;
				$both_rel_history{"$sub_fam $sub_img $ref_fam u1"} = $ref_fam;
				
				# increment relationship count or create relationship entry
				if (!exists $both_relationships{"$sub_fam $ref_fam u1"}) {
					$both_relationships{"$sub_fam $ref_fam u1"} = 1;
				}
				else{
					$both_relationships{"$sub_fam $ref_fam u1"}++;
				}
				
				# add record for u1 relationship
				$both_copies[$both_copy_ctr][0]=$sub_fam; $both_copies[$both_copy_ctr][1]="u1";
				$both_copies[$both_copy_ctr][2]=$ref_fam; $both_copies[$both_copy_ctr][3]='B';
				$both_copies[$both_copy_ctr][4]=$sub_start; $both_copies[$both_copy_ctr][5]=$sub_end;
				$both_copies[$both_copy_ctr][6]=$ref_start; $both_copies[$both_copy_ctr++][7]=$ref_end;
			}
		}
		# Upstream: u2 (500-1000 bases)
		elsif(($sub_end <= $ref_start-500) && ($sub_end > $ref_start-1000)){
			if ($ref_strand eq $sub_strand){
				if($ref_strand eq '+'){# pos strand
					# check if the ref image has this relation with this family already
					if ((!(exists $pos_rel_history{"$ref_fam $ref_img $sub_fam u2"}) || ($pos_rel_history{"$ref_fam $ref_img $sub_fam u2"} ne $sub_fam)) && (!(exists $pos_rel_history{"$sub_fam $sub_img $ref_fam u2"}) || ($pos_rel_history{"$sub_fam $sub_img $ref_fam u2"} ne $ref_fam))) {
						# create history entry
						$pos_rel_history{"$ref_fam $ref_img $sub_fam u2"} = $sub_fam;
						$pos_rel_history{"$sub_fam $sub_img $ref_fam u2"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $pos_relationships{"$sub_fam $ref_fam u2"}) {
							$pos_relationships{"$sub_fam $ref_fam u2"} = 1;
						}
						else{
							$pos_relationships{"$sub_fam $ref_fam u2"}++;
						}
						# add record for u2 relationship
						$pos_copies[$pos_copy_ctr][0]=$sub_fam; $pos_copies[$pos_copy_ctr][1]="u2";
						$pos_copies[$pos_copy_ctr][2]=$ref_fam; $pos_copies[$pos_copy_ctr][3]='+';
						$pos_copies[$pos_copy_ctr][4]=$sub_start; $pos_copies[$pos_copy_ctr][5]=$sub_end;
						$pos_copies[$pos_copy_ctr][6]=$ref_start; $pos_copies[$pos_copy_ctr++][7]=$ref_end;
					}
				}
				elsif($ref_strand eq 'C'){# comp strand
					# check if the ref image has this relation with this family already
					if ((!(exists $comp_rel_history{"$ref_fam $ref_img $sub_fam u2"}) || ($comp_rel_history{"$ref_fam $ref_img $sub_fam u2"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_img $ref_fam u2"}) || ($comp_rel_history{"$sub_fam $sub_img $ref_fam u2"} ne $ref_fam))) {
						# create history entry
						$comp_rel_history{"$ref_fam $ref_img $sub_fam u2"} = $sub_fam;
						$comp_rel_history{"$sub_fam $sub_img $ref_fam u2"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $comp_relationships{"$ref_fam $sub_fam u2"}) {
							# now ref fam is upstream of sub fam as we are 
							# counting from right
							$comp_relationships{"$ref_fam $sub_fam u2"} = 1;
						}
						else{
							$comp_relationships{"$ref_fam $sub_fam u2"}++;
						}
						
						# add record for u2 relationship
						$comp_copies[$comp_copy_ctr][0]=$ref_fam; $comp_copies[$comp_copy_ctr][1]="u2";
						$comp_copies[$comp_copy_ctr][2]=$sub_fam; $comp_copies[$comp_copy_ctr][3]='C';
						$comp_copies[$comp_copy_ctr][4]=$ref_start; $comp_copies[$comp_copy_ctr][5]=$ref_end;
						$comp_copies[$comp_copy_ctr][6]=$sub_start; $comp_copies[$comp_copy_ctr++][7]=$sub_end;
					}
				}
			}
			# irrespective of strand
			# check if the ref image has this relation with this family already
			if ((!(exists $both_rel_history{"$ref_fam $ref_img $sub_fam u2"}) || ($both_rel_history{"$ref_fam $ref_img $sub_fam u2"} ne $sub_fam)) && (!(exists $both_rel_history{"$sub_fam $sub_img $ref_fam u2"}) || ($both_rel_history{"$sub_fam $sub_img $ref_fam u2"} ne $ref_fam))) {
				# create history entry
				$both_rel_history{"$ref_fam $ref_img $sub_fam u2"} = $sub_fam;
				$both_rel_history{"$sub_fam $sub_img $ref_fam u2"} = $ref_fam;
				
				# increment relationship count or create relationship entry
				if (!exists $both_relationships{"$sub_fam $ref_fam u2"}) {
					$both_relationships{"$sub_fam $ref_fam u2"} = 1;
				}
				else{
					$both_relationships{"$sub_fam $ref_fam u2"}++;
				}
				# add record for u2 relationship
				$both_copies[$both_copy_ctr][0]=$sub_fam; $both_copies[$both_copy_ctr][1]="u2";
				$both_copies[$both_copy_ctr][2]=$ref_fam; $both_copies[$both_copy_ctr][3]='B';
				$both_copies[$both_copy_ctr][4]=$sub_start; $both_copies[$both_copy_ctr][5]=$sub_end;
				$both_copies[$both_copy_ctr][6]=$ref_start; $both_copies[$both_copy_ctr++][7]=$ref_end;
			}
		}
		# Upstream: u3 (1000-5000 bases)
		elsif(($sub_end <= $ref_start-1000) && ($sub_end > $ref_start-5000)){
			if ($ref_strand eq $sub_strand){
				if($ref_strand eq '+'){# pos strand
					# check if the ref image has this relation with this family already
					if ((!(exists $pos_rel_history{"$ref_fam $ref_img $sub_fam u3"}) || ($pos_rel_history{"$ref_fam $ref_img $sub_fam u3"} ne $sub_fam)) && (!(exists $pos_rel_history{"$sub_fam $sub_img $ref_fam u3"}) || ($pos_rel_history{"$sub_fam $sub_img $ref_fam u3"} ne $ref_fam))) {
						# create history entry
						$pos_rel_history{"$ref_fam $ref_img $sub_fam u3"} = $sub_fam;
						$pos_rel_history{"$sub_fam $sub_img $ref_fam u3"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $pos_relationships{"$sub_fam $ref_fam u3"}) {
							$pos_relationships{"$sub_fam $ref_fam u3"} = 1;
						}
						else{
							$pos_relationships{"$sub_fam $ref_fam u3"}++;
						}
						# add record for u3 relationship
						$pos_copies[$pos_copy_ctr][0]=$sub_fam; $pos_copies[$pos_copy_ctr][1]="u3";
						$pos_copies[$pos_copy_ctr][2]=$ref_fam; $pos_copies[$pos_copy_ctr][3]='+';
						$pos_copies[$pos_copy_ctr][4]=$sub_start; $pos_copies[$pos_copy_ctr][5]=$sub_end;
						$pos_copies[$pos_copy_ctr][6]=$ref_start; $pos_copies[$pos_copy_ctr++][7]=$ref_end;
					}
				}
				elsif($ref_strand eq 'C'){# comp strand
					# check if the ref image has this relation with this family already
					if ((!(exists $comp_rel_history{"$ref_fam $ref_img $sub_fam u3"}) || ($comp_rel_history{"$ref_fam $ref_img $sub_fam u3"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_img $ref_fam u3"}) || ($comp_rel_history{"$sub_fam $sub_img $ref_fam u3"} ne $ref_fam))) {
						# create history entry
						$comp_rel_history{"$ref_fam $ref_img $sub_fam u3"} = $sub_fam;
						$comp_rel_history{"$sub_fam $sub_img $ref_fam u3"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $comp_relationships{"$ref_fam $sub_fam u3"}) {
							# now ref fam is upstream of sub fam as we are 
							# counting from right
							$comp_relationships{"$ref_fam $sub_fam u3"} = 1;
						}
						else{
							$comp_relationships{"$ref_fam $sub_fam u3"}++;
						}
						# add record for u3 relationship
						$comp_copies[$comp_copy_ctr][0]=$ref_fam; $comp_copies[$comp_copy_ctr][1]="u3";
						$comp_copies[$comp_copy_ctr][2]=$sub_fam; $comp_copies[$comp_copy_ctr][3]='C';
						$comp_copies[$comp_copy_ctr][4]=$ref_start; $comp_copies[$comp_copy_ctr][5]=$ref_end;
						$comp_copies[$comp_copy_ctr][6]=$sub_start; $comp_copies[$comp_copy_ctr++][7]=$sub_end;
					}
				}
			}
			# irrespective of strand
			# check if the ref image has this relation with this family already
			if ((!(exists $both_rel_history{"$ref_fam $ref_img $sub_fam u3"}) || ($both_rel_history{"$ref_fam $ref_img $sub_fam u3"} ne $sub_fam)) && (!(exists $both_rel_history{"$sub_fam $sub_img $ref_fam u3"}) || ($both_rel_history{"$sub_fam $sub_img $ref_fam u3"} ne $ref_fam))) {
				# create history entry
				$both_rel_history{"$ref_fam $ref_img $sub_fam u3"} = $sub_fam;
				$both_rel_history{"$sub_fam $sub_img $ref_fam u3"} = $ref_fam;
				
				# increment relationship count or create relationship entry
				if (!exists $both_relationships{"$sub_fam $ref_fam u3"}) {
					$both_relationships{"$sub_fam $ref_fam u3"} = 1;
				}
				else{
					$both_relationships{"$sub_fam $ref_fam u3"}++;
				}
				# add record for u3 relationship
				$both_copies[$both_copy_ctr][0]=$sub_fam; $both_copies[$both_copy_ctr][1]="u3";
				$both_copies[$both_copy_ctr][2]=$ref_fam; $both_copies[$both_copy_ctr][3]='B';
				$both_copies[$both_copy_ctr][4]=$sub_start; $both_copies[$both_copy_ctr][5]=$sub_end;
				$both_copies[$both_copy_ctr][6]=$ref_start; $both_copies[$both_copy_ctr++][7]=$ref_end;
			}
		}
		# Upstream: u4 (5000-10000 bases)
		elsif(($sub_end <= $ref_start-5000) && ($sub_end > $ref_start-10000)){
			if ($ref_strand eq $sub_strand){
				if($ref_strand eq '+'){# pos strand
					# check if the ref image has this relation with this family already
					if ((!(exists $pos_rel_history{"$ref_fam $ref_img $sub_fam u4"}) || ($pos_rel_history{"$ref_fam $ref_img $sub_fam u4"} ne $sub_fam)) && (!(exists $pos_rel_history{"$sub_fam $sub_img $ref_fam u4"}) || ($pos_rel_history{"$sub_fam $sub_img $ref_fam u4"} ne $ref_fam))) {
						# create history entry
						$pos_rel_history{"$ref_fam $ref_img $sub_fam u4"} = $sub_fam;
						$pos_rel_history{"$sub_fam $sub_img $ref_fam u4"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $pos_relationships{"$sub_fam $ref_fam u4"}) {
							$pos_relationships{"$sub_fam $ref_fam u4"} = 1;
						}
						else{
							$pos_relationships{"$sub_fam $ref_fam u4"}++;
						}
						# add record for u4 relationship
						$pos_copies[$pos_copy_ctr][0]=$sub_fam; $pos_copies[$pos_copy_ctr][1]="u4";
						$pos_copies[$pos_copy_ctr][2]=$ref_fam; $pos_copies[$pos_copy_ctr][3]='+';
						$pos_copies[$pos_copy_ctr][4]=$sub_start; $pos_copies[$pos_copy_ctr][5]=$sub_end;
						$pos_copies[$pos_copy_ctr][6]=$ref_start; $pos_copies[$pos_copy_ctr++][7]=$ref_end;
					}
				}
				elsif($ref_strand eq 'C'){# comp strand
					# check if the ref image has this relation with this family already
					if ((!(exists $comp_rel_history{"$ref_fam $ref_img $sub_fam u4"}) || ($comp_rel_history{"$ref_fam $ref_img $sub_fam u4"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_img $ref_fam u4"}) || ($comp_rel_history{"$sub_fam $sub_img $ref_fam u4"} ne $ref_fam))) {
						# create history entry
						$comp_rel_history{"$ref_fam $ref_img $sub_fam u4"} = $sub_fam;
						$comp_rel_history{"$sub_fam $sub_img $ref_fam u4"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $comp_relationships{"$ref_fam $sub_fam u4"}) {
							# now ref fam is upstream of sub fam as we are 
							# counting from right
							$comp_relationships{"$ref_fam $sub_fam u4"} = 1;
						}
						else{
							$comp_relationships{"$ref_fam $sub_fam u4"}++;
						}
						# add record for u4 relationship
						$comp_copies[$comp_copy_ctr][0]=$ref_fam; $comp_copies[$comp_copy_ctr][1]="u4";
						$comp_copies[$comp_copy_ctr][2]=$sub_fam; $comp_copies[$comp_copy_ctr][3]='C';
						$comp_copies[$comp_copy_ctr][4]=$ref_start; $comp_copies[$comp_copy_ctr][5]=$ref_end;
						$comp_copies[$comp_copy_ctr][6]=$sub_start; $comp_copies[$comp_copy_ctr++][7]=$sub_end;
					}
				}
			}
			# irrespective of strand
			# check if the ref image has this relation with this family already
			if ((!(exists $both_rel_history{"$ref_fam $ref_img $sub_fam u4"}) || ($both_rel_history{"$ref_fam $ref_img $sub_fam u4"} ne $sub_fam)) && (!(exists $both_rel_history{"$sub_fam $sub_img $ref_fam u4"}) || ($both_rel_history{"$sub_fam $sub_img $ref_fam u4"} ne $ref_fam))) {
				# create history entry
				$both_rel_history{"$ref_fam $ref_img $sub_fam u4"} = $sub_fam;
				$both_rel_history{"$sub_fam $sub_img $ref_fam u4"} = $ref_fam;

				
				# increment relationship count or create relationship entry
				if (!exists $both_relationships{"$sub_fam $ref_fam u4"}) {
					$both_relationships{"$sub_fam $ref_fam u4"} = 1;
				}
				else{
					$both_relationships{"$sub_fam $ref_fam u4"}++;
				}
				# add record for u4 relationship
				$both_copies[$both_copy_ctr][0]=$sub_fam; $both_copies[$both_copy_ctr][1]="u4";
				$both_copies[$both_copy_ctr][2]=$ref_fam; $both_copies[$both_copy_ctr][3]='B';
				$both_copies[$both_copy_ctr][4]=$sub_start; $both_copies[$both_copy_ctr][5]=$sub_end;
				$both_copies[$both_copy_ctr][6]=$ref_start; $both_copies[$both_copy_ctr++][7]=$ref_end;
			}
		}
		# Upstream: u5 (10000-15000 bases)
		elsif(($sub_end <= $ref_start-10000) && ($sub_end > $ref_start-15000)){
			if ($ref_strand eq $sub_strand){
				if($ref_strand eq '+'){# pos strand
					# check if the ref image has this relation with this family already
					if ((!(exists $pos_rel_history{"$ref_fam $ref_img $sub_fam u5"}) || ($pos_rel_history{"$ref_fam $ref_img $sub_fam u5"} ne $sub_fam)) && (!(exists $pos_rel_history{"$sub_fam $sub_img $ref_fam u5"}) || ($pos_rel_history{"$sub_fam $sub_img $ref_fam u5"} ne $ref_fam))) {
						# create history entry
						$pos_rel_history{"$ref_fam $ref_img $sub_fam u5"} = $sub_fam;
						$pos_rel_history{"$sub_fam $sub_img $ref_fam u5"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $pos_relationships{"$sub_fam $ref_fam u5"}) {
							$pos_relationships{"$sub_fam $ref_fam u5"} = 1;
						}
						else{
							$pos_relationships{"$sub_fam $ref_fam u5"}++;
						}
						# add record for u5 relationship
						$pos_copies[$pos_copy_ctr][0]=$sub_fam; $pos_copies[$pos_copy_ctr][1]="u5";
						$pos_copies[$pos_copy_ctr][2]=$ref_fam; $pos_copies[$pos_copy_ctr][3]='+';
						$pos_copies[$pos_copy_ctr][4]=$sub_start; $pos_copies[$pos_copy_ctr][5]=$sub_end;
						$pos_copies[$pos_copy_ctr][6]=$ref_start; $pos_copies[$pos_copy_ctr++][7]=$ref_end;
					}
				}
				elsif($ref_strand eq 'C'){# comp strand
					# check if the ref image has this relation with this family already
					if ((!(exists $comp_rel_history{"$ref_fam $ref_img $sub_fam u5"}) || ($comp_rel_history{"$ref_fam $ref_img $sub_fam u5"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_img $ref_fam u5"}) || ($comp_rel_history{"$sub_fam $sub_img $ref_fam u5"} ne $ref_fam))) {
						# create history entry
						$comp_rel_history{"$ref_fam $ref_img $sub_fam u5"} = $sub_fam;
						$comp_rel_history{"$sub_fam $sub_img $ref_fam u5"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $comp_relationships{"$ref_fam $sub_fam u5"}) {
							# now ref fam is upstream of sub fam as we are 
							# counting from right
							$comp_relationships{"$ref_fam $sub_fam u5"} = 1;
						}
						else{
							$comp_relationships{"$ref_fam $sub_fam u5"}++;
						}
						# add record for u5 relationship
						$comp_copies[$comp_copy_ctr][0]=$ref_fam; $comp_copies[$comp_copy_ctr][1]="u5";
						$comp_copies[$comp_copy_ctr][2]=$sub_fam; $comp_copies[$comp_copy_ctr][3]='B';
						$comp_copies[$comp_copy_ctr][4]=$ref_start; $comp_copies[$comp_copy_ctr][5]=$ref_end;
						$comp_copies[$comp_copy_ctr][6]=$sub_start; $comp_copies[$comp_copy_ctr++][7]=$sub_end;
					}
				}
			}
			# irrespective of strand
			# check if the ref image has this relation with this family already
			if ((!(exists $both_rel_history{"$ref_fam $ref_img $sub_fam u5"}) || ($both_rel_history{"$ref_fam $ref_img $sub_fam u5"} ne $sub_fam)) && (!(exists $both_rel_history{"$sub_fam $sub_img $ref_fam u5"}) || ($both_rel_history{"$sub_fam $sub_img $ref_fam u5"} ne $ref_fam))) {
				# create history entry
				$both_rel_history{"$ref_fam $ref_img $sub_fam u5"} = $sub_fam;
				$both_rel_history{"$sub_fam $sub_img $ref_fam u5"} = $ref_fam;
				
				# increment relationship count or create relationship entry
				if (!exists $both_relationships{"$sub_fam $ref_fam u5"}) {
					$both_relationships{"$sub_fam $ref_fam u5"} = 1;
				}
				else{
					$both_relationships{"$sub_fam $ref_fam u5"}++;
				}
				# add record for u5 relationship
				$both_copies[$both_copy_ctr][0]=$sub_fam; $both_copies[$both_copy_ctr][1]="u5";
				$both_copies[$both_copy_ctr][2]=$ref_fam; $both_copies[$both_copy_ctr][3]='B';
				$both_copies[$both_copy_ctr][4]=$sub_start; $both_copies[$both_copy_ctr][5]=$sub_end;
				$both_copies[$both_copy_ctr][6]=$ref_start; $both_copies[$both_copy_ctr++][7]=$ref_end;
			}
		}
		if($copy_file_flag){
			# temporary fix to reduce memory consumption when 
			# copies are not needed
			@pos_copies=();
			@comp_copies=();
			@both_copies=();# deallocating memory
			$pos_copy_ctr=$comp_copy_ctr=$both_copy_ctr=0;# resetting counters
		}
# 		print STDERR '.';
	}# end while

	$j=$i;
	
	# only look for relationships with images located after it
	# and starting within 15k bases after ref_end (enforced by condition above)
	# or anytime after ref_start (enforced by sorting the list on start pos)
	while(($j!=$#table) && ($table[$j+1][5] < $ref_end+15000)){
		$dns_ctr++;
		$j++;
		$sub_start=$table[$j][5]; $sub_end=$table[$j][6];
		$sub_strand=$table[$j][8]; $sub_fam=$table[$j][9];
		$sub_img=$table[$j][14];
		
		# cleaning up
		$sub_start=~ s/\s//g; $sub_end=~ s/\s//g;
		$sub_strand=~ s/\s//g; $sub_fam=~ s/\s//g;
		$sub_img=~ s/\s//g;
		
		
		# Note: since all relationship are exclusive, I have used elsif
		
		# In: Location of ref fam is entirely within sub fam (IN)
		# IN should be first bcos if sub start is near the ref start, it will
		# be listed right after the ref record in the list
		if(($sub_start == $ref_start) && ($sub_end >= $ref_end)){
			if ($ref_strand eq $sub_strand){
				if($ref_strand eq '+'){# pos strand
					# check if the ref image has this relation with this family already
					if ((!(exists $pos_rel_history{"$ref_fam $ref_img $sub_fam IN"}) || ($pos_rel_history{"$ref_fam $ref_img $sub_fam IN"} ne $sub_fam)) && (!(exists $pos_rel_history{"$sub_fam $sub_img $ref_fam IN"}) || ($pos_rel_history{"$sub_fam $sub_img $ref_fam IN"} ne $ref_fam))) {
						# create history entry
						$pos_rel_history{"$ref_fam $ref_img $sub_fam IN"} = $sub_fam;
						$pos_rel_history{"$sub_fam $sub_img $ref_fam IN"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $pos_relationships{"$ref_fam $sub_fam IN"}) {
							$pos_relationships{"$ref_fam $sub_fam IN"} = 1;
						}
						else{
							$pos_relationships{"$ref_fam $sub_fam IN"}++;
						}
						# add record for IN relationship
						$pos_copies[$pos_copy_ctr][0]=$ref_fam; $pos_copies[$pos_copy_ctr][1]="IN";
						$pos_copies[$pos_copy_ctr][2]=$sub_fam; $pos_copies[$pos_copy_ctr][3]='+';
						$pos_copies[$pos_copy_ctr][4]=$ref_start; $pos_copies[$pos_copy_ctr][5]=$ref_end;
						$pos_copies[$pos_copy_ctr][6]=$sub_start; $pos_copies[$pos_copy_ctr++][7]=$sub_end;
						
						# increment reciprocal relationship count or create relationship entry
						if (!exists $pos_relationships{"$sub_fam $ref_fam CONT"}) {
							$pos_relationships{"$sub_fam $ref_fam CONT"} = 1;
						}
						else{
							$pos_relationships{"$sub_fam $ref_fam CONT"}++;
						}
						# add record for CONT relationship
						$pos_copies[$pos_copy_ctr][0]=$sub_fam; $pos_copies[$pos_copy_ctr][1]="CONT";
						$pos_copies[$pos_copy_ctr][2]=$ref_fam; $pos_copies[$pos_copy_ctr][3]='+';
						$pos_copies[$pos_copy_ctr][4]=$sub_start; $pos_copies[$pos_copy_ctr][5]=$sub_end;
						$pos_copies[$pos_copy_ctr][6]=$ref_start; $pos_copies[$pos_copy_ctr++][7]=$ref_end;
					}
				}
				elsif($ref_strand eq 'C'){# comp strand
					# check if the ref image has this relation with this family already
					if ((!(exists $comp_rel_history{"$ref_fam $ref_img $sub_fam IN"}) || ($comp_rel_history{"$ref_fam $ref_img $sub_fam IN"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_img $ref_fam IN"}) || ($comp_rel_history{"$sub_fam $sub_img $ref_fam IN"} ne $ref_fam))) {
						# create history entry
						$comp_rel_history{"$ref_fam $ref_img $sub_fam IN"} = $sub_fam;
						$comp_rel_history{"$sub_fam $sub_img $ref_fam IN"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $comp_relationships{"$ref_fam $sub_fam IN"}) {
							$comp_relationships{"$ref_fam $sub_fam IN"} = 1;
						}
						else{
							$comp_relationships{"$ref_fam $sub_fam IN"}++;
						}
						# add record for IN relationship
						$comp_copies[$comp_copy_ctr][0]=$ref_fam; $comp_copies[$comp_copy_ctr][1]="IN";
						$comp_copies[$comp_copy_ctr][2]=$sub_fam; $comp_copies[$comp_copy_ctr][3]='C';
						$comp_copies[$comp_copy_ctr][4]=$ref_start; $comp_copies[$comp_copy_ctr][5]=$ref_end;
						$comp_copies[$comp_copy_ctr][6]=$sub_start; $comp_copies[$comp_copy_ctr++][7]=$sub_end;
						
						# increment reciprocal relationship count or create relationship entry
						if (!exists $comp_relationships{"$sub_fam $ref_fam CONT"}) {
							$comp_relationships{"$sub_fam $ref_fam CONT"} = 1;
						}
						else{
							$comp_relationships{"$sub_fam $ref_fam CONT"}++;
						}
						# add record for CONT relationship
						$comp_copies[$comp_copy_ctr][0]=$sub_fam; $comp_copies[$comp_copy_ctr][1]="CONT";
						$comp_copies[$comp_copy_ctr][2]=$ref_fam; $comp_copies[$comp_copy_ctr][3]='C';
						$comp_copies[$comp_copy_ctr][4]=$sub_start; $comp_copies[$comp_copy_ctr][5]=$sub_end;
						$comp_copies[$comp_copy_ctr][6]=$ref_start; $comp_copies[$comp_copy_ctr++][7]=$ref_end;
					}
				}
			}
			# irrespective of strand
			# check if the ref image has this relation with this family already
			if ((!(exists $both_rel_history{"$ref_fam $ref_img $sub_fam IN"}) || ($both_rel_history{"$ref_fam $ref_img $sub_fam IN"} ne $sub_fam)) && (!(exists $both_rel_history{"$sub_fam $sub_img $ref_fam IN"}) || ($both_rel_history{"$sub_fam $sub_img $ref_fam IN"} ne $ref_fam))) {
				# create history entry
				$both_rel_history{"$ref_fam $ref_img $sub_fam IN"} = $sub_fam;
				$both_rel_history{"$sub_fam $sub_img $ref_fam IN"} = $ref_fam;
				
				# increment relationship count or create relationship entry
				if (!exists $both_relationships{"$ref_fam $sub_fam IN"}) {
					$both_relationships{"$ref_fam $sub_fam IN"} = 1;
				}
				else{
					$both_relationships{"$ref_fam $sub_fam IN"}++;
				}
				# add record for IN relationship
				$both_copies[$both_copy_ctr][0]=$ref_fam; $both_copies[$both_copy_ctr][1]="IN";
				$both_copies[$both_copy_ctr][2]=$sub_fam; $both_copies[$both_copy_ctr][3]='B';
				$both_copies[$both_copy_ctr][4]=$ref_start; $both_copies[$both_copy_ctr][5]=$ref_end;
				$both_copies[$both_copy_ctr][6]=$sub_start; $both_copies[$both_copy_ctr++][7]=$sub_end;

				# increment reciprocal relationship count or create relationship entry
				if (!exists $both_relationships{"$sub_fam $ref_fam CONT"}) {
					$both_relationships{"$sub_fam $ref_fam CONT"} = 1;
				}
				else{
					$both_relationships{"$sub_fam $ref_fam CONT"}++;
				}
				# add record for CONT relationship
				$both_copies[$both_copy_ctr][0]=$sub_fam; $both_copies[$both_copy_ctr][1]="CONT";
				$both_copies[$both_copy_ctr][2]=$ref_fam; $both_copies[$both_copy_ctr][3]='B';
				$both_copies[$both_copy_ctr][4]=$sub_start; $both_copies[$both_copy_ctr][5]=$sub_end;
				$both_copies[$both_copy_ctr][6]=$ref_start; $both_copies[$both_copy_ctr++][7]=$ref_end;

			}
		}# IN end
		
		# Overlap: If overlap is more than 10% of length of either family (Ovlap)
		# now if subject fam ends within the reference fam
		elsif (($sub_start > $ref_start) &&  ($sub_start < $ref_end)) {
			my ($ovlap, $ref_ovlap, $sub_ovlap);
			$ovlap = $ref_end - $sub_start;
			
			$ref_ovlap = ($ovlap / ($ref_end - $ref_start)) * 100;
			$sub_ovlap = ($ovlap / ($sub_end - $sub_start)) * 100;
			# Overlap :10% to 30% (Ovlap-10to30)
			if ((($ref_ovlap > 10.00) && ($ref_ovlap <= 30.00)) || 
			(($sub_ovlap > 10.00) && ($sub_ovlap <= 30.00))) {
				if ($ref_strand eq $sub_strand){
					if($ref_strand eq '+'){# pos strand
						# check if the ref image has this relation with this family already
						if ((!(exists $pos_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-10to30"}) || ($pos_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-10to30"} ne $sub_fam)) && (!(exists $pos_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-10to30"}) || ($pos_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-10to30"} ne $ref_fam))) {
							# create history entry
							$pos_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-10to30"} = $sub_fam;
							$pos_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-10to30"} = $ref_fam;
							
							# increment relationship count or create relationship entry
							if (!exists $pos_relationships{"$sub_fam $ref_fam Ovlap-10to30"}) {
								$pos_relationships{"$sub_fam $ref_fam Ovlap-10to30"} = 1;
							}
							else{
								$pos_relationships{"$sub_fam $ref_fam Ovlap-10to30"}++;
							}
							# add record for Ovlap-10to30 relationship
							$pos_copies[$pos_copy_ctr][0]=$sub_fam; $pos_copies[$pos_copy_ctr][1]="Ovlap-10to30";
							$pos_copies[$pos_copy_ctr][2]=$ref_fam; $pos_copies[$pos_copy_ctr][3]='+';
							$pos_copies[$pos_copy_ctr][4]=$sub_start; $pos_copies[$pos_copy_ctr][5]=$sub_end;
							$pos_copies[$pos_copy_ctr][6]=$ref_start; $pos_copies[$pos_copy_ctr++][7]=$ref_end;
						}
					}
					elsif($ref_strand eq 'C'){# comp strand
						# check if the ref image has this relation with this family already
						if ((!(exists $comp_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-10to30"}) || ($comp_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-10to30"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-10to30"}) || ($comp_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-10to30"} ne $ref_fam))) {
							# create history entry
							$comp_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-10to30"} = $sub_fam;
							$comp_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-10to30"} = $ref_fam;
							
							# increment relationship count or create relationship entry
							if (!exists $comp_relationships{"$sub_fam $ref_fam Ovlap-10to30"}) {
								$comp_relationships{"$sub_fam $ref_fam Ovlap-10to30"} = 1;
							}
							else{
								$comp_relationships{"$sub_fam $ref_fam Ovlap-10to30"}++;
							}
							# add record for Ovlap-10to30 relationship
							$comp_copies[$comp_copy_ctr][0]=$sub_fam; $comp_copies[$comp_copy_ctr][1]="Ovlap-10to30";
							$comp_copies[$comp_copy_ctr][2]=$ref_fam; $comp_copies[$comp_copy_ctr][3]='C';
							$comp_copies[$comp_copy_ctr][4]=$sub_start; $comp_copies[$comp_copy_ctr][5]=$sub_end;
							$comp_copies[$comp_copy_ctr][6]=$ref_start; $comp_copies[$comp_copy_ctr++][7]=$ref_end;
						}
					}
				}
				# irrespective of strand
				if ((!(exists $both_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-10to30"}) || ($both_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-10to30"} ne $sub_fam)) && (!(exists $both_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-10to30"}) || ($both_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-10to30"} ne $ref_fam))) {
					# create history entry
					$both_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-10to30"} = $sub_fam;
					$both_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-10to30"} = $ref_fam;

					
					# increment relationship count or create relationship entry
					if (!exists $both_relationships{"$sub_fam $ref_fam Ovlap-10to30"}) {
						$both_relationships{"$sub_fam $ref_fam Ovlap-10to30"} = 1;
					}
					else{
						$both_relationships{"$sub_fam $ref_fam Ovlap-10to30"}++;
					}
					# add record for Ovlap-10to30 relationship
					$both_copies[$both_copy_ctr][0]=$sub_fam; $both_copies[$both_copy_ctr][1]="Ovlap-10to30";
					$both_copies[$both_copy_ctr][2]=$ref_fam; $both_copies[$both_copy_ctr][3]='B';
					$both_copies[$both_copy_ctr][4]=$sub_start; $both_copies[$both_copy_ctr][5]=$sub_end;
					$both_copies[$both_copy_ctr][6]=$ref_start; $both_copies[$both_copy_ctr++][7]=$ref_end;

				}
			}
			# Overlap :30% to 70% (Ovlap-30to70)
			elsif ((($ref_ovlap > 30.00) && ($ref_ovlap <= 70.00)) || 
			(($sub_ovlap > 30.00) && ($sub_ovlap <= 70.00))) {
				if ($ref_strand eq $sub_strand){
					if($ref_strand eq '+'){# pos strand
						# check if the ref image has this relation with this family already
						if ((!(exists $pos_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-30to70"}) || ($pos_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-30to70"} ne $sub_fam)) && (!(exists $pos_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-30to70"}) || ($pos_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-30to70"} ne $ref_fam))) {
							# create history entry
							$pos_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-30to70"} = $sub_fam;
							$pos_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-30to70"} = $ref_fam;
							
							# increment relationship count or create relationship entry
							if (!exists $pos_relationships{"$sub_fam $ref_fam Ovlap-30to70"}) {
								$pos_relationships{"$sub_fam $ref_fam Ovlap-30to70"} = 1;
							}
							else{
								$pos_relationships{"$sub_fam $ref_fam Ovlap-30to70"}++;
							}
							# add record for Ovlap-30to70 relationship
							$pos_copies[$pos_copy_ctr][0]=$sub_fam; $pos_copies[$pos_copy_ctr][1]="Ovlap-30to70";
							$pos_copies[$pos_copy_ctr][2]=$ref_fam; $pos_copies[$pos_copy_ctr][3]='+';
							$pos_copies[$pos_copy_ctr][4]=$sub_start; $pos_copies[$pos_copy_ctr][5]=$sub_end;
							$pos_copies[$pos_copy_ctr][6]=$ref_start; $pos_copies[$pos_copy_ctr++][7]=$ref_end;

						}
					}
					elsif($ref_strand eq 'C'){# comp strand
						# check if the ref image has this relation with this family already
						if ((!(exists $comp_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-30to70"}) || ($comp_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-30to70"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-30to70"}) || ($comp_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-30to70"} ne $ref_fam))) {
							# create history entry
							$comp_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-30to70"} = $sub_fam;
							$comp_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-30to70"} = $ref_fam;
							
							# increment relationship count or create relationship entry
							if (!exists $comp_relationships{"$sub_fam $ref_fam Ovlap-30to70"}) {
								$comp_relationships{"$sub_fam $ref_fam Ovlap-30to70"} = 1;
							}
							else{
								$comp_relationships{"$sub_fam $ref_fam Ovlap-30to70"}++;
							}
							# add record for Ovlap-30to70 relationship
							$comp_copies[$comp_copy_ctr][0]=$sub_fam; $comp_copies[$comp_copy_ctr][1]="Ovlap-30to70";
							$comp_copies[$comp_copy_ctr][2]=$ref_fam; $comp_copies[$comp_copy_ctr][3]='C';
							$comp_copies[$comp_copy_ctr][4]=$sub_start; $comp_copies[$comp_copy_ctr][5]=$sub_end;
							$comp_copies[$comp_copy_ctr][6]=$ref_start; $comp_copies[$comp_copy_ctr++][7]=$ref_end;
						}
					}
				}
				# irrespective of strand
				# check if the ref image has this relation with this family already
				if ((!(exists $both_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-30to70"}) || ($both_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-30to70"} ne $sub_fam)) && (!(exists $both_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-30to70"}) || ($both_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-30to70"} ne $ref_fam))) {
					# create history entry
					$both_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-30to70"} = $sub_fam;
					$both_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-30to70"} = $ref_fam;
					
					# increment relationship count or create relationship entry
					if (!exists $both_relationships{"$sub_fam $ref_fam Ovlap-30to70"}) {
						$both_relationships{"$sub_fam $ref_fam Ovlap-30to70"} = 1;
					}
					else{
						$both_relationships{"$sub_fam $ref_fam Ovlap-30to70"}++;
					}
					# add record for Ovlap-30to70 relationship
					$both_copies[$both_copy_ctr][0]=$sub_fam; $both_copies[$both_copy_ctr][1]="Ovlap-30to70";
					$both_copies[$both_copy_ctr][2]=$ref_fam; $both_copies[$both_copy_ctr][3]='B';
					$both_copies[$both_copy_ctr][4]=$sub_start; $both_copies[$both_copy_ctr][5]=$sub_end;
					$both_copies[$both_copy_ctr][6]=$ref_start; $both_copies[$both_copy_ctr++][7]=$ref_end;
				}
			}
			
			# Overlap : >70% (Ovlap-70plus)
			elsif (($ref_ovlap > 70.00) || ($sub_ovlap > 70.00)) {
				if ($ref_strand eq $sub_strand){
					if($ref_strand eq '+'){# pos strand
						# check if the ref image has this relation with this family already
						if ((!(exists $pos_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-70plus"}) || ($pos_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-70plus"} ne $sub_fam)) && (!(exists $pos_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-70plus"}) || ($pos_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-70plus"} ne $ref_fam))) {
							# create history entry
							$pos_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-70plus"} = $sub_fam;
							$pos_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-70plus"} = $ref_fam;
							
							# increment relationship count or create relationship entry
							if (!exists $pos_relationships{"$sub_fam $ref_fam Ovlap-70plus"}) {
								$pos_relationships{"$sub_fam $ref_fam Ovlap-70plus"} = 1;
							}
							else{
								$pos_relationships{"$sub_fam $ref_fam Ovlap-70plus"}++;
							}
							# add record for Ovlap-70plus relationship
							$pos_copies[$pos_copy_ctr][0]=$sub_fam; $pos_copies[$pos_copy_ctr][1]="Ovlap-70plus";
							$pos_copies[$pos_copy_ctr][2]=$ref_fam; $pos_copies[$pos_copy_ctr][3]='+';
							$pos_copies[$pos_copy_ctr][4]=$sub_start; $pos_copies[$pos_copy_ctr][5]=$sub_end;
							$pos_copies[$pos_copy_ctr][6]=$ref_start; $pos_copies[$pos_copy_ctr++][7]=$ref_end;
						}
					}
					elsif($ref_strand eq 'C'){# comp strand
						# check if the ref image has this relation with this family already
						if ((!(exists $comp_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-70plus"}) || ($comp_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-70plus"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-70plus"}) || ($comp_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-70plus"} ne $ref_fam))) {
							# create history entry
							$comp_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-70plus"} = $sub_fam;
							$comp_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-70plus"} = $ref_fam;
							
							# increment relationship count or create relationship entry
							if (!exists $comp_relationships{"$sub_fam $ref_fam Ovlap-70plus"}) {
								$comp_relationships{"$sub_fam $ref_fam Ovlap-70plus"} = 1;
							}
							else{
								$comp_relationships{"$sub_fam $ref_fam Ovlap-70plus"}++;
							}
							# add record for Ovlap-70plus relationship
							$comp_copies[$comp_copy_ctr][0]=$sub_fam; $comp_copies[$comp_copy_ctr][1]="Ovlap-70plus";
							$comp_copies[$comp_copy_ctr][2]=$ref_fam; $comp_copies[$comp_copy_ctr][3]='C';
							$comp_copies[$comp_copy_ctr][4]=$sub_start; $comp_copies[$comp_copy_ctr][5]=$sub_end;
							$comp_copies[$comp_copy_ctr][6]=$ref_start; $comp_copies[$comp_copy_ctr++][7]=$ref_end;
						}
					}
				}
				# irrespective of strand
				# check if the ref image has this relation with this family already
				if ((!(exists $both_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-70plus"}) || ($both_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-70plus"} ne $sub_fam)) && (!(exists $both_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-70plus"}) || ($both_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-70plus"} ne $ref_fam))) {
					# create history entry
					$both_rel_history{"$ref_fam $ref_img $sub_fam Ovlap-70plus"} = $sub_fam;
					$both_rel_history{"$sub_fam $sub_img $ref_fam Ovlap-70plus"} = $ref_fam;
					
					# increment relationship count or create relationship entry
					if (!exists $both_relationships{"$sub_fam $ref_fam Ovlap-70plus"}) {
						$both_relationships{"$sub_fam $ref_fam Ovlap-70plus"} = 1;
					}
					else{
						$both_relationships{"$sub_fam $ref_fam Ovlap-70plus"}++;
					}
					# add record for Ovlap-70plus relationship
					$both_copies[$both_copy_ctr][0]=$sub_fam; $both_copies[$both_copy_ctr][1]="Ovlap-70plus";
					$both_copies[$both_copy_ctr][2]=$ref_fam; $both_copies[$both_copy_ctr][3]='B';
					$both_copies[$both_copy_ctr][4]=$sub_start; $both_copies[$both_copy_ctr][5]=$sub_end;
					$both_copies[$both_copy_ctr][6]=$ref_start; $both_copies[$both_copy_ctr++][7]=$ref_end;
				}
			}
		}# overlap end
		
		# Downstream: d1 (0-500 bases)
		elsif(($sub_start >= $ref_end) && ($sub_start < $ref_end+500)){
			if ($ref_strand eq $sub_strand){
				if($ref_strand eq '+'){# pos strand
					# check if the ref image has this relation with this family already
					if ((!(exists $pos_rel_history{"$ref_fam $ref_img $sub_fam d1"}) || ($pos_rel_history{"$ref_fam $ref_img $sub_fam d1"} ne $sub_fam)) && (!(exists $pos_rel_history{"$sub_fam $sub_img $ref_fam d1"}) || ($pos_rel_history{"$sub_fam $sub_img $ref_fam d1"} ne $ref_fam))) {
						# create history entry
						$pos_rel_history{"$ref_fam $ref_img $sub_fam d1"} = $sub_fam;
						$pos_rel_history{"$sub_fam $sub_img $ref_fam d1"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $pos_relationships{"$sub_fam $ref_fam d1"}) {
							$pos_relationships{"$sub_fam $ref_fam d1"} = 1;
						}
						else{
							$pos_relationships{"$sub_fam $ref_fam d1"}++;
						}
						# add record for d1 relationship
						$pos_copies[$pos_copy_ctr][0]=$sub_fam; $pos_copies[$pos_copy_ctr][1]="d1";
						$pos_copies[$pos_copy_ctr][2]=$ref_fam; $pos_copies[$pos_copy_ctr][3]='+';
						$pos_copies[$pos_copy_ctr][4]=$sub_start; $pos_copies[$pos_copy_ctr][5]=$sub_end;
						$pos_copies[$pos_copy_ctr][6]=$ref_start; $pos_copies[$pos_copy_ctr++][7]=$ref_end;
					}
				}
				elsif($ref_strand eq 'C'){# comp strand
					# check if the ref image has this relation with this family already
					if ((!(exists $comp_rel_history{"$ref_fam $ref_img $sub_fam d1"}) || ($comp_rel_history{"$ref_fam $ref_img $sub_fam d1"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_img $ref_fam d1"}) || ($comp_rel_history{"$sub_fam $sub_img $ref_fam d1"} ne $ref_fam))) {
						# create history entry
						$comp_rel_history{"$ref_fam $ref_img $sub_fam d1"} = $sub_fam;
						$comp_rel_history{"$sub_fam $sub_img $ref_fam d1"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $comp_relationships{"$ref_fam $sub_fam d1"}) {
							# now ref fam is upstream of sub fam as we are 
							# counting from right
							$comp_relationships{"$ref_fam $sub_fam d1"} = 1;
						}
						else{
							$comp_relationships{"$ref_fam $sub_fam d1"}++;
						}
						# add record for d1 relationship
						$comp_copies[$comp_copy_ctr][0]=$ref_fam; $comp_copies[$comp_copy_ctr][1]="d1";
						$comp_copies[$comp_copy_ctr][2]=$sub_fam; $comp_copies[$comp_copy_ctr][3]='C';
						$comp_copies[$comp_copy_ctr][4]=$ref_start; $comp_copies[$comp_copy_ctr][5]=$ref_end;
						$comp_copies[$comp_copy_ctr][6]=$sub_start; $comp_copies[$comp_copy_ctr++][7]=$sub_end;
					}
				}
			}
			# irrespective of strand
			# check if the ref image has this relation with this family already
			if ((!(exists $both_rel_history{"$ref_fam $ref_img $sub_fam d1"}) || ($both_rel_history{"$ref_fam $ref_img $sub_fam d1"} ne $sub_fam)) && (!(exists $both_rel_history{"$sub_fam $sub_img $ref_fam d1"}) || ($both_rel_history{"$sub_fam $sub_img $ref_fam d1"} ne $ref_fam))) {
				# create history entry
				$both_rel_history{"$ref_fam $ref_img $sub_fam d1"} = $sub_fam;
				$both_rel_history{"$sub_fam $sub_img $ref_fam d1"} = $ref_fam;
				
				# increment relationship count or create relationship entry
				if (!exists $both_relationships{"$sub_fam $ref_fam d1"}) {
					$both_relationships{"$sub_fam $ref_fam d1"} = 1;
				}
				else{
					$both_relationships{"$sub_fam $ref_fam d1"}++;
				}
				# add record for d1 relationship
				$both_copies[$both_copy_ctr][0]=$sub_fam; $both_copies[$both_copy_ctr][1]="d1";
				$both_copies[$both_copy_ctr][2]=$ref_fam; $both_copies[$both_copy_ctr][3]='B';
				$both_copies[$both_copy_ctr][4]=$sub_start; $both_copies[$both_copy_ctr][5]=$sub_end;
				$both_copies[$both_copy_ctr][6]=$ref_start; $both_copies[$both_copy_ctr++][7]=$ref_end;
			}
		}

		# Downstream: d2 (500-1000 bases)
		elsif(($sub_start >= $ref_end+500) && ($sub_start < $ref_end+1000)){
			if ($ref_strand eq $sub_strand){
				if($ref_strand eq '+'){# pos strand
					# check if the ref image has this relation with this family already
					if ((!(exists $pos_rel_history{"$ref_fam $ref_img $sub_fam d2"}) || ($pos_rel_history{"$ref_fam $ref_img $sub_fam d2"} ne $sub_fam)) && (!(exists $pos_rel_history{"$sub_fam $sub_img $ref_fam d2"}) || ($pos_rel_history{"$sub_fam $sub_img $ref_fam d2"} ne $ref_fam))) {
						# create history entry
						$pos_rel_history{"$ref_fam $ref_img $sub_fam d2"} = $sub_fam;
						$pos_rel_history{"$sub_fam $sub_img $ref_fam d2"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $pos_relationships{"$sub_fam $ref_fam d2"}) {
							$pos_relationships{"$sub_fam $ref_fam d2"} = 1;
						}
						else{
							$pos_relationships{"$sub_fam $ref_fam d2"}++;
						}
						# add record for d2 relationship
						$pos_copies[$pos_copy_ctr][0]=$sub_fam; $pos_copies[$pos_copy_ctr][1]="d2";
						$pos_copies[$pos_copy_ctr][2]=$ref_fam; $pos_copies[$pos_copy_ctr][3]='+';
						$pos_copies[$pos_copy_ctr][4]=$sub_start; $pos_copies[$pos_copy_ctr][5]=$sub_end;
						$pos_copies[$pos_copy_ctr][6]=$ref_start; $pos_copies[$pos_copy_ctr++][7]=$ref_end;
					}
				}
				elsif($ref_strand eq 'C'){# comp strand
					# check if the ref image has this relation with this family already
					if ((!(exists $comp_rel_history{"$ref_fam $ref_img $sub_fam d2"}) || ($comp_rel_history{"$ref_fam $ref_img $sub_fam d2"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_img $ref_fam d2"}) || ($comp_rel_history{"$sub_fam $sub_img $ref_fam d2"} ne $ref_fam))) {
						# create history entry
						$comp_rel_history{"$ref_fam $ref_img $sub_fam d2"} = $sub_fam;
						$comp_rel_history{"$sub_fam $sub_img $ref_fam d2"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $comp_relationships{"$ref_fam $sub_fam d2"}) {
							# now ref fam is upstream of sub fam as we are 
							# counting from right
							$comp_relationships{"$ref_fam $sub_fam d2"} = 1;
						}
						else{
							$comp_relationships{"$ref_fam $sub_fam d2"}++;
						}
						# add record for d2 relationship
						$comp_copies[$comp_copy_ctr][0]=$ref_fam; $comp_copies[$comp_copy_ctr][1]="d2";
						$comp_copies[$comp_copy_ctr][2]=$sub_fam; $comp_copies[$comp_copy_ctr][3]='C';
						$comp_copies[$comp_copy_ctr][4]=$ref_start; $comp_copies[$comp_copy_ctr][5]=$ref_end;
						$comp_copies[$comp_copy_ctr][6]=$sub_start; $comp_copies[$comp_copy_ctr++][7]=$sub_end;
					}
				}
			}
			# irrespective of strand
			# check if the ref image has this relation with this family already
			if ((!(exists $both_rel_history{"$ref_fam $ref_img $sub_fam d2"}) || ($both_rel_history{"$ref_fam $ref_img $sub_fam d2"} ne $sub_fam)) && (!(exists $both_rel_history{"$sub_fam $sub_img $ref_fam d2"}) || ($both_rel_history{"$sub_fam $sub_img $ref_fam d2"} ne $ref_fam))) {
				# create history entry
				$both_rel_history{"$ref_fam $ref_img $sub_fam d2"} = $sub_fam;
				$both_rel_history{"$sub_fam $sub_img $ref_fam d2"} = $ref_fam;
				
				# increment relationship count or create relationship entry
				if (!exists $both_relationships{"$sub_fam $ref_fam d2"}) {
					$both_relationships{"$sub_fam $ref_fam d2"} = 1;
				}
				else{
					$both_relationships{"$sub_fam $ref_fam d2"}++;
				}
				# add record for d2 relationship
				$both_copies[$both_copy_ctr][0]=$sub_fam; $both_copies[$both_copy_ctr][1]="d2";
				$both_copies[$both_copy_ctr][2]=$ref_fam; $both_copies[$both_copy_ctr][3]='B';
				$both_copies[$both_copy_ctr][4]=$sub_start; $both_copies[$both_copy_ctr][5]=$sub_end;
				$both_copies[$both_copy_ctr][6]=$ref_start; $both_copies[$both_copy_ctr++][7]=$ref_end;
			}
		}
		# Downstream: d3 (1000-5000 bases)
		elsif(($sub_start >= $ref_end+1000) && ($sub_start < $ref_end+5000)){
			if ($ref_strand eq $sub_strand){
				if($ref_strand eq '+'){# pos strand
					# check if the ref image has this relation with this family already
					if ((!(exists $pos_rel_history{"$ref_fam $ref_img $sub_fam d3"}) || ($pos_rel_history{"$ref_fam $ref_img $sub_fam d3"} ne $sub_fam)) && (!(exists $pos_rel_history{"$sub_fam $sub_img $ref_fam d3"}) || ($pos_rel_history{"$sub_fam $sub_img $ref_fam d3"} ne $ref_fam))) {
						# create history entry
						$pos_rel_history{"$ref_fam $ref_img $sub_fam d3"} = $sub_fam;
						$pos_rel_history{"$sub_fam $sub_img $ref_fam d3"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $pos_relationships{"$sub_fam $ref_fam d3"}) {
							$pos_relationships{"$sub_fam $ref_fam d3"} = 1;
						}
						else{
							$pos_relationships{"$sub_fam $ref_fam d3"}++;
						}
						# add record for d3 relationship
						$pos_copies[$pos_copy_ctr][0]=$sub_fam; $pos_copies[$pos_copy_ctr][1]="d3";
						$pos_copies[$pos_copy_ctr][2]=$ref_fam; $pos_copies[$pos_copy_ctr][3]='+';
						$pos_copies[$pos_copy_ctr][4]=$sub_start; $pos_copies[$pos_copy_ctr][5]=$sub_end;
						$pos_copies[$pos_copy_ctr][6]=$ref_start; $pos_copies[$pos_copy_ctr++][7]=$ref_end;
					}
				}
				elsif($ref_strand eq 'C'){# comp strand
					# check if the ref image has this relation with this family already
					if ((!(exists $comp_rel_history{"$ref_fam $ref_img $sub_fam d3"}) || ($comp_rel_history{"$ref_fam $ref_img $sub_fam d3"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_img $ref_fam d3"}) || ($comp_rel_history{"$sub_fam $sub_img $ref_fam d3"} ne $ref_fam))) {
						# create history entry
						$comp_rel_history{"$ref_fam $ref_img $sub_fam d3"} = $sub_fam;
						$comp_rel_history{"$sub_fam $sub_img $ref_fam d3"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $comp_relationships{"$ref_fam $sub_fam d3"}) {
							# now ref fam is upstream of sub fam as we are 
							# counting from right
							$comp_relationships{"$ref_fam $sub_fam d3"} = 1;
						}
						else{
							$comp_relationships{"$ref_fam $sub_fam d3"}++;
						}
						# add record for d3 relationship
						$comp_copies[$comp_copy_ctr][0]=$ref_fam; $comp_copies[$comp_copy_ctr][1]="d3";
						$comp_copies[$comp_copy_ctr][2]=$sub_fam; $comp_copies[$comp_copy_ctr][3]='C';
						$comp_copies[$comp_copy_ctr][4]=$ref_start; $comp_copies[$comp_copy_ctr][5]=$ref_end;
						$comp_copies[$comp_copy_ctr][6]=$sub_start; $comp_copies[$comp_copy_ctr++][7]=$sub_end;
					}
				}
			}
			# irrespective of strand
			# check if the ref image has this relation with this family already
			if ((!(exists $both_rel_history{"$ref_fam $ref_img $sub_fam d3"}) || ($both_rel_history{"$ref_fam $ref_img $sub_fam d3"} ne $sub_fam)) && (!(exists $both_rel_history{"$sub_fam $sub_img $ref_fam d3"}) || ($both_rel_history{"$sub_fam $sub_img $ref_fam d3"} ne $ref_fam))) {
				# create history entry
				$both_rel_history{"$ref_fam $ref_img $sub_fam d3"} = $sub_fam;
				$both_rel_history{"$sub_fam $sub_img $ref_fam d3"} = $ref_fam;

				
				# increment relationship count or create relationship entry
				if (!exists $both_relationships{"$sub_fam $ref_fam d3"}) {
					$both_relationships{"$sub_fam $ref_fam d3"} = 1;
				}
				else{
					$both_relationships{"$sub_fam $ref_fam d3"}++;
				}
				# add record for d3 relationship
				$both_copies[$both_copy_ctr][0]=$sub_fam; $both_copies[$both_copy_ctr][1]="d3";
				$both_copies[$both_copy_ctr][2]=$ref_fam; $both_copies[$both_copy_ctr][3]='B';
				$both_copies[$both_copy_ctr][4]=$sub_start; $both_copies[$both_copy_ctr][5]=$sub_end;
				$both_copies[$both_copy_ctr][6]=$ref_start; $both_copies[$both_copy_ctr++][7]=$ref_end;
			}
		}
		# Downstream: d4 (5000-10000 bases)
		elsif(($sub_start >= $ref_end+5000) && ($sub_start < $ref_end+10000)){
			if ($ref_strand eq $sub_strand){
				if($ref_strand eq '+'){# pos strand
					# check if the ref image has this relation with this family already
					if ((!(exists $pos_rel_history{"$ref_fam $ref_img $sub_fam d4"}) || ($pos_rel_history{"$ref_fam $ref_img $sub_fam d4"} ne $sub_fam)) && (!(exists $pos_rel_history{"$sub_fam $sub_img $ref_fam d4"}) || ($pos_rel_history{"$sub_fam $sub_img $ref_fam d4"} ne $ref_fam))) {
						# create history entry
						$pos_rel_history{"$ref_fam $ref_img $sub_fam d4"} = $sub_fam;
						$pos_rel_history{"$sub_fam $sub_img $ref_fam d4"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $pos_relationships{"$sub_fam $ref_fam d4"}) {
							$pos_relationships{"$sub_fam $ref_fam d4"} = 1;
						}
						else{
							$pos_relationships{"$sub_fam $ref_fam d4"}++;
						}
						# add record for d4 relationship
						$pos_copies[$pos_copy_ctr][0]=$ref_fam; $pos_copies[$pos_copy_ctr][1]="d4";
						$pos_copies[$pos_copy_ctr][2]=$sub_fam; $pos_copies[$pos_copy_ctr][3]='+';
						$pos_copies[$pos_copy_ctr][4]=$ref_start; $pos_copies[$pos_copy_ctr][5]=$ref_end;
						$pos_copies[$pos_copy_ctr][6]=$sub_start; $pos_copies[$pos_copy_ctr++][7]=$sub_end;
					}
				}
				elsif($ref_strand eq 'C'){# comp strand
					# check if the ref image has this relation with this family already
					if ((!(exists $comp_rel_history{"$ref_fam $ref_img $sub_fam d4"}) || ($comp_rel_history{"$ref_fam $ref_img $sub_fam d4"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_img $ref_fam d4"}) || ($comp_rel_history{"$sub_fam $sub_img $ref_fam d4"} ne $ref_fam))) {
						# create history entry
						$comp_rel_history{"$ref_fam $ref_img $sub_fam d4"} = $sub_fam;
						$comp_rel_history{"$sub_fam $sub_img $ref_fam d4"} = $ref_fam;
						
# 						# debugging
# 						if($ref_fam eq "R=759" && $sub_fam eq "R=759"){
# 							print ERRFILE "\nREF image data:\n";
# 							foreach(0..14){ print ERRFILE $table[$i][$_],' ';}
# 							print ERRFILE "\n";
# 							print ERRFILE "SUB image data:\n";
# 							foreach(0..14){ print ERRFILE $table[$j][$_],' ';}
# 							print ERRFILE "\n\n";
# 						}
						
						# increment relationship count or create relationship entry
						if (!exists $comp_relationships{"$ref_fam $sub_fam d4"}) {
							# now ref fam is upstream of sub fam as we are 
							# counting from right
							$comp_relationships{"$ref_fam $sub_fam d4"} = 1;
						}
						else{
							$comp_relationships{"$ref_fam $sub_fam d4"}++;
						}
						# add record for d4 relationship
						$comp_copies[$comp_copy_ctr][0]=$ref_fam; $comp_copies[$comp_copy_ctr][1]="d4";
						$comp_copies[$comp_copy_ctr][2]=$sub_fam; $comp_copies[$comp_copy_ctr][3]='C';
						$comp_copies[$comp_copy_ctr][4]=$ref_start; $comp_copies[$comp_copy_ctr][5]=$ref_end;
						$comp_copies[$comp_copy_ctr][6]=$sub_start; $comp_copies[$comp_copy_ctr++][7]=$sub_end;
					}
				}
			}
			# irrespective of strand
			# check if the ref image has this relation with this family already
			if ((!(exists $both_rel_history{"$ref_fam $ref_img $sub_fam d4"}) || ($both_rel_history{"$ref_fam $ref_img $sub_fam d4"} ne $sub_fam)) && (!(exists $both_rel_history{"$sub_fam $sub_img $ref_fam d4"}) || ($both_rel_history{"$sub_fam $sub_img $ref_fam d4"} ne $ref_fam))) {
				# create history entry
				$both_rel_history{"$ref_fam $ref_img $sub_fam d4"} = $sub_fam;
				$both_rel_history{"$sub_fam $sub_img $ref_fam d4"} = $ref_fam;
				
				# increment relationship count or create relationship entry
				if (!exists $both_relationships{"$sub_fam $ref_fam d4"}) {
					$both_relationships{"$sub_fam $ref_fam d4"} = 1;
				}
				else{
					$both_relationships{"$sub_fam $ref_fam d4"}++;
				}
				# add record for d4 relationship
				$both_copies[$both_copy_ctr][0]=$sub_fam; $both_copies[$both_copy_ctr][1]="d4";
				$both_copies[$both_copy_ctr][2]=$ref_fam; $both_copies[$both_copy_ctr][3]='B';
				$both_copies[$both_copy_ctr][4]=$sub_start; $both_copies[$both_copy_ctr][5]=$sub_end;
				$both_copies[$both_copy_ctr][6]=$ref_start; $both_copies[$both_copy_ctr++][7]=$ref_end;
			}
		}
		# Downstream: d5 (10000-15000 bases)
		elsif(($sub_start >= $ref_end+10000) && ($sub_start < $ref_end+15000)){
			if ($ref_strand eq $sub_strand){
				if($ref_strand eq '+'){# pos strand
					# check if the ref image has this relation with this family already
					if ((!(exists $pos_rel_history{"$ref_fam $ref_img $sub_fam d5"}) || ($pos_rel_history{"$ref_fam $ref_img $sub_fam d5"} ne $sub_fam)) && (!(exists $pos_rel_history{"$sub_fam $sub_img $ref_fam d5"}) || ($pos_rel_history{"$sub_fam $sub_img $ref_fam d5"} ne $ref_fam))) {
						# create history entry
						$pos_rel_history{"$ref_fam $ref_img $sub_fam d5"} = $sub_fam;
						$pos_rel_history{"$sub_fam $sub_img $ref_fam d5"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $pos_relationships{"$sub_fam $ref_fam d5"}) {
							$pos_relationships{"$sub_fam $ref_fam d5"} = 1;
						}
						else{
							$pos_relationships{"$sub_fam $ref_fam d5"}++;
						}
						# add record for d5 relationship
						$pos_copies[$pos_copy_ctr][0]=$sub_fam; $pos_copies[$pos_copy_ctr][1]="d5";
						$pos_copies[$pos_copy_ctr][2]=$ref_fam; $pos_copies[$pos_copy_ctr][3]='+';
						$pos_copies[$pos_copy_ctr][4]=$sub_start; $pos_copies[$pos_copy_ctr][5]=$sub_end;
						$pos_copies[$pos_copy_ctr][6]=$ref_start; $pos_copies[$pos_copy_ctr++][7]=$ref_end;
					}
				}
				elsif($ref_strand eq 'C'){# comp strand
					# check if the ref image has this relation with this family already
					if ((!(exists $comp_rel_history{"$ref_fam $ref_img $sub_fam d5"}) || ($comp_rel_history{"$ref_fam $ref_img $sub_fam d5"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_img $ref_fam d5"}) || ($comp_rel_history{"$sub_fam $sub_img $ref_fam d5"} ne $ref_fam))) {
						# create history entry
						$comp_rel_history{"$ref_fam $ref_img $sub_fam d5"} = $sub_fam;
						$comp_rel_history{"$sub_fam $sub_img $ref_fam d5"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $comp_relationships{"$ref_fam $sub_fam d5"}) {
							# now ref fam is upstream of sub fam as we are 
							# counting from right
							$comp_relationships{"$ref_fam $sub_fam d5"} = 1;
						}
						else{
							$comp_relationships{"$ref_fam $sub_fam d5"}++;
						}
						# add record for d5 relationship
						$comp_copies[$comp_copy_ctr][0]=$ref_fam; $comp_copies[$comp_copy_ctr][1]="d5";
						$comp_copies[$comp_copy_ctr][2]=$sub_fam; $comp_copies[$comp_copy_ctr][3]='C';
						$comp_copies[$comp_copy_ctr][4]=$ref_start; $comp_copies[$comp_copy_ctr][5]=$ref_end;
						$comp_copies[$comp_copy_ctr][6]=$sub_start; $comp_copies[$comp_copy_ctr++][7]=$sub_end;
					}
				}
			}
			# irrespective of strand
			# check if the ref image has this relation with this family already
			if ((!(exists $both_rel_history{"$ref_fam $ref_img $sub_fam d5"}) || ($both_rel_history{"$ref_fam $ref_img $sub_fam d5"} ne $sub_fam)) && (!(exists $both_rel_history{"$sub_fam $sub_img $ref_fam d5"}) || ($both_rel_history{"$sub_fam $sub_img $ref_fam d5"} ne $ref_fam))) {
				# create history entry
				$both_rel_history{"$ref_fam $ref_img $sub_fam d5"} = $sub_fam;
				$both_rel_history{"$sub_fam $sub_img $ref_fam d5"} = $ref_fam;
				
				# increment relationship count or create relationship entry
				if (!exists $both_relationships{"$sub_fam $ref_fam d5"}) {
					$both_relationships{"$sub_fam $ref_fam d5"} = 1;
				}
				else{
					$both_relationships{"$sub_fam $ref_fam d5"}++;
				}
				# add record for d5 relationship
				$both_copies[$both_copy_ctr][0]=$sub_fam; $both_copies[$both_copy_ctr][1]="d5";
				$both_copies[$both_copy_ctr][2]=$ref_fam; $both_copies[$both_copy_ctr][3]='B';
				$both_copies[$both_copy_ctr][4]=$sub_start; $both_copies[$both_copy_ctr][5]=$sub_end;
				$both_copies[$both_copy_ctr][6]=$ref_start; $both_copies[$both_copy_ctr++][7]=$ref_end;
			}
		}
		if($copy_file_flag){
			# temporary fix to reduce memory consumption when 
			# copies are not needed
			@pos_copies=();
			@comp_copies=();
			@both_copies=();# deallocating memory
			$pos_copy_ctr=$comp_copy_ctr=$both_copy_ctr=0;# resetting counters
		}
# 		print STDERR '.';
	}#end while
}# end relationship finding

#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\n\# Runtime details after finding relationships: \n";
print OUTFILEDATA "\# System time for process: $system_t\n";
print OUTFILEDATA "\# User time for process: $user_t\n\n\n";

print STDERR "\n\# Runtime details after finding relationships: \n";
print STDERR "\# System time for process: $system_t\n";
print STDERR "\# User time for process: $user_t\n\n";

# PRINTING THE ITEMSETS
# %pos_relationships : [fam1 fam2 category] = count
# %comp_relationships : [fam1 fam2 category] = count
# %both_relationships : [fam1 fam2 category] = count
# @count: fam occurences avg-len imagenum

# Creating the frequent itemsets in the format
# fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, Occurence, Strand, Category
print OUTFILEDATA "\# Total records in OUT file: $tot_recs\n";
print OUTFILEDATA "\# Total number of families: $tot_fams\n\n";
print OUTFILEDATA "\# Note: If dns rec ~ ups rec, then the regions were located uniformly\n";
print OUTFILEDATA "\# Average number of upstream OUT records processed per image: ".ceil($ups_ctr/$tot_recs)."\n";
print OUTFILEDATA "\# Average number of downstream OUT records processed per image: ".ceil($dns_ctr/$tot_recs)."\n";
print OUTFILEDATA "\# Average number of OUT records processed per image: ".ceil(($ups_ctr+$dns_ctr)/$tot_recs)."\n\n";
print OUTFILEDATA "\# Total relationships on pos strand:".keys(%pos_relationships)."\n";
print OUTFILEDATA "\# Total copies/clusters on pos strand:".$pos_copy_ctr."\n";
print OUTFILEDATA "\# Total relationships on comp strand:".keys(%comp_relationships)."\n";
print OUTFILEDATA "\# Total copies/clusters on comp strand:".$comp_copy_ctr."\n";
print OUTFILEDATA "\# Total relationships on both strands:".keys(%both_relationships)."\n";
print OUTFILEDATA "\# Total copies/clusters on both strands:".$both_copy_ctr."\n\n\n";

# TESTING
# relationships on the positive strand
# while( ($i,$j) = each %pos_relationships){
# 	@temp=split(' ',$i);
# 	print OUTFILEDATA "$temp[0]\t$temp[2]\t$temp[1]\t$j\t+\n";
# }
# 
# # relationships on the comp strand
# while( ($i,$j) = each %comp_relationships){
# 	@temp=split(' ',$i);
# 	print OUTFILEDATA "$temp[0]\t$temp[2]\t$temp[1]\t$j\tC\n";
# }


# relationships on the positive strand
while( ($i,$j) = each %pos_relationships){
	@temp=split(' ',$i);
	$rec=&get_index($temp[0]);
	print OUTFILEDATA "$temp[0]\t$counts[$rec][1]\t$counts[$rec][2]\t";
	$rec=&get_index($temp[1]);
	print OUTFILEDATA "$temp[1]\t$counts[$rec][1]\t$counts[$rec][2]\t";
	print OUTFILEDATA "$j\t+\t$temp[2]\n";
}

# relationships on the comp strand
while( ($i,$j) = each %comp_relationships){
	@temp=split(' ',$i);
	$rec=&get_index($temp[0]);
	print OUTFILEDATA "$temp[0]\t$counts[$rec][1]\t$counts[$rec][2]\t";
	$rec=&get_index($temp[1]);
	print OUTFILEDATA "$temp[1]\t$counts[$rec][1]\t$counts[$rec][2]\t";
	print OUTFILEDATA "$j\tC\t$temp[2]\n";
}
# relationships on both strands
while( ($i,$j) = each %both_relationships){
	@temp=split(' ',$i);
	$rec=&get_index($temp[0]);
	print OUTFILEDATA "$temp[0]\t$counts[$rec][1]\t$counts[$rec][2]\t";
	$rec=&get_index($temp[1]);
	print OUTFILEDATA "$temp[1]\t$counts[$rec][1]\t$counts[$rec][2]\t";
	print OUTFILEDATA "$j\tB\t$temp[2]\n";
}


#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\n\# Runtime details after printing itemsets: \n";
print OUTFILEDATA "\# System time for process: $system_t\n";
print OUTFILEDATA "\# User time for process: $user_t\n";

print STDERR "\n\# Runtime details after printing itemsets: \n";
print STDERR "\# System time for process: $system_t\n";
print STDERR "\# User time for process: $user_t\n";

# PRINTING  @copies 
# @copies : fam1 rel fam2 Strand fam1_st fam1_end fam2_st fam2_end
if($copy_file_flag){
	print OUTFILECOPIESPOS "#fam1\trel\tfam2\tStrand\tfam1-st\tf1-end\tf2-st\tf2-end\n";
	foreach $i (@pos_copies){
		print OUTFILECOPIESPOS "$i->[0]\t$i->[1]\t$i->[2]\t$i->[3]\t$i->[4]\t$i->[5]\t$i->[6]\t$i->[7]\n";
	}
	
	print OUTFILECOPIESCOMP "#fam1\trel\tfam2\tStrand\tfam1-st\tf1-end\tf2-st\tf2-end\n";
	foreach $i (@comp_copies){
		print OUTFILECOPIESCOMP "$i->[0]\t$i->[1]\t$i->[2]\t$i->[3]\t$i->[4]\t$i->[5]\t$i->[6]\t$i->[7]\n";
	}
	
	print OUTFILECOPIESBOTH "#fam1\trel\tfam2\tStrand\tfam1-st\tf1-end\tf2-st\tf2-end\n";
	foreach $i (@both_copies){
		print OUTFILECOPIESBOTH "$i->[0]\t$i->[1]\t$i->[2]\t$i->[3]\t$i->[4]\t$i->[5]\t$i->[6]\t$i->[7]\n";
	}


	#calculating time taken
	($user_t,$system_t,$cuser_t,$csystem_t) = times;
	print OUTFILECOPIESPOS "\n\# Runtime details after printing copy info: \n";
	print OUTFILECOPIESPOS "\# System time for process: $system_t\n";
	print OUTFILECOPIESPOS "\# User time for process: $user_t\n";
	print OUTFILECOPIESCOMP "\n\# Runtime details after printing copy info: \n";
	print OUTFILECOPIESCOMP "\# System time for process: $system_t\n";
	print OUTFILECOPIESCOMP "\# User time for process: $user_t\n";
	print OUTFILECOPIESBOTH "\n\# Runtime details after printing copy info: \n";
	print OUTFILECOPIESBOTH "\# System time for process: $system_t\n";
	print OUTFILECOPIESBOTH "\# User time for process: $user_t\n";
	
	print STDERR "\n\# Runtime details after printing copy info: \n";
	print STDERR "\# System time for process: $system_t\n";
	print STDERR "\# User time for process: $user_t\n";
}

close (INFILEDATA);
close (OUTFILEDATA);
if($copy_file_flag){
	close (OUTFILECOPIESPOS);
	close (OUTFILECOPIESCOMP);
	close (OUTFILECOPIESBOTH);
}
# # debugging
# close (ERRFILE);
exit;
