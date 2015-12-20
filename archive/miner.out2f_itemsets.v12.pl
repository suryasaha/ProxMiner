#!/usr/bin/perl -w
# MGEL
# Surya Saha 3/15/07
# reading cmd line input .out file which is sorted on the start position
# and finds the relationship among elements and families
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

# v5: count relations for elements INSTEAD of families
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

# v11: Modified get_index to use a hash instead of iterating thru an array
# v11: Improved runtime on chr12 from 25 mins to 2 mins
# v11: Fixed it so no information is recorded for copies unless required

# v12 : Read in the all ranges and relation names from config file, even for Ovlap

my $ver=12;

use strict;
use warnings;
use POSIX;

unless (@ARGV == 3){
	print "USAGE: $0 <config file> <RM .out file> <write copies??(0/1)>\n";
	exit;
}


my  ($ifname,$rec,@temp,%temphash,$ctr,$i,$j,$copy_file_flag);

my (@table,@famnames,@counts,%counts_index,$ups_ctr,$dns_ctr, $ref_img,
$ref_fam,$ref_start,$ref_end,$ref_strand,$sub_fam,$sub_img,$sub_start,
$sub_end,$sub_strand,%pos_relationships, %comp_relationships, %both_relationships, 
%pos_rel_history, %comp_rel_history, %both_rel_history, $tot_fams, 
$tot_recs,$user_t,$system_t, $cuser_t,$csystem_t);

my ($pos_copy_ctr,$comp_copy_ctr, $both_copy_ctr,
@pos_copies,@comp_copies,@both_copies, $ups_lim, $dns_lim,
@ups_rels, @dns_rels, @ovlap_rels, $rel);


$ifname=$ARGV[0];
chomp $ifname;
unless(open(INFILECONFIG,$ifname)){print "not able to open ".$ifname."\n\n";exit;}

$ifname=$ARGV[1];
chomp $ifname;
unless(open(INFILEDATA,$ifname)){print "not able to open ".$ifname."\n\n";exit;}
unless(open(OUTFILEDATA,">$ifname.f_itemsets.tab")){print "not able to open ".$ifname."f_itemsets.tab\n\n";exit;}

$copy_file_flag=$ARGV[2];
chomp $copy_file_flag;
if(!($copy_file_flag == 0 || $copy_file_flag == 1)){ print STDERR "flag can be only 0 or 1\nExiting..\n"; exit;}

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
# reimplementing the subroutine to use a hash to return the location of 
# the family, help in speedup??
#params: $fam
sub get_index{
	$_[0]=~ s/\s*//g;
	return $counts_index{$_[0]};
}

# Getting the relationships from the config file
# validate the ranges 
# get the U and D limits
# @ups_rels, @dns_rels, @ovlap_rels: Dir	Rel	Rng Strt	Rng End 
$ups_lim=$dns_lim=0;
while($rec=<INFILECONFIG>){
	if($rec =~ /#/){next;}
	if(length ($rec) < 10){next;}#for avoiding last line
	
	@temp = split(' ',$rec);
	
	# checking direction
	if($temp[0] ne 'U' && $temp[0] ne 'D' && $temp[0] ne 'O'){ 
		print STDERR "Direction can be only U, O or D\nExiting..\n"; exit;
	}
	# checking range
	if($temp[2] > $temp[3]){ 
		print STDERR "Incorrect coordinates for $temp[1] relationship \nExiting..\n"; exit;
	}
	
	# Errors not checked : 
	# if rel names are unique or not
	# if ranges are overlapping or not
	
	if($temp[0] eq 'U'){
		push @ups_rels, [@temp];
		if($temp[3] > $ups_lim){ $ups_lim=$temp[3];}
	}
	elsif($temp[0] eq 'D'){
		push @dns_rels, [@temp];
		if($temp[3] > $dns_lim){ $dns_lim=$temp[3];}
	}
	else{
		push @ovlap_rels, [@temp];
	}
}
close (INFILECONFIG);



# sort the rel arrays in ascending order
@temp = sort {($a->[3] <=> $b->[3])} @ups_rels;
@ups_rels = @temp;

@temp = sort {($a->[3] <=> $b->[3])} @dns_rels;
@dns_rels = @temp;

@temp = sort {($a->[3] <=> $b->[3])} @ovlap_rels;
@ovlap_rels = @temp;

# printing out the relationships its searching for
print STDERR "Upstream relationships\n";
foreach(@ups_rels) {print STDERR "$_->[0] $_->[1] $_->[2] $_->[3]\n";}
print STDERR "\n";
print STDERR "Downstream relationships\n";
foreach(@dns_rels) {print STDERR "$_->[0] $_->[1] $_->[2] $_->[3]\n";}
print STDERR "\n";
print STDERR "Overlap relationships\n";
foreach(@ovlap_rels) {print STDERR "$_->[0] $_->[1] $_->[2] $_->[3]\n";}
print STDERR "\n";

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
close (INFILEDATA);

print OUTFILEDATA "\# Version: $ver\n";
if($copy_file_flag){
	print OUTFILECOPIESPOS "\# Version: $ver\n";
	print OUTFILECOPIESCOMP "\# Version: $ver\n";
	print OUTFILECOPIESBOTH "\# Version: $ver\n";
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
print OUTFILEDATA "\# System time for process: ",ceil($system_t/60)," mins\n";
print OUTFILEDATA "\# User time for process: ",ceil($user_t/60)," mins\n";
print OUTFILEDATA "\n";


print STDERR "\# Runtime details after reading $tot_recs from file: \n";
print STDERR "\# System time for process: ",ceil($system_t/60)," mins\n";
print STDERR "\# User time for process: ",ceil($user_t/60)," mins\n";


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
# @count: fam occurences avg-len elementnum
$ctr=0;
foreach(@famnames){
	$counts[$ctr][0]=$_;
	#adding a value into the hash for family pos
	$counts_index{"$_"} = $ctr;
	#initializing all counters to 0
	$counts[$ctr][1]=0;#occurences
	$counts[$ctr][2]=0;#avg length
	$counts[$ctr++][3]=0;#number of elements (mined till now)
}
$tot_fams=$ctr;

# populating the @counts array
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
# where @table[][14]=element number
foreach (@table){
	# since @counts[][3] is initialized to 0
	$_->[14]=1+$counts[&get_index($_->[9])][3]++; 
}


#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\# Runtime details after preparing \@counts and appending \@table: \n";
print OUTFILEDATA "\# System time for process: ",ceil($system_t/60)," mins\n";
print OUTFILEDATA "\# User time for process: ",ceil($user_t/60)," mins\n";
print OUTFILEDATA "\n";

print STDERR "\# Runtime details after preparing \@counts and appending \@table: \n";
print STDERR "\# System time for process: ",ceil($system_t/60)," mins\n";
print STDERR "\# User time for process: ",ceil($user_t/60)," mins\n";
print STDERR "\n";

# FINDING ALL RELATIONS
# @table sorted on start position
# 1935 10.6  0.0  2.8 chr12 8936  9225 27748096 C  R=286 Unknown (0) 283   2   3
# 0    1     2    3   4     5     6    7        8  9     10      11  12    13 14
# @count: fam occurences avg-len elementnum
# finding the relationships
# %pos_relationships : [fam1 fam2 category] = count
# %comp_relationships : [fam1 fam2 category] = count
# %both_relationships : [fam1 fam2 category] = count
# %pos_rel_history : [fam1 fam1-img category] = fam2
# %comp_rel_history : [fam1 fam1-img category] = fam2
# %both_rel_history : [fam1 fam1-img category] = fam2

$ups_ctr=$dns_ctr=0;

if($copy_file_flag){
	$pos_copy_ctr=$comp_copy_ctr=$both_copy_ctr=0;
}

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
	
	# only look for relationships with elements located before it and
	# ending within $ups_lim(15k) bases before ref_start or anytime after it
	while(($j!=0) && ($table[$j-1][6] > $ref_start-$ups_lim)) {
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
					# check if the ref element has this relation with this family already
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
						if($copy_file_flag){
							$pos_copies[$pos_copy_ctr][0]=$ref_fam; $pos_copies[$pos_copy_ctr][1]="IN";
							$pos_copies[$pos_copy_ctr][2]=$sub_fam; $pos_copies[$pos_copy_ctr][3]='+';
							$pos_copies[$pos_copy_ctr][4]=$ref_start; $pos_copies[$pos_copy_ctr][5]=$ref_end;
							$pos_copies[$pos_copy_ctr][6]=$sub_start; $pos_copies[$pos_copy_ctr++][7]=$sub_end;
						}

						# increment reciprocal relationship count or create relationship entry
						if (!exists $pos_relationships{"$sub_fam $ref_fam CONT"}) {
							$pos_relationships{"$sub_fam $ref_fam CONT"} = 1;
						}
						else{
							$pos_relationships{"$sub_fam $ref_fam CONT"}++;
						}
						
						if($copy_file_flag){
							# add record for CONT relationship
							$pos_copies[$pos_copy_ctr][0]=$sub_fam; $pos_copies[$pos_copy_ctr][1]="CONT";
							$pos_copies[$pos_copy_ctr][2]=$ref_fam; $pos_copies[$pos_copy_ctr][3]='+';
							$pos_copies[$pos_copy_ctr][4]=$sub_start; $pos_copies[$pos_copy_ctr][5]=$sub_end;
							$pos_copies[$pos_copy_ctr][6]=$ref_start; $pos_copies[$pos_copy_ctr++][7]=$ref_end;
						}
					}
				}
				elsif($ref_strand eq 'C'){# comp strand
					# check if the ref element has this relation with this family already
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
						if($copy_file_flag){
							$comp_copies[$comp_copy_ctr][0]=$ref_fam; $comp_copies[$comp_copy_ctr][1]="IN";
							$comp_copies[$comp_copy_ctr][2]=$sub_fam; $comp_copies[$comp_copy_ctr][3]='C';
							$comp_copies[$comp_copy_ctr][4]=$ref_start; $comp_copies[$comp_copy_ctr][5]=$ref_end;
							$comp_copies[$comp_copy_ctr][6]=$sub_start; $comp_copies[$comp_copy_ctr++][7]=$sub_end;
						}
						# increment reciprocal relationship count or create relationship entry
						if (!exists $comp_relationships{"$sub_fam $ref_fam CONT"}) {
							$comp_relationships{"$sub_fam $ref_fam CONT"} = 1;
						}
						else{
							$comp_relationships{"$sub_fam $ref_fam CONT"}++;
						}
						
						# add record for CONT relationship
						if($copy_file_flag){
							$comp_copies[$comp_copy_ctr][0]=$sub_fam; $comp_copies[$comp_copy_ctr][1]="CONT";
							$comp_copies[$comp_copy_ctr][2]=$ref_fam; $comp_copies[$comp_copy_ctr][3]='C';
							$comp_copies[$comp_copy_ctr][4]=$sub_start; $comp_copies[$comp_copy_ctr][5]=$sub_end;
							$comp_copies[$comp_copy_ctr][6]=$ref_start; $comp_copies[$comp_copy_ctr++][7]=$ref_end;
						}
					}
				}
			}
			# irrespective of strand
			# check if the ref element has this relation with this family already
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
				if($copy_file_flag){
					$both_copies[$both_copy_ctr][0]=$ref_fam; $both_copies[$both_copy_ctr][1]="IN";
					$both_copies[$both_copy_ctr][2]=$sub_fam; $both_copies[$both_copy_ctr][3]='B';
					$both_copies[$both_copy_ctr][4]=$ref_start; $both_copies[$both_copy_ctr][5]=$ref_end;
					$both_copies[$both_copy_ctr][6]=$sub_start; $both_copies[$both_copy_ctr++][7]=$sub_end;
				}
				# increment reciprocal relationship count or create relationship entry
				if (!exists $both_relationships{"$sub_fam $ref_fam CONT"}) {
					$both_relationships{"$sub_fam $ref_fam CONT"} = 1;
				}
				else{
					$both_relationships{"$sub_fam $ref_fam CONT"}++;
				}
				
				# add record for CONT relationship
				if($copy_file_flag){
					$both_copies[$both_copy_ctr][0]=$sub_fam; $both_copies[$both_copy_ctr][1]="CONT";
					$both_copies[$both_copy_ctr][2]=$ref_fam; $both_copies[$both_copy_ctr][3]='B';
					$both_copies[$both_copy_ctr][4]=$sub_start; $both_copies[$both_copy_ctr][5]=$sub_end;
					$both_copies[$both_copy_ctr][6]=$ref_start; $both_copies[$both_copy_ctr++][7]=$ref_end;
				}
			}
		}# IN end
		# Overlap: If overlap is more than 10%?? of length of either family (Ovlap)
		# now if subject fam ends within the reference fam
		elsif (($sub_end > $ref_start) &&  ($sub_end < $ref_end)) {
			my ($ovlap, $ref_ovlap, $sub_ovlap);
			$ovlap = $sub_end - $ref_start;
			
			$ref_ovlap = ($ovlap / ($ref_end - $ref_start)) * 100;
			$sub_ovlap = ($ovlap / ($sub_end - $sub_start)) * 100;
			
			# for all given overlap relations
			# O Ovlap-10to30 10 30
			# O Ovlap-30to70 30 70
			# O Ovlap-70plus 70 100
			
			# debugging $ctr=0;
			foreach $rel (@ovlap_rels){
				my ($rel_name,$ubound, $lbound);
				$rel_name = $rel->[1];
				$lbound = $rel->[2];
				$ubound = $rel->[3];
				
				# cleaning up
				$rel_name=~ s/\s//g; 
				$lbound=~ s/\s//g;
				$ubound=~ s/\s//g;
				
				if ((($ref_ovlap > $lbound) && ($ref_ovlap <= $ubound)) || 
				(($sub_ovlap > $lbound) && ($sub_ovlap <= $ubound))) {
					# && $ctr==0) {
					# debugging $ctr++;
					if ($ref_strand eq $sub_strand){
						if($ref_strand eq '+'){# pos strand
							# check if the ref element has this relation with this family already
							if ((!(exists $pos_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"}) || ($pos_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"} ne $sub_fam)) && (!(exists $pos_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"}) || ($pos_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"} ne $ref_fam))) {
								# create history entry
								$pos_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"} = $sub_fam;
								$pos_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"} = $ref_fam;
								
								# increment relationship count or create relationship entry
								if (!exists $pos_relationships{"$sub_fam $ref_fam $rel_name"}) {
									$pos_relationships{"$sub_fam $ref_fam $rel_name"} = 1;
								}
								else{
									$pos_relationships{"$sub_fam $ref_fam $rel_name"}++;
								}
								
								# add record for relationship
								if($copy_file_flag){
									$pos_copies[$pos_copy_ctr][0]=$sub_fam; $pos_copies[$pos_copy_ctr][1]=$rel_name;
									$pos_copies[$pos_copy_ctr][2]=$ref_fam; $pos_copies[$pos_copy_ctr][3]='+';
									$pos_copies[$pos_copy_ctr][4]=$sub_start; $pos_copies[$pos_copy_ctr][5]=$sub_end;
									$pos_copies[$pos_copy_ctr][6]=$ref_start; $pos_copies[$pos_copy_ctr++][7]=$ref_end;
								}
							}
						}
						elsif($ref_strand eq 'C'){# comp strand
							# check if the ref element has this relation with this family already
							if ((!(exists $comp_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"}) || ($comp_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"}) || ($comp_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"} ne $ref_fam))) {
								# create history entry
								$comp_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"} = $sub_fam;
								$comp_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"} = $ref_fam;
								
								# increment relationship count or create relationship entry
								if (!exists $comp_relationships{"$sub_fam $ref_fam $rel_name"}) {
									$comp_relationships{"$sub_fam $ref_fam $rel_name"} = 1;
								}
								else{
									$comp_relationships{"$sub_fam $ref_fam $rel_name"}++;
								}
								
								# add record for relationship
								if($copy_file_flag){
									$comp_copies[$comp_copy_ctr][0]=$sub_fam; $comp_copies[$comp_copy_ctr][1]=$rel_name;
									$comp_copies[$comp_copy_ctr][2]=$ref_fam; $comp_copies[$comp_copy_ctr][3]='C';
									$comp_copies[$comp_copy_ctr][4]=$sub_start; $comp_copies[$comp_copy_ctr][5]=$sub_end;
									$comp_copies[$comp_copy_ctr][6]=$ref_start; $comp_copies[$comp_copy_ctr++][7]=$ref_end;
								}
							}
						}
					}
					# irrespective of strand
					# check if the ref element has this relation with this family already
					if ((!(exists $both_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"}) || ($both_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"} ne $sub_fam)) && (!(exists $both_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"}) || ($both_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"} ne $ref_fam))) {
						# create history entry
						$both_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"} = $sub_fam;
						$both_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"} = $ref_fam;
	
						
						# increment relationship count or create relationship entry
						if (!exists $both_relationships{"$sub_fam $ref_fam $rel_name"}) {
							$both_relationships{"$sub_fam $ref_fam $rel_name"} = 1;
						}
						else{
							$both_relationships{"$sub_fam $ref_fam $rel_name"}++;
						}
						
						# add record for Ovlap-10to30 relationship
						if($copy_file_flag){
							$both_copies[$both_copy_ctr][0]=$sub_fam; $both_copies[$both_copy_ctr][1]=$rel_name;
							$both_copies[$both_copy_ctr][2]=$ref_fam; $both_copies[$both_copy_ctr][3]='B';
							$both_copies[$both_copy_ctr][4]=$sub_start; $both_copies[$both_copy_ctr][5]=$sub_end;
							$both_copies[$both_copy_ctr][6]=$ref_start; $both_copies[$both_copy_ctr++][7]=$ref_end;
						}
					}
				}
			}# foreach ovlap end
		}# Ovlap end
		# All Upstream relationships
		elsif(($sub_end <= $ref_start) && ($sub_end > $ref_start-$ups_lim)){
			foreach $rel (@ups_rels){
				my ($rel_name,$ubound, $lbound);
				$rel_name = $rel->[1];
				$lbound = $rel->[2];
				$ubound = $rel->[3];
				
				# cleaning up
				$rel_name=~ s/\s//g; 
				$lbound=~ s/\s//g;
				$ubound=~ s/\s//g;
				
				# for all given overlap relations
				# U u1 0 500
				# U u2 500 1000
				# U u3 1000 5000
				# U u4 5000 10000
				# U u5 10000 15000

				if(($sub_end <= $ref_start-$lbound) && ($sub_end > $ref_start-$ubound)){
					if ($ref_strand eq $sub_strand){
						if($ref_strand eq '+'){# pos strand
							# check if the ref element has this relation with this family already
							if ((!(exists $pos_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"}) || ($pos_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"} ne $sub_fam)) && (!(exists $pos_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"}) || ($pos_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"} ne $ref_fam))) {
								# create history entry
								$pos_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"} = $sub_fam;
								$pos_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"} = $ref_fam;
								
								# increment relationship count or create relationship entry
								if (!exists $pos_relationships{"$sub_fam $ref_fam $rel_name"}) {
									$pos_relationships{"$sub_fam $ref_fam $rel_name"} = 1;
								}
								else{
									$pos_relationships{"$sub_fam $ref_fam $rel_name"}++;
								}
								# add record for u1 relationship
								if($copy_file_flag){
									$pos_copies[$pos_copy_ctr][0]=$sub_fam; $pos_copies[$pos_copy_ctr][1]=$rel_name;
									$pos_copies[$pos_copy_ctr][2]=$ref_fam; $pos_copies[$pos_copy_ctr][3]='+';
									$pos_copies[$pos_copy_ctr][4]=$sub_start; $pos_copies[$pos_copy_ctr][5]=$sub_end;
									$pos_copies[$pos_copy_ctr][6]=$ref_start; $pos_copies[$pos_copy_ctr++][7]=$ref_end;
								}
							}
						}
						elsif($ref_strand eq 'C'){# comp strand
							# check if the ref element has this relation with this family already
							if ((!(exists $comp_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"}) || ($comp_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"}) || ($comp_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"} ne $ref_fam))) {
								# create history entry
								$comp_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"} = $sub_fam;
								$comp_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"} = $ref_fam;
								
								# increment relationship count or create relationship entry
								if (!exists $comp_relationships{"$ref_fam $sub_fam $rel_name"}) {
									# now ref fam is upstream of sub fam as we are 
									# counting from right
									$comp_relationships{"$ref_fam $sub_fam $rel_name"} = 1;
								}
								else{
									$comp_relationships{"$ref_fam $sub_fam $rel_name"}++;
								}
								
								# add record for u1 relationship
								if($copy_file_flag){
									$comp_copies[$comp_copy_ctr][0]=$ref_fam; $comp_copies[$comp_copy_ctr][1]=$rel_name;
									$comp_copies[$comp_copy_ctr][2]=$sub_fam; $comp_copies[$comp_copy_ctr][3]='C';
									$comp_copies[$comp_copy_ctr][4]=$ref_start; $comp_copies[$comp_copy_ctr][5]=$ref_end;
									$comp_copies[$comp_copy_ctr][6]=$sub_start; $comp_copies[$comp_copy_ctr++][7]=$sub_end;
								}
							}
						}
					}
					# irrespective of strand
					# check if the ref element has this relation with this family already
					if ((!(exists $both_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"}) || ($both_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"} ne $sub_fam)) && (!(exists $both_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"}) || ($both_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"} ne $ref_fam))) {
						# create history entry
						$both_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"} = $sub_fam;
						$both_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $both_relationships{"$sub_fam $ref_fam $rel_name"}) {
							$both_relationships{"$sub_fam $ref_fam $rel_name"} = 1;
						}
						else{
							$both_relationships{"$sub_fam $ref_fam $rel_name"}++;
						}
						
						# add record for u1 relationship
						if($copy_file_flag){
							$both_copies[$both_copy_ctr][0]=$sub_fam; $both_copies[$both_copy_ctr][1]=$rel_name;
							$both_copies[$both_copy_ctr][2]=$ref_fam; $both_copies[$both_copy_ctr][3]='B';
							$both_copies[$both_copy_ctr][4]=$sub_start; $both_copies[$both_copy_ctr][5]=$sub_end;
							$both_copies[$both_copy_ctr][6]=$ref_start; $both_copies[$both_copy_ctr++][7]=$ref_end;
						}
					}
				}# endif range
			}# foreach upstream end
		}# upstream rel end
	}# end while	
	
	$j=$i;
	
	# only look for relationships with elements located after it
	# and starting within $dns_lim(15k) bases after ref_end (enforced by condition above)
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
		
		# In & Contains: Location of ref fam is entirely within sub fam (IN,CONT)
		# IN should be first bcos if sub start is near the ref start, it will
		# be listed right after the ref record in the list
		if(($sub_start == $ref_start) && ($sub_end >= $ref_end)){
			if ($ref_strand eq $sub_strand){
				if($ref_strand eq '+'){# pos strand
					# check if the ref element has this relation with this family already
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
						if($copy_file_flag){
							$pos_copies[$pos_copy_ctr][0]=$ref_fam; $pos_copies[$pos_copy_ctr][1]="IN";
							$pos_copies[$pos_copy_ctr][2]=$sub_fam; $pos_copies[$pos_copy_ctr][3]='+';
							$pos_copies[$pos_copy_ctr][4]=$ref_start; $pos_copies[$pos_copy_ctr][5]=$ref_end;
							$pos_copies[$pos_copy_ctr][6]=$sub_start; $pos_copies[$pos_copy_ctr++][7]=$sub_end;
						}
						
						# increment reciprocal relationship count or create relationship entry
						if (!exists $pos_relationships{"$sub_fam $ref_fam CONT"}) {
							$pos_relationships{"$sub_fam $ref_fam CONT"} = 1;
						}
						else{
							$pos_relationships{"$sub_fam $ref_fam CONT"}++;
						}
						
						# add record for CONT relationship
						if($copy_file_flag){
							$pos_copies[$pos_copy_ctr][0]=$sub_fam; $pos_copies[$pos_copy_ctr][1]="CONT";
							$pos_copies[$pos_copy_ctr][2]=$ref_fam; $pos_copies[$pos_copy_ctr][3]='+';
							$pos_copies[$pos_copy_ctr][4]=$sub_start; $pos_copies[$pos_copy_ctr][5]=$sub_end;
							$pos_copies[$pos_copy_ctr][6]=$ref_start; $pos_copies[$pos_copy_ctr++][7]=$ref_end;
						}
					}
				}
				elsif($ref_strand eq 'C'){# comp strand
					# check if the ref element has this relation with this family already
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
						if($copy_file_flag){
							$comp_copies[$comp_copy_ctr][0]=$ref_fam; $comp_copies[$comp_copy_ctr][1]="IN";
							$comp_copies[$comp_copy_ctr][2]=$sub_fam; $comp_copies[$comp_copy_ctr][3]='C';
							$comp_copies[$comp_copy_ctr][4]=$ref_start; $comp_copies[$comp_copy_ctr][5]=$ref_end;
							$comp_copies[$comp_copy_ctr][6]=$sub_start; $comp_copies[$comp_copy_ctr++][7]=$sub_end;
						}
						
						# increment reciprocal relationship count or create relationship entry
						if (!exists $comp_relationships{"$sub_fam $ref_fam CONT"}) {
							$comp_relationships{"$sub_fam $ref_fam CONT"} = 1;
						}
						else{
							$comp_relationships{"$sub_fam $ref_fam CONT"}++;
						}
						
						# add record for CONT relationship
						if($copy_file_flag){
							$comp_copies[$comp_copy_ctr][0]=$sub_fam; $comp_copies[$comp_copy_ctr][1]="CONT";
							$comp_copies[$comp_copy_ctr][2]=$ref_fam; $comp_copies[$comp_copy_ctr][3]='C';
							$comp_copies[$comp_copy_ctr][4]=$sub_start; $comp_copies[$comp_copy_ctr][5]=$sub_end;
							$comp_copies[$comp_copy_ctr][6]=$ref_start; $comp_copies[$comp_copy_ctr++][7]=$ref_end;
						}
					}
				}
			}
			# irrespective of strand
			# check if the ref element has this relation with this family already
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
				if($copy_file_flag){
					$both_copies[$both_copy_ctr][0]=$ref_fam; $both_copies[$both_copy_ctr][1]="IN";
					$both_copies[$both_copy_ctr][2]=$sub_fam; $both_copies[$both_copy_ctr][3]='B';
					$both_copies[$both_copy_ctr][4]=$ref_start; $both_copies[$both_copy_ctr][5]=$ref_end;
					$both_copies[$both_copy_ctr][6]=$sub_start; $both_copies[$both_copy_ctr++][7]=$sub_end;
				}

				# increment reciprocal relationship count or create relationship entry
				if (!exists $both_relationships{"$sub_fam $ref_fam CONT"}) {
					$both_relationships{"$sub_fam $ref_fam CONT"} = 1;
				}
				else{
					$both_relationships{"$sub_fam $ref_fam CONT"}++;
				}
				
				# add record for CONT relationship
				if($copy_file_flag){
					$both_copies[$both_copy_ctr][0]=$sub_fam; $both_copies[$both_copy_ctr][1]="CONT";
					$both_copies[$both_copy_ctr][2]=$ref_fam; $both_copies[$both_copy_ctr][3]='B';
					$both_copies[$both_copy_ctr][4]=$sub_start; $both_copies[$both_copy_ctr][5]=$sub_end;
					$both_copies[$both_copy_ctr][6]=$ref_start; $both_copies[$both_copy_ctr++][7]=$ref_end;
				}
			}
		}# IN end
		# Overlap: If overlap is more than 10%?? of length of either family (Ovlap)
		# now if subject fam ends within the reference fam
		elsif (($sub_start > $ref_start) &&  ($sub_start < $ref_end)) {
			my ($ovlap, $ref_ovlap, $sub_ovlap);
			$ovlap = $ref_end - $sub_start;
			
			$ref_ovlap = ($ovlap / ($ref_end - $ref_start)) * 100;
			$sub_ovlap = ($ovlap / ($sub_end - $sub_start)) * 100;

			# for all given overlap relations
			# O Ovlap-10to30 10 30
			# O Ovlap-30to70 30 70
			# O Ovlap-70plus 70 100
			
			# debugging $ctr=0;
			foreach $rel (@ovlap_rels){
				my ($rel_name,$ubound, $lbound);
				$rel_name = $rel->[1];
				$lbound = $rel->[2];
				$ubound = $rel->[3];
				
				# cleaning up
				$rel_name=~ s/\s//g; 
				$lbound=~ s/\s//g;
				$ubound=~ s/\s//g;
				
				if ((($ref_ovlap > $lbound) && ($ref_ovlap <= $ubound)) || 
				(($sub_ovlap > $lbound) && ($sub_ovlap <= $ubound))) {
					# && $ctr==0) {
					# $ctr++;
					if ($ref_strand eq $sub_strand){
						if($ref_strand eq '+'){# pos strand
							# check if the ref element has this relation with this family already
							if ((!(exists $pos_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"}) || ($pos_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"} ne $sub_fam)) && (!(exists $pos_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"}) || ($pos_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"} ne $ref_fam))) {
								# create history entry
								$pos_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"} = $sub_fam;
								$pos_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"} = $ref_fam;
								
								# increment relationship count or create relationship entry
								if (!exists $pos_relationships{"$sub_fam $ref_fam $rel_name"}) {
									$pos_relationships{"$sub_fam $ref_fam $rel_name"} = 1;
								}
								else{
									$pos_relationships{"$sub_fam $ref_fam $rel_name"}++;
								}
								
								# add record for relationship
								if($copy_file_flag){
									$pos_copies[$pos_copy_ctr][0]=$sub_fam; $pos_copies[$pos_copy_ctr][1]=$rel_name;
									$pos_copies[$pos_copy_ctr][2]=$ref_fam; $pos_copies[$pos_copy_ctr][3]='+';
									$pos_copies[$pos_copy_ctr][4]=$sub_start; $pos_copies[$pos_copy_ctr][5]=$sub_end;
									$pos_copies[$pos_copy_ctr][6]=$ref_start; $pos_copies[$pos_copy_ctr++][7]=$ref_end;
								}
							}
						}
						elsif($ref_strand eq 'C'){# comp strand
							# check if the ref element has this relation with this family already
							if ((!(exists $comp_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"}) || ($comp_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"}) || ($comp_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"} ne $ref_fam))) {
								# create history entry
								$comp_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"} = $sub_fam;
								$comp_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"} = $ref_fam;
								
								# increment relationship count or create relationship entry
								if (!exists $comp_relationships{"$sub_fam $ref_fam $rel_name"}) {
									$comp_relationships{"$sub_fam $ref_fam $rel_name"} = 1;
								}
								else{
									$comp_relationships{"$sub_fam $ref_fam $rel_name"}++;
								}
								
								# add record for relationship
								if($copy_file_flag){
									$comp_copies[$comp_copy_ctr][0]=$sub_fam; $comp_copies[$comp_copy_ctr][1]=$rel_name;
									$comp_copies[$comp_copy_ctr][2]=$ref_fam; $comp_copies[$comp_copy_ctr][3]='C';
									$comp_copies[$comp_copy_ctr][4]=$sub_start; $comp_copies[$comp_copy_ctr][5]=$sub_end;
									$comp_copies[$comp_copy_ctr][6]=$ref_start; $comp_copies[$comp_copy_ctr++][7]=$ref_end;
								}
							}
						}
					}
					# irrespective of strand
					# check if the ref element has this relation with this family already
					if ((!(exists $both_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"}) || ($both_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"} ne $sub_fam)) && (!(exists $both_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"}) || ($both_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"} ne $ref_fam))) {
						# create history entry
						$both_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"} = $sub_fam;
						$both_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"} = $ref_fam;
							
						# increment relationship count or create relationship entry
						if (!exists $both_relationships{"$sub_fam $ref_fam $rel_name"}) {
							$both_relationships{"$sub_fam $ref_fam $rel_name"} = 1;
						}
						else{
							$both_relationships{"$sub_fam $ref_fam $rel_name"}++;
						}
						
						# add record for Ovlap-10to30 relationship
						if($copy_file_flag){
							$both_copies[$both_copy_ctr][0]=$sub_fam; $both_copies[$both_copy_ctr][1]=$rel_name;
							$both_copies[$both_copy_ctr][2]=$ref_fam; $both_copies[$both_copy_ctr][3]='B';
							$both_copies[$both_copy_ctr][4]=$sub_start; $both_copies[$both_copy_ctr][5]=$sub_end;
							$both_copies[$both_copy_ctr][6]=$ref_start; $both_copies[$both_copy_ctr++][7]=$ref_end;
						}
					}
				}
			}# foreach ovlap end
		}# Ovlap rel end
		# All Downstream relationships
		elsif(($sub_start >= $ref_end) && ($sub_start < $ref_end+$dns_lim)){
			foreach $rel (@dns_rels){
				my ($rel_name,$ubound, $lbound);
				$rel_name = $rel->[1];
				$lbound = $rel->[2];
				$ubound = $rel->[3];
				
				# cleaning up
				$rel_name=~ s/\s//g; 
				$lbound=~ s/\s//g;
				$ubound=~ s/\s//g;
				
				# for all given relationships
				# D d1 0 500
				# D d2 500 1000
				# D d3 1000 5000
				# D d4 5000 10000
				# D d5 10000 15000
				
				# check if within range
				if(($sub_start >= $ref_end+$lbound) && ($sub_start < $ref_end+$ubound)){
					if ($ref_strand eq $sub_strand){
						if($ref_strand eq '+'){# pos strand
							# check if the ref element has this relation with this family already
							if ((!(exists $pos_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"}) || ($pos_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"} ne $sub_fam)) && (!(exists $pos_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"}) || ($pos_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"} ne $ref_fam))) {
								# create history entry
								$pos_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"} = $sub_fam;
								$pos_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"} = $ref_fam;
								
								# increment relationship count or create relationship entry
								if (!exists $pos_relationships{"$sub_fam $ref_fam $rel_name"}) {
									$pos_relationships{"$sub_fam $ref_fam $rel_name"} = 1;
								}
								else{
									$pos_relationships{"$sub_fam $ref_fam $rel_name"}++;
								}
								
								# add record for d1 relationship
								if($copy_file_flag){
									$pos_copies[$pos_copy_ctr][0]=$sub_fam; $pos_copies[$pos_copy_ctr][1]=$rel_name;
									$pos_copies[$pos_copy_ctr][2]=$ref_fam; $pos_copies[$pos_copy_ctr][3]='+';
									$pos_copies[$pos_copy_ctr][4]=$sub_start; $pos_copies[$pos_copy_ctr][5]=$sub_end;
									$pos_copies[$pos_copy_ctr][6]=$ref_start; $pos_copies[$pos_copy_ctr++][7]=$ref_end;
								}
							}
						}
						elsif($ref_strand eq 'C'){# comp strand
							# check if the ref element has this relation with this family already
							if ((!(exists $comp_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"}) || ($comp_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"}) || ($comp_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"} ne $ref_fam))) {
								# create history entry
								$comp_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"} = $sub_fam;
								$comp_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"} = $ref_fam;
								
								# increment relationship count or create relationship entry
								if (!exists $comp_relationships{"$ref_fam $sub_fam $rel_name"}) {
									# now ref fam is upstream of sub fam as we are 
									# counting from right
									$comp_relationships{"$ref_fam $sub_fam $rel_name"} = 1;
								}
								else{
									$comp_relationships{"$ref_fam $sub_fam $rel_name"}++;
								}
								
								# add record for d1 relationship
								if($copy_file_flag){
									$comp_copies[$comp_copy_ctr][0]=$ref_fam; $comp_copies[$comp_copy_ctr][1]=$rel_name;
									$comp_copies[$comp_copy_ctr][2]=$sub_fam; $comp_copies[$comp_copy_ctr][3]='C';
									$comp_copies[$comp_copy_ctr][4]=$ref_start; $comp_copies[$comp_copy_ctr][5]=$ref_end;
									$comp_copies[$comp_copy_ctr][6]=$sub_start; $comp_copies[$comp_copy_ctr++][7]=$sub_end;
								}
							}
						}
					}
					# irrespective of strand
					# check if the ref element has this relation with this family already
					if ((!(exists $both_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"}) || ($both_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"} ne $sub_fam)) && (!(exists $both_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"}) || ($both_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"} ne $ref_fam))) {
						# create history entry
						$both_rel_history{"$ref_fam $ref_img $sub_fam $rel_name"} = $sub_fam;
						$both_rel_history{"$sub_fam $sub_img $ref_fam $rel_name"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $both_relationships{"$sub_fam $ref_fam $rel_name"}) {
							$both_relationships{"$sub_fam $ref_fam $rel_name"} = 1;
						}
						else{
							$both_relationships{"$sub_fam $ref_fam $rel_name"}++;
						}
						
						# add record for d1 relationship
						if($copy_file_flag){
							$both_copies[$both_copy_ctr][0]=$sub_fam; $both_copies[$both_copy_ctr][1]=$rel_name;
							$both_copies[$both_copy_ctr][2]=$ref_fam; $both_copies[$both_copy_ctr][3]='B';
							$both_copies[$both_copy_ctr][4]=$sub_start; $both_copies[$both_copy_ctr][5]=$sub_end;
							$both_copies[$both_copy_ctr][6]=$ref_start; $both_copies[$both_copy_ctr++][7]=$ref_end;
						}
					}
				}# endif range check 
			}# end foreach Downstream rels
		}# end Downstream relationships
	}# end while	
}# end relationship finding


#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\n\# Runtime details after finding relationships: \n";
print OUTFILEDATA "\# System time for process: ",ceil($system_t/60)," mins\n";
print OUTFILEDATA "\# User time for process: ",ceil($user_t/60)," mins\n\n\n";

print STDERR "\n\# Runtime details after finding relationships: \n";
print STDERR "\# System time for process: ",ceil($system_t/60)," mins\n";
print STDERR "\# User time for process: ",ceil($user_t/60)," mins\n\n";

# PRINTING THE ITEMSETS
# %pos_relationships : [fam1 fam2 category] = count
# %comp_relationships : [fam1 fam2 category] = count
# %both_relationships : [fam1 fam2 category] = count
# @count: fam occurences avg-len elementnum

# Creating the frequent itemsets in the format
# fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, Occurence, Strand, Category
print OUTFILEDATA "\# Total records in OUT file: $tot_recs\n";
print OUTFILEDATA "\# Total number of families: $tot_fams\n\n";
print OUTFILEDATA "\# Note: If dns rec ~ ups rec, then the regions were located uniformly\n";
print OUTFILEDATA "\# Average number of upstream OUT records processed per element: ".ceil($ups_ctr/$tot_recs)."\n";
print OUTFILEDATA "\# Average number of downstream OUT records processed per element: ".ceil($dns_ctr/$tot_recs)."\n";
print OUTFILEDATA "\# Average number of OUT records processed per element: ".ceil(($ups_ctr+$dns_ctr)/$tot_recs)."\n\n";
print OUTFILEDATA "\# Total relationships on pos strand:".keys(%pos_relationships)."\n";
if($copy_file_flag){ print OUTFILEDATA "\# Total copies/clusters on pos strand:".$pos_copy_ctr."\n";}
print OUTFILEDATA "\# Total relationships on comp strand:".keys(%comp_relationships)."\n";
if($copy_file_flag){ print OUTFILEDATA "\# Total copies/clusters on comp strand:".$comp_copy_ctr."\n";}
print OUTFILEDATA "\# Total relationships on both strands:".keys(%both_relationships)."\n";
if($copy_file_flag){ print OUTFILEDATA "\# Total copies/clusters on both strands:".$both_copy_ctr."\n\n\n";}
else{print OUTFILEDATA "\n\n";}

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
#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\n\# Runtime details after printing positive itemsets: \n";
print OUTFILEDATA "\# System time for process: ",ceil($system_t/60)," mins\n";
print OUTFILEDATA "\# User time for process: ",ceil($user_t/60)," mins\n";

print STDERR "\n\# Runtime details after printing positive itemsets: \n";
print STDERR "\# System time for process: ",ceil($system_t/60)," mins\n";
print STDERR "\# User time for process: ",ceil($user_t/60)," mins\n";


# relationships on the comp strand
while( ($i,$j) = each %comp_relationships){
	@temp=split(' ',$i);
	$rec=&get_index($temp[0]);
	print OUTFILEDATA "$temp[0]\t$counts[$rec][1]\t$counts[$rec][2]\t";
	$rec=&get_index($temp[1]);
	print OUTFILEDATA "$temp[1]\t$counts[$rec][1]\t$counts[$rec][2]\t";
	print OUTFILEDATA "$j\tC\t$temp[2]\n";
}

#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\n\# Runtime details after printing negative itemsets: \n";
print OUTFILEDATA "\# System time for process: ",ceil($system_t/60)," mins\n";
print OUTFILEDATA "\# User time for process: ",ceil($user_t/60)," mins\n";

print STDERR "\n\# Runtime details after printing negative itemsets: \n";
print STDERR "\# System time for process: ",ceil($system_t/60)," mins\n";
print STDERR "\# User time for process: ",ceil($user_t/60)," mins\n";

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
print OUTFILEDATA "\n\# Runtime details after printing both itemsets: \n";
print OUTFILEDATA "\# System time for process: ",ceil($system_t/60)," mins\n";
print OUTFILEDATA "\# User time for process: ",ceil($user_t/60)," mins\n";

print STDERR "\n\# Runtime details after printing both itemsets: \n";
print STDERR "\# System time for process: ",ceil($system_t/60)," mins\n";
print STDERR "\# User time for process: ",ceil($user_t/60)," mins\n";

# PRINTING  @copies 
# @copies : fam1 rel fam2 Strand fam1_st fam1_end fam2_st fam2_end
if($copy_file_flag){
	print OUTFILECOPIESPOS "#fam1\trel\tfam2\tstrand\tfam1-st\tf1-end\tf2-st\tf2-end\n";
	foreach $i (@pos_copies){
		print OUTFILECOPIESPOS "$i->[0]\t$i->[1]\t$i->[2]\t$i->[3]\t$i->[4]\t$i->[5]\t$i->[6]\t$i->[7]\n";
	}
	
	print OUTFILECOPIESCOMP "#fam1\trel\tfam2\tstrand\tfam1-st\tf1-end\tf2-st\tf2-end\n";
	foreach $i (@comp_copies){
		print OUTFILECOPIESCOMP "$i->[0]\t$i->[1]\t$i->[2]\t$i->[3]\t$i->[4]\t$i->[5]\t$i->[6]\t$i->[7]\n";
	}
	
	print OUTFILECOPIESBOTH "#fam1\trel\tfam2\tstrand\tfam1-st\tf1-end\tf2-st\tf2-end\n";
	foreach $i (@both_copies){
		print OUTFILECOPIESBOTH "$i->[0]\t$i->[1]\t$i->[2]\t$i->[3]\t$i->[4]\t$i->[5]\t$i->[6]\t$i->[7]\n";
	}


	#calculating time taken
	($user_t,$system_t,$cuser_t,$csystem_t) = times;
	print OUTFILECOPIESPOS "\n\# Runtime details after printing copy info: \n";
	print OUTFILECOPIESPOS "\# System time for process: ",ceil($system_t/60)," mins\n";
	print OUTFILECOPIESPOS "\# User time for process: ",ceil($user_t/60)," mins\n";
	print OUTFILECOPIESCOMP "\n\# Runtime details after printing copy info: \n";
	print OUTFILECOPIESCOMP "\# System time for process: ",ceil($system_t/60)," mins\n";
	print OUTFILECOPIESCOMP "\# User time for process: ",ceil($user_t/60)," mins\n";
	print OUTFILECOPIESBOTH "\n\# Runtime details after printing copy info: \n";
	print OUTFILECOPIESBOTH "\# System time for process: ",ceil($system_t/60)," mins\n";
	print OUTFILECOPIESBOTH "\# User time for process: ",ceil($user_t/60)," mins\n";
	
	print STDERR "\n\# Runtime details after printing copy info: \n";
	print STDERR "\# System time for process: ",ceil($system_t/60)," mins\n";
	print STDERR "\# User time for process: ",ceil($user_t/60)," mins\n";
}

close (OUTFILEDATA);
if($copy_file_flag){
	close (OUTFILECOPIESPOS);
	close (OUTFILECOPIESCOMP);
	close (OUTFILECOPIESBOTH);
}
# # debugging
# close (ERRFILE);
exit;
