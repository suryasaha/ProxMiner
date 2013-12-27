#!/usr/bin/perl -w
# MGEL
# Surya Saha 02/02/08
# reading cmd line input .out file which is sorted on the start position
# and finds the relationship among elements and families
# Relationship types: 
# Upstream: 
#	u (0-MAX bases)
#	spl case Contains: Location of fam2 is entirely within fam1 (Cont)
# Downstream: 
#	d (0-MAX bases)
#	spl case In: Location of fam2 is entirely within fam1 (IN)
# Overlap: 

# OUTPUT: Creating the frequent itemsets in the format
# fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, Occurence, Strand, Category, Avg. dist., Std. dev.
# Printing "NA" for Avg. dist., Std. dev. in case of IN/CONT relations

# v1: Using only 1 MAX range for U and D relation. No other kind of relation.
#   : No config file required. Only 1 output file for combined pos and comp strand information. 
#   : No output for both strands (irrespective of strand)
#   : Added the avg and std. dev. for each family pair in a relationship
#   : Not accounting for IN/CONT type rels in std devs
#   NOTE: Deviations and rel. counts are different if any IN relationships are present, but only for the + strand. 
#   No diff's for the Comp strand as a new UPS relation is created and the copies are barred from creating any more #   relations. See doc./paper.


my $ver=1;

use strict;
use warnings;
use POSIX;
# use Math::Complex;

unless (@ARGV == 3){
	print "USAGE: $0 <MAX_RANGE> <RM .out file> <write copies??(-nocop/-cop)>\n";
	exit;
}


my  ($ifname,$rec,@temp,%temphash,$ctr,$i,$j,$k,$copy_file_flag,$MAX_RANGE);

my (@table,@famnames,@counts,%counts_index,$ups_ctr,$dns_ctr, $ref_elem,$ref_fam,$ref_start,$ref_end,$ref_strand,
$sub_fam,$sub_elem,$sub_start,$sub_end,$sub_strand,%pos_relationships,%comp_relationships,$pos_copy_ctr,
$comp_copy_ctr,@pos_copies,@comp_copies,%pos_rel_history,%comp_rel_history,$tot_fams,$tot_recs,
$user_t,$system_t,$cuser_t,$csystem_t);

my (%reci_rels, $rel, @devs, @old_devs, %std_dev, $avg);


$MAX_RANGE=$ARGV[0];
chomp $MAX_RANGE;

$ifname=$ARGV[1];
chomp $ifname;
unless(open(INFILEDATA,$ifname)){print "not able to open ".$ifname."\n\n";exit;}
unless(open(OUTFILEDATA,">$ifname.f_itemsets.stage1.tab")){print "not able to open ".$ifname."f_itemsets.stage1.tab\n\n";exit;}

$copy_file_flag=$ARGV[2];
chomp $copy_file_flag;
if(!($copy_file_flag eq "-nocop" || $copy_file_flag eq "-cop")){ print STDERR "flag can be only -cop or -nocop\nExiting..\n"; exit;}
# convert to 0 or 1
if($copy_file_flag eq "-nocop"){ $copy_file_flag=0;}
elsif ($copy_file_flag eq "-cop"){ $copy_file_flag=1;}

if($copy_file_flag){
	unless(open(OUTFILECOPIESPOS,">$ifname.pos.copies.tab")){print "not able to open ".$ifname.".copies.tab \n\n";exit;}
	unless(open(OUTFILECOPIESCOMP,">$ifname.comp.copies.tab")){print "not able to open ".$ifname."comp.copies.tab \n\n"; exit;}
}

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


sub round2{
	my ($num);
	$num=$_[0];
	$num=$num*100;
	$num=int($num);
	$num=$num/100;
	return $num;
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
close (INFILEDATA);

print OUTFILEDATA "\# Version: $ver\n";
if($copy_file_flag){
	print OUTFILECOPIES "\# Version: $ver\n";
}
$i=localtime();
print OUTFILEDATA "\# Time: $i\n";
if($copy_file_flag){
	print OUTFILECOPIES "\# Time: $i\n";
}

#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\# Runtime details after reading in the file: \n";
print OUTFILEDATA "\# System time for process: ",ceil($system_t/60)," mins\n";
print OUTFILEDATA "\# User time for process: ",ceil($user_t/60)," mins\n";
print OUTFILEDATA "\n";


print STDERR "\# Runtime details after reading $tot_recs records from file: \n";
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
# %pos_rel_history : [fam1 fam1-img fam2 category] = fam2
# %comp_rel_history : [fam1 fam1-img fam2 category] = fam2
# Feb 1, 08 : I do not remember why we have 
# %pos_rel_history : [fam1 fam1-img fam2 category] = fam2 in the code???
# but testing shows that only [fam1 fam1-img category] produces conf > 1.0
# %std_dev : ["fam1 fam2 rel strand"] => [array of deviations]

$ups_ctr=$dns_ctr=0;

for $i (0 .. $#table){
	$ref_start=$table[$i][5]; $ref_end=$table[$i][6];
	$ref_strand=$table[$i][8]; $ref_fam=$table[$i][9];
	$ref_elem=$table[$i][14];
	
	# cleaning up
	$ref_start=~ s/\s//g; $ref_end=~ s/\s//g;
	$ref_strand=~ s/\s//g; $ref_fam=~ s/\s//g;
	$ref_elem=~ s/\s//g;

	print STDERR '.';
	$j=$i;
	
	# only look for relationships with elements located before (pos strand) it and
	# ending within $MAX_RANGE bases before ref_start or anytime after it
	while(($j!=0) && ($table[$j-1][6] > $ref_start-$MAX_RANGE)) {
		$ups_ctr++;
		$j--;# for next iteration to look for element located before ref element
		$sub_start=$table[$j][5]; $sub_end=$table[$j][6];
		$sub_strand=$table[$j][8]; $sub_fam=$table[$j][9];
		$sub_elem=$table[$j][14];
		
		# cleaning up
		$sub_start=~ s/\s//g; $sub_end=~ s/\s//g;
		$sub_strand=~ s/\s//g; $sub_fam=~ s/\s//g;
		$sub_elem=~ s/\s//g;
		
		# Note: since all relationship are exclusive, I have used elsif
		
		# Spl case of U and D rel
		# Location of ref fam is entirely WITHIN sub fam
		# Should be first bcos if sub start is near the ref start, it will
		# be listed right before the ref record in the list
		if(($sub_start <= $ref_start) && ($sub_end >= $ref_end)){
			# Keeping track of bases for CONT and IN relations
			# for calculating std. dev's
			# NOTE: NOT DOING ANYTHING FOR THESE RIGHT NOW
		
			if ($ref_strand eq $sub_strand){
				if($ref_strand eq '+'){# pos strand
					# check if the ref element has this relation with this sub family already
					# check if the sub element has this relation with this ref family already
					if ((!(exists $pos_rel_history{"$ref_fam $ref_elem $sub_fam D"}) || ($pos_rel_history{"$ref_fam $ref_elem $sub_fam D"} ne $sub_fam)) && (!(exists $pos_rel_history{"$sub_fam $sub_elem $ref_fam U"}) || ($pos_rel_history{"$sub_fam $sub_elem $ref_fam U"} ne $ref_fam))) {
						# create history entry
						$pos_rel_history{"$ref_fam $ref_elem $sub_fam D"} =$sub_fam;
						$pos_rel_history{"$sub_fam $sub_elem $ref_fam U"} = $ref_fam;# changed from D to U
						
						# increment relationship count or create relationship entry
						if (!exists $pos_relationships{"$ref_fam $sub_fam D"}) {
							$pos_relationships{"$ref_fam $sub_fam D"} = 1;
						}
						else{
							$pos_relationships{"$ref_fam $sub_fam D"}++;
						}
						
						# add record for D relationship
						if($copy_file_flag){
							$pos_copies[$pos_copy_ctr][0]=$ref_fam; $pos_copies[$pos_copy_ctr][1]='D';
							$pos_copies[$pos_copy_ctr][2]=$sub_fam; $pos_copies[$pos_copy_ctr][3]='+';
							$pos_copies[$pos_copy_ctr][4]=$ref_start; $pos_copies[$pos_copy_ctr][5]=$ref_end;
							$pos_copies[$pos_copy_ctr][6]=$sub_start; $pos_copies[$pos_copy_ctr++][7]=$sub_end;
						}
						
						#NOTE: DEBUGGING
						# recording data for std dev
						# add to existing list or create relationship entry
# 						if (!exists $std_dev{"$ref_fam $sub_fam D +"}) {
# 							@temp=();
# 							$temp[0]=0;# put value in array and then assign
# 							$std_dev{"$ref_fam $sub_fam D +"}=[@temp];
# 						}
# 						else{
# 							# add in number to the array in the value list
# 							$k=$std_dev{"$ref_fam $sub_fam D +"};
# 							@temp=@$k;# necessary bcos value of hash is pointer to array
# 							$temp[$#temp+1]=0;
# 							$std_dev{"$ref_fam $sub_fam D +"}=[@temp];
# 							@temp=();#emptying array
# 						}
					}
				}
				elsif($ref_strand eq 'C'){# comp strand
					# check if the ref element has this relation with this family already
					if ((!(exists $comp_rel_history{"$ref_fam $ref_elem $sub_fam D"}) || ($comp_rel_history{"$ref_fam $ref_elem $sub_fam D"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_elem $ref_fam U"}) || ($comp_rel_history{"$sub_fam $sub_elem $ref_fam U"} ne $ref_fam))) {
						# create history entry
						$comp_rel_history{"$ref_fam $ref_elem $sub_fam D"} = $sub_fam;
						$comp_rel_history{"$sub_fam $sub_elem $ref_fam U"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $comp_relationships{"$ref_fam $sub_fam D"}) {
							$comp_relationships{"$ref_fam $sub_fam D"} = 1;
						}
						else{
							$comp_relationships{"$ref_fam $sub_fam D"}++;
						}
						
						# add record for IN relationship
						if($copy_file_flag){
							$comp_copies[$comp_copy_ctr][0]=$ref_fam; $comp_copies[$comp_copy_ctr][1]='D';
							$comp_copies[$comp_copy_ctr][2]=$sub_fam; $comp_copies[$comp_copy_ctr][3]='C';
							$comp_copies[$comp_copy_ctr][4]=$ref_start; $comp_copies[$comp_copy_ctr][5]=$ref_end;
							$comp_copies[$comp_copy_ctr][6]=$sub_start; $comp_copies[$comp_copy_ctr++][7]=$sub_end;
						}
						
						# NOTE: DEBUGGING
# 						print STDERR "\nAdded IN rel for comp strand for :$ref_fam D $sub_fam\n"; 
					}
				}
			}
		}# IN end
		
		# Upstream and Overlap: If sub fam starts before ref fam or overlap is more than 1 bp (Ovlap)
		# now if subject fam starts before the reference fam
		# and ends anytime after the start of ref fam
		# NOTE: Upstream and Ovlap are the same now
		elsif (($sub_start < $ref_start) &&  ($sub_end <= $ref_end)) {
			# Overlap case
			my ($dist);
			# will be neg value for overlaps
			$dist = $ref_start - $sub_end;
			
			# do direct rel
			if ($ref_strand eq $sub_strand){
				if($ref_strand eq '+'){# pos strand
					# check if the ref element has this relation with this family already
					if ((!(exists $pos_rel_history{"$ref_fam $ref_elem $sub_fam D"}) || ($pos_rel_history{"$ref_fam $ref_elem $sub_fam D"} ne $sub_fam)) && (!(exists $pos_rel_history{"$sub_fam $sub_elem $ref_fam U"}) || ($pos_rel_history{"$sub_fam $sub_elem $ref_fam U"} ne $ref_fam))) {
						# create history entry
						$pos_rel_history{"$ref_fam $ref_elem $sub_fam D"} = $sub_fam;
						$pos_rel_history{"$sub_fam $sub_elem $ref_fam U"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $pos_relationships{"$ref_fam $sub_fam D"}) {
							$pos_relationships{"$ref_fam $sub_fam D"} = 1;
						}
						else{
							$pos_relationships{"$ref_fam $sub_fam D"}++;
						}
						
						# add record for relationship
						if($copy_file_flag){
							$pos_copies[$pos_copy_ctr][0]=$ref_fam; $pos_copies[$pos_copy_ctr][1]='D';
							$pos_copies[$pos_copy_ctr][2]=$sub_fam; $pos_copies[$pos_copy_ctr][3]='+';
							$pos_copies[$pos_copy_ctr][4]=$ref_start; $pos_copies[$pos_copy_ctr][5]=$ref_end;
							$pos_copies[$pos_copy_ctr][6]=$sub_start; $pos_copies[$pos_copy_ctr++][7]=$sub_end;
						}
						
						# recording data for std dev
						# add to existing list or create relationship entry
						if (!exists $std_dev{"$ref_fam $sub_fam D +"}) {
							@temp=();
							$temp[0]=$dist;# put value in array and then assign
							$std_dev{"$ref_fam $sub_fam D +"}=[@temp];
						}
						else{
							# add in number to the array in the value list
							$k=$std_dev{"$ref_fam $sub_fam D +"};
							@temp=@$k;# necessary bcos value of hash is pointer to array
							$temp[$#temp+1]=$dist;
							$std_dev{"$ref_fam $sub_fam D +"}=[@temp];
							@temp=();#emptying array
						}
					}
				}
				elsif($ref_strand eq 'C'){# comp strand
					# check if the ref element has this relation with this family already
					if ((!(exists $comp_rel_history{"$ref_fam $ref_elem $sub_fam U"}) || ($comp_rel_history{"$ref_fam $ref_elem $sub_fam U"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_elem $ref_fam D"}) || ($comp_rel_history{"$sub_fam $sub_elem $ref_fam D"} ne $ref_fam))) {
						# create history entry
						$comp_rel_history{"$ref_fam $ref_elem $sub_fam U"} = $sub_fam;
						$comp_rel_history{"$sub_fam $sub_elem $ref_fam D"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $comp_relationships{"$ref_fam $sub_fam U"}) {
							$comp_relationships{"$ref_fam $sub_fam U"} = 1;
						}
						else{
							$comp_relationships{"$ref_fam $sub_fam U"}++;
						}
						
						# add record for relationship
						if($copy_file_flag){
							$comp_copies[$comp_copy_ctr][0]=$ref_fam; $comp_copies[$comp_copy_ctr][1]='U';
							$comp_copies[$comp_copy_ctr][2]=$sub_fam; $comp_copies[$comp_copy_ctr][3]='C';
							$comp_copies[$comp_copy_ctr][4]=$ref_start; $comp_copies[$comp_copy_ctr][5]=$ref_end;
							$comp_copies[$comp_copy_ctr][6]=$sub_start; $comp_copies[$comp_copy_ctr++][7]=$sub_end;
						}
						# recording data for std dev
						# add to existing list or create relationship entry
						if (!exists $std_dev{"$ref_fam $sub_fam U C"}) {
							@temp=();
							$temp[0]=$dist;# put value in array and then assign
							$std_dev{"$ref_fam $sub_fam U C"}=[@temp];
						}
						else{
							# add in number to the array in the value list
							$k=$std_dev{"$ref_fam $sub_fam U C"};
							@temp=@$k;# necessary bcos value of hash is pointer to array
							$temp[$#temp+1]=$dist;
							$std_dev{"$ref_fam $sub_fam U C"}=[@temp];
							@temp=();#emptying array
						}
					}
				}
			}
		}
	}# end while	
	
	$j=$i;
	
	# only look for relationships with elements located after it
	# and starting within $MAX_RANGE bases after ref_end (enforced by condition above)
	# or anytime after ref_start (enforced by sorting the list on start pos)
	while(($j!=$#table) && ($table[$j+1][5] < $ref_end+$MAX_RANGE)){
		$dns_ctr++;
		$j++;
		$sub_start=$table[$j][5]; $sub_end=$table[$j][6];
		$sub_strand=$table[$j][8]; $sub_fam=$table[$j][9];
		$sub_elem=$table[$j][14];
		
		# cleaning up
		$sub_start=~ s/\s//g; $sub_end=~ s/\s//g;
		$sub_strand=~ s/\s//g; $sub_fam=~ s/\s//g;
		$sub_elem=~ s/\s//g;

		# Note: since all relationship are exclusive, I have used elsif
		
		# Spl case of U and D rel
		# Location of sub fam is entirely WITHIN ref fam
		# IN should be first bcos if sub start is near the ref start, it will
		# be listed right after the ref record in the list
		if(($sub_start >= $ref_start) && ($sub_end <= $ref_end)){
			# Keeping track of bases for CONT and IN relations
			# for calculating std. dev's
			# NOTE: NOT DOING ANYTHING FOR THESE RIGHT NOW
		
			if ($ref_strand eq $sub_strand){
				if($ref_strand eq '+'){# pos strand
					# check if the ref element has this relation with this sub family already
					# check if the sub element has this relation with this ref family already
					if ((!(exists $pos_rel_history{"$ref_fam $ref_elem $sub_fam U"}) || ($pos_rel_history{"$ref_fam $ref_elem $sub_fam U"} ne $sub_fam)) && (!(exists $pos_rel_history{"$sub_fam $sub_elem $ref_fam D"}) || ($pos_rel_history{"$sub_fam $sub_elem $ref_fam D"} ne $ref_fam))) {
						# create history entry
						$pos_rel_history{"$ref_fam $ref_elem $sub_fam U"} =$sub_fam;
						$pos_rel_history{"$sub_fam $sub_elem $ref_fam D"} = $ref_fam;# changed from U to D
						
						# increment relationship count or create relationship entry
						if (!exists $pos_relationships{"$ref_fam $sub_fam U"}) {
							$pos_relationships{"$ref_fam $sub_fam U"} = 1;
						}
						else{
							$pos_relationships{"$ref_fam $sub_fam U"}++;
						}

						# add record for U relationship
						if($copy_file_flag){
							$pos_copies[$pos_copy_ctr][0]=$ref_fam; $pos_copies[$pos_copy_ctr][1]='U';
							$pos_copies[$pos_copy_ctr][2]=$sub_fam; $pos_copies[$pos_copy_ctr][3]='+';
							$pos_copies[$pos_copy_ctr][4]=$ref_start; $pos_copies[$pos_copy_ctr][5]=$ref_end;
							$pos_copies[$pos_copy_ctr][6]=$sub_start; $pos_copies[$pos_copy_ctr++][7]=$sub_end;
						}
						
						#NOTE: DEBUGGING
						# recording data for std dev
						# add to existing list or create relationship entry
# 						if (!exists $std_dev{"$ref_fam $sub_fam U +"}) {
# 							@temp=();
# 							$temp[0]=0;# put value in array and then assign
# 							$std_dev{"$ref_fam $sub_fam U +"}=[@temp];
# 						}
# 						else{
# 							# add in number to the array in the value list
# 							$k=$std_dev{"$ref_fam $sub_fam U +"};
# 							@temp=@$k;# necessary bcos value of hash is pointer to array
# 							$temp[$#temp+1]=0;
# 							$std_dev{"$ref_fam $sub_fam U +"}=[@temp];
# 							@temp=();#emptying array
# 						}
					}
				}
				elsif($ref_strand eq 'C'){# comp strand
					# check if the ref element has this relation with this family already
					if ((!(exists $comp_rel_history{"$ref_fam $ref_elem $sub_fam U"}) || ($comp_rel_history{"$ref_fam $ref_elem $sub_fam U"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_elem $ref_fam D"}) || ($comp_rel_history{"$sub_fam $sub_elem $ref_fam D"} ne $ref_fam))) {
						# create history entry
						$comp_rel_history{"$ref_fam $ref_elem $sub_fam U"} = $sub_fam;
						$comp_rel_history{"$sub_fam $sub_elem $ref_fam D"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $comp_relationships{"$ref_fam $sub_fam U"}) {
							$comp_relationships{"$ref_fam $sub_fam U"} = 1;
						}
						else{
							$comp_relationships{"$ref_fam $sub_fam U"}++;
						}
						
						# add record for IN relationship
						if($copy_file_flag){
							$comp_copies[$comp_copy_ctr][0]=$ref_fam; $comp_copies[$comp_copy_ctr][1]='U';
							$comp_copies[$comp_copy_ctr][2]=$sub_fam; $comp_copies[$comp_copy_ctr][3]='C';
							$comp_copies[$comp_copy_ctr][4]=$ref_start; $comp_copies[$comp_copy_ctr][5]=$ref_end;
							$comp_copies[$comp_copy_ctr][6]=$sub_start; $comp_copies[$comp_copy_ctr++][7]=$sub_end;
						}
						
						# NOTE: DEBUGGING
# 						print STDERR "\nAdded IN rel for comp strand for :$ref_fam U $sub_fam\n"; 
					}
				}
			}
		}# IN end
		# Downstream and Overlap: If sub fam starts after ref fam or overlap is more than 1 bp (Ovlap)
		# now if subject fam starts after the reference fam
		# and ends anytime after the end of ref fam
		# NOTE: Downstream and Ovlap are the same now
		elsif (($sub_start >= $ref_start) &&  ($sub_end >= $ref_end)) {
			# Overlap case
			my ($dist);
			# will be neg value for overlaps
			$dist = $sub_start - $ref_end;
			
			# do direct rel
			if ($ref_strand eq $sub_strand){
				if($ref_strand eq '+'){# pos strand
					# check if the ref element has this relation with this family already
					if ((!(exists $pos_rel_history{"$ref_fam $ref_elem $sub_fam U"}) || ($pos_rel_history{"$ref_fam $ref_elem $sub_fam U"} ne $sub_fam)) && (!(exists $pos_rel_history{"$sub_fam $sub_elem $ref_fam D"}) || ($pos_rel_history{"$sub_fam $sub_elem $ref_fam D"} ne $ref_fam))) {
						# create history entry
						$pos_rel_history{"$ref_fam $ref_elem $sub_fam U"} = $sub_fam;
						$pos_rel_history{"$sub_fam $sub_elem $ref_fam D"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $pos_relationships{"$ref_fam $sub_fam U"}) {
							$pos_relationships{"$ref_fam $sub_fam U"} = 1;
						}
						else{
							$pos_relationships{"$ref_fam $sub_fam U"}++;
						}
						
						# add record for relationship
						if($copy_file_flag){
							$pos_copies[$pos_copy_ctr][0]=$ref_fam; $pos_copies[$pos_copy_ctr][1]='U';
							$pos_copies[$pos_copy_ctr][2]=$sub_fam; $pos_copies[$pos_copy_ctr][3]='+';
							$pos_copies[$pos_copy_ctr][4]=$ref_start; $pos_copies[$pos_copy_ctr][5]=$ref_end;
							$pos_copies[$pos_copy_ctr][6]=$sub_start; $pos_copies[$pos_copy_ctr++][7]=$sub_end;
						}
						
						# recording data for std dev
						# add to existing list or create relationship entry
						if (!exists $std_dev{"$ref_fam $sub_fam U +"}) {
							@temp=();
							$temp[0]=$dist;# put value in array and then assign
							$std_dev{"$ref_fam $sub_fam U +"}=[@temp];
						}
						else{
							# add in number to the array in the value list
							$k=$std_dev{"$ref_fam $sub_fam U +"};
							@temp=@$k;# necessary bcos value of hash is pointer to array
							$temp[$#temp+1]=$dist;
							$std_dev{"$ref_fam $sub_fam U +"}=[@temp];
							@temp=();#emptying array
						}
					}
				}
				elsif($ref_strand eq 'C'){# comp strand
					# check if the ref element has this relation with this family already
					if ((!(exists $comp_rel_history{"$ref_fam $ref_elem $sub_fam D"}) || ($comp_rel_history{"$ref_fam $ref_elem $sub_fam D"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_elem $ref_fam U"}) || ($comp_rel_history{"$sub_fam $sub_elem $ref_fam U"} ne $ref_fam))) {
						# create history entry
						$comp_rel_history{"$ref_fam $ref_elem $sub_fam D"} = $sub_fam;
						$comp_rel_history{"$sub_fam $sub_elem $ref_fam U"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						if (!exists $comp_relationships{"$ref_fam $sub_fam D"}) {
							$comp_relationships{"$ref_fam $sub_fam D"} = 1;
						}
						else{
							$comp_relationships{"$ref_fam $sub_fam D"}++;
						}
						
						# add record for relationship
						if($copy_file_flag){
							$comp_copies[$comp_copy_ctr][0]=$ref_fam; $comp_copies[$comp_copy_ctr][1]='D';
							$comp_copies[$comp_copy_ctr][2]=$sub_fam; $comp_copies[$comp_copy_ctr][3]='C';
							$comp_copies[$comp_copy_ctr][4]=$ref_start; $comp_copies[$comp_copy_ctr][5]=$ref_end;
							$comp_copies[$comp_copy_ctr][6]=$sub_start; $comp_copies[$comp_copy_ctr++][7]=$sub_end;
						}
						# recording data for std dev
						# add to existing list or create relationship entry
						if (!exists $std_dev{"$ref_fam $sub_fam D C"}) {
							@temp=();
							$temp[0]=$dist;# put value in array and then assign
							$std_dev{"$ref_fam $sub_fam D C"}=[@temp];
						}
						else{
							# add in number to the array in the value list
							$k=$std_dev{"$ref_fam $sub_fam D C"};
							@temp=@$k;# necessary bcos value of hash is pointer to array
							$temp[$#temp+1]=$dist;
							$std_dev{"$ref_fam $sub_fam D C"}=[@temp];
							@temp=();#emptying array
						}
					}
				}
			}
		}

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


# UPDATING THE F_ITEMSETS
# RECIPROCAL COUNTS, IF ANY FOUND ON THE OPPOSITE STRAND
# A u1 B on pos, B d1 A on neg

# %pos_relationships : [fam1 fam2 category] = count
# %comp_relationships : [fam1 fam2 category] = count
# %both_relationships : [fam1 fam2 category] = count
# %std_dev : ["fam1 fam2 rel strand"] => [array of deviations]

# updating relationships on the positive strand
# saving orig %pos_relationships for processing %comp_relationships
# so that the orig counts can be used instead of the consolidated ones
%temphash=%pos_relationships;
while( ($i,$j) = each %pos_relationships){
	
	@temp=split(' ',$i);
	
	# only for U rels
	if($temp[2] eq 'U'){
		# check if dev vlues have been recorded
		# since this can be a IN only relationship with 
		# no recorded deviations
		if (exists $std_dev{"$temp[0] $temp[1] U +"}){
			@old_devs=@devs=();#init arrays
			
			#add the deviations to the std_dev hash
			$k=$std_dev{"$temp[0] $temp[1] U +"};# get old devs
			@old_devs=@$k;# necessary bcos value of hash is pointer to array
			
			#NOTE: DEBUGGING
# 			if($j != ($#old_devs+1)) { 
# 				print STDERR "Dev count mismatch error $temp[0] U $temp[1] + , count:$j, nof devs:",$#old_devs+1,"\n";
# 			}
		
			# if a reci rel exists for this family pair
			if (exists $comp_relationships{"$temp[0] $temp[1] U"}) {
				$j+=$comp_relationships{"$temp[0] $temp[1] U"};
				
				# increment the current count
				$pos_relationships{$i}=$j;
				
				# check if dev vlues have been recorded
				# since this can be a IN only relationship with 
				# no recorded deviations
				if (exists $std_dev{"$temp[0] $temp[1] U C"}){
					$k=$std_dev{"$temp[0] $temp[1] U C"};# get new devs
					@devs=@$k;# necessary bcos value of hash is pointer to array
					
					# add new devs to old devs
					for($k=0;$k<(scalar(@devs)-1);$k++){
						$old_devs[$#old_devs+1]=$devs[$k];
					}
					#assign updated dev array
					$std_dev{"$temp[0] $temp[1] U +"}=[@old_devs];
				}
			}
		}
	}
	# only for D rels
	elsif($temp[2] eq 'D'){
		# check if dev vlues have been recorded
		# since this can be a IN only relationship with 
		# no recorded deviations
		if (exists $std_dev{"$temp[0] $temp[1] D +"}){
			@old_devs=@devs=();#init arrays
			
			#add the deviations to the std_dev hash
			$k=$std_dev{"$temp[0] $temp[1] D +"};# get old devs
			@old_devs=@$k;# necessary bcos value of hash is pointer to array
			
			#NOTE: DEBUGGING
# 			if($j != ($#old_devs+1)) { 
# 				print STDERR "Dev count mismatch error, $temp[0] D $temp[1] + , count:$j, nof devs:",$#old_devs+1,"\n";
# 			}
		
			# if a reci rel exists for this family pair
			if (exists $comp_relationships{"$temp[0] $temp[1] D"}) {
				$j+=$comp_relationships{"$temp[0] $temp[1] D"};
				
				# increment the current count
				$pos_relationships{$i}=$j;
				
				if (exists $std_dev{"$temp[0] $temp[1] D C"}){
					$k=$std_dev{"$temp[0] $temp[1] D C"};# get new devs
					@devs=@$k;# necessary bcos value of hash is pointer to array
					
					# add new devs to old devs
					for($k=0;$k<(scalar(@devs)-1);$k++){
						$old_devs[$#old_devs+1]=$devs[$k];
					}
					
					#assign updated dev array
					$std_dev{"$temp[0] $temp[1] D +"}=[@old_devs];
				}
			}
		}
	}
}

# updating relationships on the comp strand
while( ($i,$j) = each %comp_relationships){
	
	@temp=split(' ',$i);

	# only for U rels
	if($temp[2] eq 'U'){
		# check if dev vlues have been recorded
		# since this can be a IN only relationship with 
		# no recorded deviations
		if (exists $std_dev{"$temp[0] $temp[1] U C"}){
			@old_devs=@devs=();#init arrays
			
			#add the deviations to the std_dev hash
			$k=$std_dev{"$temp[0] $temp[1] U C"};# get old devs
			@old_devs=@$k;# necessary bcos value of hash is pointer to array
			#NOTE: DEBUGGING
# 			if($j != ($#old_devs+1)) { print STDERR "Dev count mismatch error on Comp\n";}

			# if a reci rel exists for this family pair
			if (exists $temphash{"$temp[0] $temp[1] U"}) {
				$j+=$temphash{"$temp[0] $temp[1] U"};
				
				# increment the current count
				$comp_relationships{$i}=$j;
				
				# check if dev vlues have been recorded
				# since this can be a IN only relationship with 
				# no recorded deviations
				if (exists $std_dev{"$temp[0] $temp[1] U +"}){
					$k=$std_dev{"$temp[0] $temp[1] U +"};# get new devs
					@devs=@$k;# necessary bcos value of hash is pointer to array
					
					# add new devs to old devs
					for($k=0;$k<(scalar(@devs)-1);$k++){
						$old_devs[$#old_devs+1]=$devs[$k];
					}
					#assign updated dev array
					$std_dev{"$temp[0] $temp[1] U C"}=[@old_devs];
				}
			}
		}
	}
	# only for D rels
	elsif($temp[2] eq 'D'){
		# check if dev vlues have been recorded
		# since this can be a IN only relationship with 
		# no recorded deviations
		if (exists $std_dev{"$temp[0] $temp[1] D C"}){
			@old_devs=@devs=();#init arrays
			
			#add the deviations to the std_dev hash
			$k=$std_dev{"$temp[0] $temp[1] D C"};# get old devs
			@old_devs=@$k;# necessary bcos value of hash is pointer to array
			
			#NOTE: DEBUGGING
# 			if($j != ($#old_devs+1)) { print STDERR "Dev count mismatch error on Comp\n";}
		
			# if a reci rel exists for this family pair
			if (exists $temphash{"$temp[0] $temp[1] D"}) {
				$j+=$temphash{"$temp[0] $temp[1] D"};
				
				# increment the current count
				$comp_relationships{$i}=$j;
				
				# check if dev vlues have been recorded
				# since this can be a IN only relationship with 
				# no recorded deviations
				if (exists $std_dev{"$temp[0] $temp[1] D +"}){
					$k=$std_dev{"$temp[0] $temp[1] D +"};# get new devs
					@devs=@$k;# necessary bcos value of hash is pointer to array
					
					# add new devs to old devs
					for($k=0;$k<(scalar(@devs)-1);$k++){
						$old_devs[$#old_devs+1]=$devs[$k];
					}
					# assign updated dev array
					$std_dev{"$temp[0] $temp[1] D C"}=[@old_devs];
				}
			}
		}
	}
}

#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\n\# Runtime details after updating relationships: \n";
print OUTFILEDATA "\# System time for process: ",ceil($system_t/60)," mins\n";
print OUTFILEDATA "\# User time for process: ",ceil($user_t/60)," mins\n\n\n";

print STDERR "\n\# Runtime details after updating relationships: \n";
print STDERR "\# System time for process: ",ceil($system_t/60)," mins\n";
print STDERR "\# User time for process: ",ceil($user_t/60)," mins\n\n";


# PRINTING THE ITEMSETS
# %pos_relationships : [fam1 fam2 category] = count
# %comp_relationships : [fam1 fam2 category] = count
# %both_relationships : [fam1 fam2 category] = count
# @count: fam occurences avg-len elementnum
# %std_dev : ["fam1 fam2 rel strand"] => [array of deviations]

# Creating the frequent itemsets in the format
# fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, Occurence, Strand, Category, Avg dist, Std dev
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
else{print OUTFILEDATA "\n\n";}


print OUTFILEDATA "\n\# Frequent itemsets from positive strand: \n";

# printing relationships on the positive strand
while( ($i,$j) = each %pos_relationships){
	@temp=split(' ',$i);
	$rec=&get_index($temp[0]);
	print OUTFILEDATA "$temp[0]\t$counts[$rec][1]\t$counts[$rec][2]\t";
	$rec=&get_index($temp[1]);
	print OUTFILEDATA "$temp[1]\t$counts[$rec][1]\t$counts[$rec][2]\t";
	print OUTFILEDATA "$j\t+\t$temp[2]\t";
	
	#printing out the avg dist and std dev
	# check if dev vlues have been recorded
	# since this can be a IN only relationship with 
	# no recorded deviations
	if (exists $std_dev{"$temp[0] $temp[1] $temp[2] +"}){
		@devs=();
		$k=$std_dev{"$temp[0] $temp[1] $temp[2] +"};
		@devs=@$k;# necessary bcos value of hash is pointer to array
		$avg=0;
		
		# get the sum
		foreach (@devs){
			$avg+=$_;
		}
		$avg=int($avg/($#devs+1));
		print OUTFILEDATA "$avg\t";# printing average
		
		$j=0;# reusing $j
		foreach (@devs){
			# find the dev and square it
			$j+= ($_ - $avg) * ($_ - $avg);
		}
		
		$j=int($j/($#devs+1));
		
		$j = sqrt ($j);
		print OUTFILEDATA &round2($j),"\n";# printing std dev
	}
	else{
		print OUTFILEDATA "NA\tNA\n";# printing std dev for IN type rels
	}
}
#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\n\# Runtime details after printing positive itemsets: \n";
print OUTFILEDATA "\# System time for process: ",ceil($system_t/60)," mins\n";
print OUTFILEDATA "\# User time for process: ",ceil($user_t/60)," mins\n\n";

print STDERR "\n\# Runtime details after printing positive itemsets: \n";
print STDERR "\# System time for process: ",ceil($system_t/60)," mins\n";
print STDERR "\# User time for process: ",ceil($user_t/60)," mins\n";

print OUTFILEDATA "\n\# Frequent itemsets from complimentary strand: \n";

# printing relationships on the comp strand
while( ($i,$j) = each %comp_relationships){
	@temp=split(' ',$i);
	$rec=&get_index($temp[0]);
	print OUTFILEDATA "$temp[0]\t$counts[$rec][1]\t$counts[$rec][2]\t";
	$rec=&get_index($temp[1]);
	print OUTFILEDATA "$temp[1]\t$counts[$rec][1]\t$counts[$rec][2]\t";
	print OUTFILEDATA "$j\tC\t$temp[2]\t";
	
	#printing out the avg dist and std dev
	# check if dev vlues have been recorded
	# since this can be a IN only relationship with 
	# no recorded deviations
	if (exists $std_dev{"$temp[0] $temp[1] $temp[2] C"}){
		@devs=();
		$k=$std_dev{"$temp[0] $temp[1] $temp[2] C"};
		@devs=@$k;# necessary bcos value of hash is pointer to array
		$avg=0;
		
		# get the sum
		foreach (@devs){
			$avg+=$_;
		}
		$avg=int($avg/($#devs+1));
		print OUTFILEDATA "$avg\t";# printing average
		
		$j=0;# reusing $j
		foreach (@devs){
			# find the dev and square it
			$j+= ($_ - $avg) * ($_ - $avg);
		}
		
		$j=int($j/($#devs+1));
		
		$j = sqrt ($j);
		print OUTFILEDATA &round2($j),"\n";# printing std dev
	}
	else{
		print OUTFILEDATA "NA\tNA\n";# printing std dev for IN type rels
	}
}
#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\n\# Runtime details after printing negative itemsets: \n";
print OUTFILEDATA "\# System time for process: ",ceil($system_t/60)," mins\n";
print OUTFILEDATA "\# User time for process: ",ceil($user_t/60)," mins\n";

print STDERR "\n\# Runtime details after printing negative itemsets: \n";
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
	
	#calculating time taken
	($user_t,$system_t,$cuser_t,$csystem_t) = times;
	print OUTFILECOPIESPOS "\n\# Runtime details after printing copy info: \n";
	print OUTFILECOPIESPOS "\# System time for process: ",ceil($system_t/60)," mins\n";
	print OUTFILECOPIESPOS "\# User time for process: ",ceil($user_t/60)," mins\n";
	print OUTFILECOPIESCOMP "\n\# Runtime details after printing copy info: \n";
	print OUTFILECOPIESCOMP "\# System time for process: ",ceil($system_t/60)," mins\n";
	print OUTFILECOPIESCOMP "\# User time for process: ",ceil($user_t/60)," mins\n";

	
	print STDERR "\n\# Runtime details after printing copy info: \n";
	print STDERR "\# System time for process: ",ceil($system_t/60)," mins\n";
	print STDERR "\# User time for process: ",ceil($user_t/60)," mins\n";
}

close (OUTFILEDATA);
if($copy_file_flag){
	close (OUTFILECOPIESPOS);
	close (OUTFILECOPIESCOMP);
}
# # debugging
# close (ERRFILE);
exit;
