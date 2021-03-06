#!/usr/bin/perl -w
# MGEL
# Surya Saha 02/02/08
# reading cmd line input .gff file which is sorted on the start position
# and finds the relationship among elements and families
# Relationship types: 
# Upstream: 
#	u (0-MAX bases)
#	Overlap
#	spl case Contains: Location of fam2 is entirely within fam1 (Cont) ??
# Downstream: 
#	d (0-MAX bases)
#	Overlap
#	spl case In: Location of fam2 is entirely within fam1 (IN) ??


# OUTPUT: Creating the frequent itemsets in the format
# fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, Occurence, Strand, Category, Avg. dist., Std. dev.
# Printing "NA" for Avg. dist., Std. dev. in case of IN/CONT relations

# v1: Using only 1 MAX range for U and D relation. No other kind of relation.
#   : No config file required. 
#   : No output for both strands (irrespective of strand)
#   : Added the avg and std. dev. for each family pair in a relationship
#   : Not accounting for IN/CONT type rels in std devs. Just putting 0 in list of distances.

# v2: Just putting 0 in list of distances for IN rels.
#   : Using only upstream relationships for both strands since we are scanning records in the 
#   : downstream (positive) direction. See description for details.
#   : Only 1 output file for combined pos and comp strand information.
#   : Only 1 copy file for all elements involved in relationships from both strands.
# NOTE: we can optimize this tool by looking only on the right for elelmets to form relationships with

# v3: fixed a bug (negative deviations were messing up the calculation of std dev)

# v1: Reads GFF files
# v2: Sorts the input GFF file on the starting position

# GFF v 2.0 fields are: 
# <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]

#226 29.3  0.0  0.0 chr12|12012              10808     10865 27486349 +  R=112           Unknown               1      58     (0)
#241 26.3  0.0  0.0 chr12|12012              10836     10892 27486322 C  R=112           Unknown             (0)      58       2

#   1935 10.6  0.0  2.8 chr12|12012               8936      9225 27487989 C  R=291           Unknown             (0)     283       2
#   1009 22.1  2.2  4.2 chr12|12012               9987     10255 27486959 +  R=39            Unknown               1     264     (0)


my $ver=1;

use strict;
use warnings;
use POSIX;

my  ($ifname,$rec,@temp,@temp1,$ctr,$i,$j,$k,$copy_file_flag,$MAX_RANGE);

my (@table,@famnames,@counts,%counts_index,$ups_ctr,$dns_ctr, $ref_elem,$ref_fam,$ref_start,$ref_end,
$ref_strand,$sub_fam,$sub_elem,$sub_start,$sub_end,$sub_strand,%pos_relationships,%comp_relationships,
%relationships,$pos_copy_ctr,$comp_copy_ctr,@pos_copies,@comp_copies,%pos_rel_history,%comp_rel_history,
$tot_fams,$tot_recs,$user_t,$system_t,$cuser_t,$csystem_t);

my (%reci_rels, $rel, @devs, @old_devs, %temphash, %std_dev, $avg);

# Supporting routines

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


# Command line arguments

unless (@ARGV == 3){
	print "USAGE: $0 <MAX_RANGE> <.gff file> <-nocop or -cop>\n";
	exit;
}

$MAX_RANGE=$ARGV[0];
chomp $MAX_RANGE;

$ifname=$ARGV[1];
chomp $ifname;
unless(open(INFILEDATA,$ifname)){print "not able to open ".$ifname."\n\n";exit;}
unless(open(OUTFILEDATA,">$ifname.f_itemsets.stg1.tab")){print "not able to open ".$ifname."f_itemsets.stg1.tab\n\n";exit;}

$copy_file_flag=$ARGV[2];
chomp $copy_file_flag;
if(!($copy_file_flag eq "-nocop" || $copy_file_flag eq "-cop")){ print STDERR "flag can be only -cop or -nocop\nExiting..\n"; exit;}
# convert to 0 or 1
if($copy_file_flag eq "-nocop"){ $copy_file_flag=0;}
elsif ($copy_file_flag eq "-cop"){ $copy_file_flag=1;}

if($copy_file_flag){
	unless(open(OUTFILECOPIES,">$ifname.copies.tab")){print "not able to open ".$ifname.".copies.tab \n\n";exit;}
}


# <BODY>

# SLURPING IN THE WHOLE .OUT REPORT FILE
# <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
#  0	           1       2         3      4      5       6         7       8          9 
# <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
# chr1	RptSct-1.0.1	repeat_unit	9599047	9599208	606	+	.	R=825
# chr1	RptSct-1.0.1	repeat_unit	20233020	20233420	2669	+	.	R=825

#     0    1    2    3      4             5        6       7    8     9         10      11   12    13
#   1935 10.6  0.0  2.8 chr12|12012     8936      9225 27487989 C  R=291     Unknown    (0)  283   2
#   1009 22.1  2.2  4.2 chr12|12012     9987     10255 27486959 +  R=39      Unknown     1   264  (0)


$ctr=0;
while($rec=<INFILEDATA>){
	if($rec =~ /#/){next;}
	if(length ($rec) < 10){next;}#for avoiding last line
	@temp=split(' ',$rec);
	
	# putting into .out file format since rest of hte code was built for .out files
# 	@temp1={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	$temp1[5]=$temp[3];#start
	$temp1[6]=$temp[4];#end
	$temp1[8]=$temp[6];#strand
	$temp1[9]=$temp[8];#famname
	
	push @table, [@temp1];
	
	@temp=();
	@temp1=();
	$ctr++;
}
# record tot recs
$tot_recs = $ctr;
close (INFILEDATA);

# sort both tables on starting location
# required in casees where our GFF file has been made by combining 2 diff GFF
# files, one from CD and one from repeats
@temp = sort {$a->[5] <=> $b->[5]} @table;
@table=@temp;


print OUTFILEDATA "\# Version: $ver\n";
print OUTFILEDATA "\# Input GFF file: $ifname\n";
if($copy_file_flag){
	print OUTFILECOPIES "\# Version: $ver\n";
	print OUTFILECOPIES "\# Input GFF file: $ifname\n";
}
$i=localtime();
print OUTFILEDATA "\# Time: $i\n";
if($copy_file_flag){
	print OUTFILECOPIES "\# Time: $i\n";
}

#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\n\# Runtime details after reading in/sorting the file: \n";
print OUTFILEDATA "\# System time for process: ",ceil($system_t/60)," mins\n";
print OUTFILEDATA "\# User time for process: ",ceil($user_t/60)," mins\n";
print OUTFILEDATA "\n";


print STDERR "\n\# Runtime details after reading $tot_recs records from file: \n";
print STDERR "\# System time for process: ",ceil($system_t/60)," mins\n";
print STDERR "\# User time for process: ",ceil($user_t/60)," mins\n";


# PREPARING DATA FOR RULE MINING
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
# finding the relationships (category=rel)
# %pos_relationships : [fam1 fam2 category] = count
# %comp_relationships : [fam1 fam2 category] = count
# %pos_rel_history : [fam1 fam1-img fam2 category] = fam2
# %comp_rel_history : [fam1 fam1-img fam2 category] = fam2
# Feb 1, 08 : I do not remember why we have 
# %pos_rel_history : [fam1 fam1-img fam2 category] = fam2 in the code???
# but testing shows that only [fam1 fam1-img category] produces conf > 1.0
# %std_dev : ["fam1 fam2 rel strand"] => [array of deviations]

$ups_ctr=$dns_ctr=0;
$pos_copy_ctr=$comp_copy_ctr=0;

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
						
						# recording data for std dev
						# putting 0 since no valid distance value
						# add to existing list or create relationship entry
						if (!exists $std_dev{"$ref_fam $sub_fam D +"}) {
							@temp=();
							$temp[0]=0;# put value in array and then assign
							$std_dev{"$ref_fam $sub_fam D +"}=[@temp];
						}
						else{
							# add in number to the array in the value list
							$k=$std_dev{"$ref_fam $sub_fam D +"};
							@temp=@$k;# necessary bcos value of hash is pointer to array
							$temp[$#temp+1]=0;
							$std_dev{"$ref_fam $sub_fam D +"}=[@temp];
							@temp=();#emptying array
						}
					}
				}
				# elsif($ref_strand eq 'C'){# comp strand
				elsif($ref_strand eq '-'){# comp strand, changed for GFF format
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
						
						# recording data for std dev
						# putting 0 since no valid distance value
						# add to existing list or create relationship entry
						if (!exists $std_dev{"$ref_fam $sub_fam D C"}) {
							@temp=();
							$temp[0]=0;# put value in array and then assign
							$std_dev{"$ref_fam $sub_fam D C"}=[@temp];
						}
						else{
							# add in number to the array in the value list
							$k=$std_dev{"$ref_fam $sub_fam D C"};
							@temp=@$k;# necessary bcos value of hash is pointer to array
							$temp[$#temp+1]=0;
							$std_dev{"$ref_fam $sub_fam D C"}=[@temp];
							@temp=();#emptying array
						}
					}
				}
			}
		}# IN end
		
		# Left and Overlap: If sub fam starts before ref fam or overlap is more than 1 bp (Ovlap)
		# now if subject fam starts before the reference fam
		# and ends anytime after the start of ref fam
		# NOTE: Left and Ovlap are the same now
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
				#elsif($ref_strand eq 'C'){# comp strand
				elsif($ref_strand eq '-'){# comp strand, changed for GFF format
					# check if the ref element has this relation with this family already
					if ((!(exists $comp_rel_history{"$ref_fam $ref_elem $sub_fam U"}) || ($comp_rel_history{"$ref_fam $ref_elem $sub_fam U"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_elem $ref_fam D"}) || ($comp_rel_history{"$sub_fam $sub_elem $ref_fam D"} ne $ref_fam))) {
						# create history entry
						$comp_rel_history{"$ref_fam $ref_elem $sub_fam U"} = $sub_fam;
						$comp_rel_history{"$sub_fam $sub_elem $ref_fam D"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						# Ver 2 CHANGE
						if (!exists $comp_relationships{"$sub_fam $ref_fam D"}) {
							$comp_relationships{"$sub_fam $ref_fam D"} = 1;
						}
						else{
							$comp_relationships{"$sub_fam $ref_fam D"}++;
						}
						
						# add record for relationship
						if($copy_file_flag){
							$comp_copies[$comp_copy_ctr][0]=$sub_fam; $comp_copies[$comp_copy_ctr][1]='D';
							$comp_copies[$comp_copy_ctr][2]=$ref_fam; $comp_copies[$comp_copy_ctr][3]='C';
							$comp_copies[$comp_copy_ctr][4]=$sub_start; $comp_copies[$comp_copy_ctr][5]=$sub_end;
							$comp_copies[$comp_copy_ctr][6]=$ref_start; $comp_copies[$comp_copy_ctr++][7]=$ref_end;
						}
						# recording data for std dev
						# add to existing list or create relationship entry
						if (!exists $std_dev{"$sub_fam $ref_fam D C"}) {
							@temp=();
							$temp[0]=$dist;# put value in array and then assign
							$std_dev{"$sub_fam $ref_fam D C"}=[@temp];
						}
						else{
							# add in number to the array in the value list
							$k=$std_dev{"$sub_fam $ref_fam D C"};
							@temp=@$k;# necessary bcos value of hash is pointer to array
							$temp[$#temp+1]=$dist;
							$std_dev{"$sub_fam $ref_fam D C"}=[@temp];
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
						
						# recording data for std dev
						# putting 0 since no valid distance value
						# add to existing list or create relationship entry
						if (!exists $std_dev{"$ref_fam $sub_fam U +"}) {
							@temp=();
							$temp[0]=0;# put value in array and then assign
							$std_dev{"$ref_fam $sub_fam U +"}=[@temp];
						}
						else{
							# add in number to the array in the value list
							$k=$std_dev{"$ref_fam $sub_fam U +"};
							@temp=@$k;# necessary bcos value of hash is pointer to array
							$temp[$#temp+1]=0;
							$std_dev{"$ref_fam $sub_fam U +"}=[@temp];
							@temp=();#emptying array
						}
					}
				}
				#elsif($ref_strand eq 'C'){# comp strand
				elsif($ref_strand eq '-'){# comp strand, changed for GFF format
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
						
						# recording data for std dev
						# putting 0 since no valid distance value
						# add to existing list or create relationship entry
						if (!exists $std_dev{"$ref_fam $sub_fam U C"}) {
							@temp=();
							$temp[0]=0;# put value in array and then assign
							$std_dev{"$ref_fam $sub_fam U C"}=[@temp];
						}
						else{
							# add in number to the array in the value list
							$k=$std_dev{"$ref_fam $sub_fam U C"};
							@temp=@$k;# necessary bcos value of hash is pointer to array
							$temp[$#temp+1]=0;
							$std_dev{"$ref_fam $sub_fam U C"}=[@temp];
							@temp=();#emptying array
						}
					}
				}
			}
		}# IN end
		# Right and Overlap: If sub fam starts after ref fam or overlap is more than 1 bp (Ovlap)
		# now if subject fam starts after the reference fam
		# and ends anytime after the end of ref fam
		# NOTE: Right and Ovlap are the same now
		elsif (($sub_start >= $ref_start) &&  ($sub_end >= $ref_end)) {
			my ($dist);
			# will be neg value for overlaps
			$dist = $sub_start - $ref_end;
			
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
				#elsif($ref_strand eq 'C'){# comp strand
				elsif($ref_strand eq '-'){# comp strand, changed for GFF format
					# check if the ref element has this relation with this family already
					if ((!(exists $comp_rel_history{"$ref_fam $ref_elem $sub_fam D"}) || ($comp_rel_history{"$ref_fam $ref_elem $sub_fam D"} ne $sub_fam)) && (!(exists $comp_rel_history{"$sub_fam $sub_elem $ref_fam U"}) || ($comp_rel_history{"$sub_fam $sub_elem $ref_fam U"} ne $ref_fam))) {
						# create history entry
						$comp_rel_history{"$ref_fam $ref_elem $sub_fam D"} = $sub_fam;
						$comp_rel_history{"$sub_fam $sub_elem $ref_fam U"} = $ref_fam;
						
						# increment relationship count or create relationship entry
						# Ver 2 CHANGE
						if (!exists $comp_relationships{"$sub_fam $ref_fam U"}) {
							$comp_relationships{"$sub_fam $ref_fam U"} = 1;
						}
						else{
							$comp_relationships{"$sub_fam $ref_fam U"}++;
						}
						
						# add record for relationship
						if($copy_file_flag){
							$comp_copies[$comp_copy_ctr][0]=$sub_fam; $comp_copies[$comp_copy_ctr][1]='U';
							$comp_copies[$comp_copy_ctr][2]=$ref_fam; $comp_copies[$comp_copy_ctr][3]='C';
							$comp_copies[$comp_copy_ctr][4]=$sub_start; $comp_copies[$comp_copy_ctr][5]=$sub_end;
							$comp_copies[$comp_copy_ctr][6]=$ref_start; $comp_copies[$comp_copy_ctr++][7]=$ref_end;
						}
						# recording data for std dev
						# add to existing list or create relationship entry
						if (!exists $std_dev{"$sub_fam $ref_fam U C"}) {
							@temp=();
							$temp[0]=$dist;# put value in array and then assign
							$std_dev{"$sub_fam $ref_fam U C"}=[@temp];
						}
						else{
							# add in number to the array in the value list
							$k=$std_dev{"$sub_fam $ref_fam U C"};
							@temp=@$k;# necessary bcos value of hash is pointer to array
							$temp[$#temp+1]=$dist;
							$std_dev{"$sub_fam $ref_fam U C"}=[@temp];
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
# %relationships : [fam1 fam2 category] = count
# %std_dev : [fam1 fam2 rel strand] => [array of deviations]

# updating relationships on the positive strand
# only updating for positive strand as all differing relationships
# will be collected into %relationships
while( ($i,$j) = each %pos_relationships){
	
	@temp=split(' ',$i);
	
	# only for U rels
	if($temp[2] eq 'U'){
		# check if dev vlues have been recorded
		# since this can be a IN only relationship with 
		# no recorded deviations
		@old_devs=@devs=();#init arrays
		
		#add the deviations to the std_dev hash
		$k=$std_dev{"$temp[0] $temp[1] U +"};# get old devs
		@old_devs=@$k;# necessary bcos value of hash is pointer to array
		
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
	# only for D rels
	elsif($temp[2] eq 'D'){
		# check if dev vlues have been recorded
		# since this can be a IN only relationship with 
		# no recorded deviations
		@old_devs=@devs=();#init arrays
		
		#add the deviations to the std_dev hash
		$k=$std_dev{"$temp[0] $temp[1] D +"};# get old devs
		@old_devs=@$k;# necessary bcos value of hash is pointer to array
		
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


# Combine relationships from both strands
# Std_dev for all rels should still be in the hash
# Start with all relationships on pos strand
%relationships=%pos_relationships;
# Now add relationships from comp strand
while( ($i,$j) = each %comp_relationships){
	if (!(exists $relationships{"$i"})){ # new entry
		$relationships{"$i"} = $j;
	}
}

#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\n\# Runtime details after updating and combining relationships: \n";
print OUTFILEDATA "\# System time for process: ",ceil($system_t/60)," mins\n";
print OUTFILEDATA "\# User time for process: ",ceil($user_t/60)," mins\n\n\n";

print STDERR "\n\# Runtime details after updating and combining relationships: \n";
print STDERR "\# System time for process: ",ceil($system_t/60)," mins\n";
print STDERR "\# User time for process: ",ceil($user_t/60)," mins\n\n";

# PRINTING THE ITEMSETS (category=rels)
# %pos_relationships : [fam1 fam2 category] = count
# %comp_relationships : [fam1 fam2 category] = count
# %both_relationships : [fam1 fam2 category] = count
# @count: fam occurences avg-len elementnum
# %std_dev : ["fam1 fam2 rel strand"] => [array of deviations]

# Creating the frequent itemsets in the format
# fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, Occurence, Category, Avg dist, Std dev
print OUTFILEDATA "\# Total records in GFF file: $tot_recs\n";
print OUTFILEDATA "\# Total number of families: $tot_fams\n\n";
print OUTFILEDATA "\# Note: If dns rec ~ ups rec, then the regions were located uniformly\n";
print OUTFILEDATA "\# Average number of upstream GFF records processed per element: ".ceil($ups_ctr/$tot_recs)."\n";
print OUTFILEDATA "\# Average number of downstream GFF records processed per element: ".ceil($dns_ctr/$tot_recs)."\n";
print OUTFILEDATA "\# Average number of GFF records processed per element: ".ceil(($ups_ctr+$dns_ctr)/$tot_recs)."\n\n";
print OUTFILEDATA "\# Total relationships on pos strand:".keys(%pos_relationships)."\n";
if($copy_file_flag){ print OUTFILECOPIES "\# Total copies/clusters on pos strand:".$pos_copy_ctr."\n";}
print OUTFILEDATA "\# Total relationships on comp strand:".keys(%comp_relationships)."\n";
if($copy_file_flag){ print OUTFILECOPIES "\# Total copies/clusters on comp strand:".$comp_copy_ctr."\n";}
else{print OUTFILEDATA "\n\n";}

print OUTFILEDATA "\# F1\tCnt\tLen\tF2\tCnt\tLen\tRcnt\tRel\tAvgDst\tStdev\n";

while( ($i,$j) = each %relationships){
	@temp=split(' ',$i);
	$rec=&get_index($temp[0]);
	print OUTFILEDATA "$temp[0]\t$counts[$rec][1]\t$counts[$rec][2]\t";
	$rec=&get_index($temp[1]);
	print OUTFILEDATA "$temp[1]\t$counts[$rec][1]\t$counts[$rec][2]\t";
	print OUTFILEDATA "$j\t$temp[2]\t";
	
	#printing out the avg dist and std dev
	# check for both strands since all relationships
	# have been combined
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
		
		# DEBUGGING
# 		print "fam1:$temp[0] fam2:$temp[1] Count:$j Devs:",$#devs+1," Values:";
# 		foreach(@devs) {print abs ($_)," ";}
		
		$j=0;# reusing $j
		foreach (@devs){
			# find the dev and square it
			# $j+= ($_ - $avg) * ($_ - $avg);
			# to prevent negative values from messing up the std dev
			$j+= (abs ($_) - $avg) * (abs ($_) - $avg);
		}
		
		$j=int($j/($#devs+1));
		
		$j = sqrt ($j);
		print OUTFILEDATA &round2($j),"\n";# printing std dev
		
		# DEBUGGING
# 		print "   Avg:$avg StdDev:$j\n";

	}
	elsif (exists $std_dev{"$temp[0] $temp[1] $temp[2] C"}){
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
			# $j+= ($_ - $avg) * ($_ - $avg);
			# to prevent negative values from messing up the std dev
			$j+= (abs ($_) - $avg) * (abs ($_) - $avg);
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
print OUTFILEDATA "\n\# Runtime details after printing f_itemsets: \n";
print OUTFILEDATA "\# System time for process: ",ceil($system_t/60)," mins\n";
print OUTFILEDATA "\# User time for process: ",ceil($user_t/60)," mins\n\n";

print STDERR "\n\# Runtime details after printing f_itemsets: \n";
print STDERR "\# System time for process: ",ceil($system_t/60)," mins\n";
print STDERR "\# User time for process: ",ceil($user_t/60)," mins\n";

# PRINTING  @copies 
# @copies : fam1 rel fam2 Strand fam1_st fam1_end fam2_st fam2_end
if($copy_file_flag){
	print OUTFILECOPIES "#fam1\trel\tfam2\tstrand\tfam1-st\tf1-end\tf2-st\tf2-end\n";
	foreach $i (@pos_copies){
		print OUTFILECOPIES "$i->[0]\t$i->[1]\t$i->[2]\t$i->[3]\t$i->[4]\t$i->[5]\t$i->[6]\t$i->[7]\n";
	}
	
	foreach $i (@comp_copies){
		print OUTFILECOPIES "$i->[0]\t$i->[1]\t$i->[2]\t$i->[3]\t$i->[4]\t$i->[5]\t$i->[6]\t$i->[7]\n";
	}
	
	#calculating time taken
	($user_t,$system_t,$cuser_t,$csystem_t) = times;
	print OUTFILECOPIES "\n\# Runtime details after printing copy info: \n";
	print OUTFILECOPIES "\# System time for process: ",ceil($system_t/60)," mins\n";
	print OUTFILECOPIES "\# User time for process: ",ceil($user_t/60)," mins\n";

	
	print STDERR "\n\# Runtime details after printing copy info: \n";
	print STDERR "\# System time for process: ",ceil($system_t/60)," mins\n";
	print STDERR "\# User time for process: ",ceil($user_t/60)," mins\n";
}

close (OUTFILEDATA);
if($copy_file_flag){
	close (OUTFILECOPIES);
}

exit;

# </BODY>