#!/usr/bin/perl -w
# MGEL
# Surya Saha 04/06/07
# reading cmd line input RM .cat file for the novel repeat library
# and find the repeat families that are annotated for > 40% of their total length
# and print out their annotation

# Version 2.0 04/16/07
# Counting the exact number of bases annotated in a family

# Version 3.0 03/06/08
# Changing the format of output to make it readable by the graph mining script

# Version 3.5  06/10/08
# Problem with R=2 fixed.  Incorrect count reports for single-element array?
# Multiple cat files can be used.

# Version 4.5 06/16/08 Jonathan
# Added attribute listing source cat file for each annotation.  Bugfixed & removal of unnecessary temp. variables.

# Version 4.6 06/24/08 Surya
# Reporting rounded values for annotated ratio

# Version 4.7 03/4/09 Surya
# Adding in stats for all classes of RepBase repeats
# Only LOGIC PROBLEM is that if the CRF has multiple annotations of which I am only counting 1, e.g. 1 En-Spm and many DNA 
# annots, I will presume that the CRF is a En-Spm when it could be something else


use strict;
use warnings;
use POSIX; #for floor() and ciel()
# use List::MoreUtils qw(uniq); # not installed

sub round2{
	my ($num);
	$num=$_[0];
	$num=$num*100;
	$num=int($num);
	$num=$num/100;
	return $num;
}

my ($ver);
$ver="4.7";

unless (@ARGV >= 3){
	print "USAGE: $0 <annotation threshold%> <rep lib fasta file> <.cat files from RM>\n";
	exit;
}

my ($ifname1,$ifname2,$rec,@cat_table,@faminfo,$DNA,$thold,@temp,@rep_types,
$ctr,$i,$j,$k,$tot_recs,$fams_annotated,$tot_seqs,$user_t,$system_t,%temp_hash,
$cuser_t,$csystem_t,$tmpdata,@fam_annot);

$ifname1=$ARGV[1];
chomp $ifname1;
unless(open(INFILEFASTA,$ifname1)){print "not able to open ".$ifname1."\n\n";exit;}

$thold=$ARGV[0];
chomp $thold;
unless(open(OUTFILEDATA,">$ifname1.$thold.annot")){print "not able to open ".$ifname1.".tab \n\n";exit;}

$thold = int ($thold);
$thold = $thold * 0.01;


# slurping in the whole report file
$ctr=0;

for (2..$#ARGV) {
	my($skipflag);
	$skipflag=0;
	$ifname2=$ARGV[$_];
	chomp $ifname2;
	unless(open(INFILECAT,$ifname2)){print "not able to open ".$ifname2."\n\n";exit;}
	while($rec=<INFILECAT>){
		$ifname2=$ARGV[$_];
		if($skipflag != 0 && $skipflag) {--$skipflag;next;}
		if($rec =~ /^#/){next;}
		# RepeatMasker version open-3.1.6 , sensitive mode
		# run with blastp version 2.0MP-WashU
		# RepBase Update 20061006, RM database version 20061006
		if ($rec=~ /RepeatMasker/){$skipflag=2; next;}
	
		if(length ($rec) < 10){next;}#for avoiding last line
		$ifname2 =~ s/\s//g;
		if($ifname2 =~ m/\/([^\W]*)$/g) {$rec = $rec." $1";}
		else{$rec="$rec $ifname2";}
		push(@cat_table,[split(' ',$rec)]);
		# get the class of the annotations
		if($cat_table[$ctr][8] eq 'C'){ $i=$cat_table[$ctr][9];}
		else { $i=$cat_table[$ctr][8];}
		@temp=();
		@temp=split(/\#/,$i);
		push(@rep_types,$temp[1]); # put in LTR/Copia of COPIA1-I_OS#LTR/Copia
		$ctr++;
	}
}
# record tot recs
$tot_recs = $ctr;

#sorting @cat_table on famname
@cat_table = sort {$a->[4] cmp $b->[4]} @cat_table;

# reading in the fasta file and the lengths of each family
$ctr=0;
$DNA='';
$tmpdata='';
while($DNA=<INFILEFASTA>){
	chomp $DNA;
	
	# >R=1 (RR=2.  TRF=0.005 NSEG=0.110)
	if($DNA=~ />/ and $ctr==0) {#first name
		$DNA=~ s/\(\S*\s*\S*\s*\S*\)//g;
		#$DNA=~ s/\([\s]*//g;
		#removing leading >
		$DNA=~ s/>//g;
		# print $DNA,"\n";
		$faminfo[$ctr++][0]=$DNA;
	}	
	elsif($DNA=~ />/ and $ctr > 0){#another name
		#remove whitespace
		$tmpdata=~ s/\s//g;
		$faminfo[$ctr-1][1]= length $tmpdata;
		
		#reinit for next sequence
		$tmpdata='';
		# for next seq
		$DNA=~ s/\(\S*\s*\S*\s*\S*\)//g;
		$DNA=~ s/>//g;
		# print $DNA,"\n";
		$faminfo[$ctr++][0]=$DNA;
	}
	else{
		$tmpdata=$tmpdata.$DNA;
	}
}

#for the last sequence
#remove whitespace
$tmpdata=~ s/\s//g;
$faminfo[$ctr-1][1]= length $tmpdata;

$tot_seqs=$ctr;

# prep annotation hash
%temp_hash = map { $_, 1 } @rep_types;
@rep_types = keys %temp_hash;

my (%CRF_annot_counts);
foreach $i (@rep_types){ $CRF_annot_counts{$i}=0;}

print OUTFILEDATA "\# version: $ver\n";
print OUTFILEDATA "\# Total sequences in library: $tot_seqs\n";
print OUTFILEDATA "\# Total records from RMRB catfile: $tot_recs\n";


#    SW  perc perc perc  query     position in query     matching        repeat         
# score  div. del. ins.  sequence  begin   end (left)   repeat          class/family   
# position in  repeat
# begin   end  (left)  ID
# 0     1     2   3    4    5  6    7  8  9                 10     11   12 13
# 41000 4.27 1.54 2.10 R=19 1 5454 (0) C SPMLIKE#DNA/En-Spm (5565) 5424 1 5
# 953 13.92 0.00 0.00 R=20 1 158 (0) MERMITE18E#DNA 1369 1526 (2) 5

# iterate through $faminfo reading info for each family
# from @cat_table , print if annotated length > threshold of total length
$fams_annotated=0;
foreach $i (@faminfo){
	my ($ref_fam, $ref_len, $annot_class, $annot_ratio, $annot_len, @annot_arr, $flag);
	$ref_fam=$i->[0];
	#cleaning up
	$ref_fam=~ s/\s//g;
	$ref_len=$i->[1];
	$annot_class=$annot_len=$flag=$ctr=0;
	@fam_annot=undef;
	@annot_arr=undef;
	
	foreach $j (@cat_table){
		#cleaning up
		$j->[4]=~ s/\s//g;
		if ($ref_fam eq $j->[4]){# if rec for family
			$flag=1;
			$fam_annot[$ctr++]=$j;
			#$annot_len+=$j->[6] - $j->[5];
		}
		# if rec NOT for family
		elsif($ref_fam ne $j->[4] && $flag==1){
			last;
		}
	}
# 	17406 3.02 0.20 0.15 R=1 1 1984 (8686) C RIRE3A_I#LTR/Gypsy (2261) 1985 1 5
# 	7065 16.11 3.38 3.66 R=1 2010 3517 (7153) TRUNCATOR#LTR/Gypsy 27 1530 (1388) 5
# 	9364 10.77 0.71 1.58 R=1 3581 4985 (5685) TRUNCATOR#LTR/Gypsy 1523 2915 (3) 5
# 	48930 4.39 0.18 0.16 R=1 4984 10670 (0) C RIRE8B_I#LTR/Gypsy (0) 5981 294 5

	# now calculate the exact number of bases annotated
	# if it had annotations
	if ($fam_annot[0]){
		# initialize to 0
		for (0.. ($ref_len-1)) { $annot_arr[$_]=0};
		
		# flag for annotation and get annotation type
		my (%fam_annot_counts);
		foreach $j (@rep_types){ $fam_annot_counts{$j}=0;}
	
		foreach $j (@fam_annot){
			for (($j->[5]-1) .. ($j->[6]-1)){ $annot_arr[$_]=1;}
			
			# get type of annot
			if($j->[8] eq 'C'){
				@temp=();
				@temp=split(/\#/,$j->[9]);
				# add length to the class counter in hash
				$fam_annot_counts{$temp[1]}=$fam_annot_counts{$temp[1]}+($j->[6]-$j->[5]);
			}
			else{
				@temp=();
				@temp=split(/\#/,$j->[8]);
				# add length to the class counter in hash
				$fam_annot_counts{$temp[1]}=$fam_annot_counts{$temp[1]}+($j->[6]-$j->[5]);
			}
		}
		
		# count annotated length
		for (0 .. $#annot_arr) {
			if ($annot_arr[$_] == 1) {$annot_len++;}
		}
		$annot_ratio= &round2(($annot_len/$ref_len) * 100);
		
		# to get the max annotated length
		@temp=();
		# sorting the hash by value
		@temp = sort { $fam_annot_counts{$b} <=> $fam_annot_counts{$a} } keys %fam_annot_counts;
		$CRF_annot_counts{$temp[0]}++; 
		$annot_class=$temp[0];
		%fam_annot_counts=();
	}
	
	# if annot is > threshold then write out the faminfo and annotations
	if ($annot_len > int(ceil($thold * $ref_len))){
		# Family Len Nof annotations %annotated
		print OUTFILEDATA "$ref_fam  $ref_len  $ctr $annot_ratio% $annot_class";
		
		print OUTFILEDATA "\n";
		$fams_annotated++;# counting nof families that are annot
		#sorting on start pos of annotation
		@fam_annot = sort {$a->[5] <=> $b->[5]} @fam_annot;
		foreach (@fam_annot){
			# for records with annot from comp strand, bcos of add. col
			if ($_->[8] eq 'C'){
				print OUTFILEDATA "$ref_fam $_->[0] $_->[1] $_->[2] $_->[3] $_->[4] $_->[5] $_->[6] $_->[7] $_->[8] $_->[9] $_->[10] $_->[11] $_->[12] $_->[13] $_->[14] \n";
			}
			else{ #for the others
				print OUTFILEDATA "$ref_fam $_->[0] $_->[1] $_->[2] $_->[3] $_->[4] $_->[5] $_->[6] $_->[7] $_->[8] $_->[9] $_->[10] $_->[11] $_->[12] $_->[13] \n";
			}
		}
		print OUTFILEDATA "\n";
	}
}

print OUTFILEDATA "\# Annotation statistics : \n";
print OUTFILEDATA "\# Families annotated: $fams_annotated\n";
print OUTFILEDATA "\# Statistics for Repbase ONLY!!: \n";
print OUTFILEDATA "\# Class\tNofCRFs\n";

while ( ($i,$j) = each %CRF_annot_counts ){
		print OUTFILEDATA "\# $i\t$j\n";
}


#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\n\n\# Runtime details : \n";
print OUTFILEDATA "\# System time for process: $system_t\n";
#print OUTFILEDATA "\# System time for children: $csystem_t\n";
print OUTFILEDATA "\# User time for process: $user_t\n";
#print OUTFILEDATA "\# User time for children: $cuser_t\n";


close (INFILEFASTA);
close (INFILECAT);
close (OUTFILEDATA);



exit;
