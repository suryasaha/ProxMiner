#!/usr/bin/perl -w
# MGEL
# Surya Saha 04/06/07
# reading cmd line input RM .cat file for the novel repeat library
# and find the repeat families that are annotated for > 40% of their total length
# and print out their annotation

# Version 2.0 04/16/07
# Counting the exact number of bases annotated in a family

use strict;
use warnings;
use POSIX; #for floor() and ciel()

unless (@ARGV == 2){
	print "USAGE: $0 <rep lib fasta file> <.cat file from RM>\n";
	exit;
}


my ($ifname1,$ifname2,$rec,@temp,@cat_table,@faminfo,$DNA,
$ctr,$i,$j,$tot_recs,$tot_seqs,$user_t,$system_t,$cuser_t,$csystem_t,
$tmpdata,@fam_annot);

$ifname1=$ARGV[0];
chomp $ifname1;
unless(open(INFILEFASTA,$ifname1)){print "not able to open ".$ifname1."\n\n";exit;}
$ifname2=$ARGV[1];
chomp $ifname2;
unless(open(INFILECAT,$ifname2)){print "not able to open ".$ifname2."\n\n";exit;}

unless(open(OUTFILEDATA,">$ifname1.annot")){print "not able to open ".$ifname1.".tab \n\n";exit;}

# slurping in the whole report file
$ctr=0;
while($rec=<INFILECAT>){
	if($rec =~ /^#/){next;}
	# RepeatMasker version open-3.1.6 , sensitive mode
	# run with blastp version 2.0MP-WashU
	# RepBase Update 20061006, RM database version 20061006
	if ($rec=~ /RepeatMasker/){ last;}

	if(length ($rec) < 10){next;}#for avoiding last line
	@temp = [split(' ',$rec)];
	
	push @cat_table,@temp;
	$ctr++;
}
# record tot recs
$tot_recs = $ctr;

#sorting @cat_table on famname
@temp = sort {$a->[4] cmp $b->[4]} @cat_table;
@cat_table=@temp;

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

print OUTFILEDATA "\# Version 2.0 04/16/07\n";
print OUTFILEDATA "\# Total sequences: $tot_seqs\n";
print OUTFILEDATA "\# Total records: $tot_recs\n";


#    SW  perc perc perc  query     position in query     matching        repeat         
# score  div. del. ins.  sequence  begin   end (left)   repeat          class/family   
# position in  repeat
# begin   end  (left)  ID
# 0     1     2   3    4    5  6    7  8  9                 10     11   12 13
# 41000 4.27 1.54 2.10 R=19 1 5454 (0) C SPMLIKE#DNA/En-Spm (5565) 5424 1 5
# 953 13.92 0.00 0.00 R=20 1 158 (0) MERMITE18E#DNA 1369 1526 (2) 5

# iterate through $faminfo reading info for each family
# from @cat_table , record if annotated length > 40% of total length
foreach $i (@faminfo){
	my ($ref_fam, $ref_len, $annot_ratio, $annot_len, @annot_arr, $flag);
	$ref_fam=$i->[0];
	#cleaning up
	$ref_fam=~ s/\s//g;
	$ref_len=$i->[1];
	$annot_len=0;
	$flag=0;
	$ctr=0;
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
	if ($#fam_annot != 0){
		# initialize to 0
		for (0.. ($ref_len-1)) { $annot_arr[$_]=0};
		
		# flag for annotation
		foreach $j (@fam_annot){
			for (($j->[5]-1) .. ($j->[6]-1)){ $annot_arr[$_]=1;}
		}
		
		# count annotated length
		for (0 .. $#annot_arr) {
			if ($annot_arr[$_] == 1) {$annot_len++;}
		}
		
		$annot_ratio= ($annot_len/$ref_len) * 100;
	}
	
	# if annot is > threshold then write out the faminfo and annotations
	if ($annot_len > int(ceil(0.4 * $ref_len))){
		print OUTFILEDATA "Family:$ref_fam  Len:$ref_len  Nof annotations:$ctr \%annotated:$annot_ratio\n";
		#sorting on start pos of annotation
		@temp = sort {$a->[5] <=> $b->[5]} @fam_annot;
		@fam_annot=@temp;
		foreach (@fam_annot){
			# for records with annot from comp strand
			if ($_->[8] eq 'C'){
				print OUTFILEDATA "$_->[0] $_->[1] $_->[2] $_->[3] $_->[4] $_->[5] $_->[6] $_->[7] $_->[8] $_->[9] $_->[10] $_->[11] $_->[12] $_->[13]\n";
			}
			else{ #for the others
				print OUTFILEDATA "$_->[0] $_->[1] $_->[2] $_->[3] $_->[4] $_->[5] $_->[6] $_->[7] $_->[8] $_->[9] $_->[10] $_->[11] $_->[12]\n";
			}
		}
		print OUTFILEDATA "\n";
	}
}


#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\# Runtime details : \n";
print OUTFILEDATA "\# System time for process: $system_t\n";
#print OUTFILEDATA "\# System time for children: $csystem_t\n";
print OUTFILEDATA "\# User time for process: $user_t\n";
#print OUTFILEDATA "\# User time for children: $cuser_t\n";


close (INFILEFASTA);
close (INFILECAT);
close (OUTFILEDATA);



exit;
