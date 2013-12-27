#!/usr/bin/perl -w
# MGEL
# Surya Saha 07/25/07
# reading the copy list from out2F_itemsets for all relationships with confidence 
# > threshold and getting the extended subseq (flanking 20kb, 50kb, 100kb) from the chromosome seq
# Writing out the copy sequences for each such relationship in separate files

# f_itemsets
# #fam1, Category, fam2, Occ, Occ conf, Avg Rand Coin, Avg Rand Coin conf, Diff, Diff conf, Strand, fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen
# R=501	u2	R=468	8	0.8	0	0	8	0.8	+	R=501	10	166	R=468	262	149
# R=810	u3	R=394	8	0.72	0	0	8	0.72	+	R=810	11	224	R=394	15	289
# R=748	u3	R=225	7	0.7	0	0	7	0.7	+	R=748	10	241	R=225	22	457

# copies
# #fam1	rel	fam2	Strand	fam1-st	f1-end	f2-st	f2-end
# R=291	d2	R=39	B	8936	9225	9987	10255
# R=291	d3	R=112	B	8936	9225	10808	10865
# R=291	d3	R=112	C	8936	9225	10836	10892
# R=291	d3	R=401	B	8936	9225	10880	11115
# R=291	d3	R=688	C	8936	9225	10966	11115
# R=291	d3	R=688	B	8936	9225	10966	11115


# v1: 

use strict;
use warnings;
use POSIX;

unless (@ARGV == 5){
	print "USAGE: $0 < merged .f_itemsets.tab file> < .copies.tab file> <chr seq file> <conf limit> <strand>\n";
	exit;
}

print STDERR "Warning: Make sure merged f_itemsets file is for given strand!!\n";
print STDERR "To Add: Print out annotations for families, if any\n\n";
# get the complement
sub comp{
	my $DNA;
	$DNA=$_[0];
	$DNA=~ s/\s*//g;# clean it
	$DNA=~ tr/ACGTacgt/TGCAtgca/;
	return $DNA;
}

my ($ifname,$rec,@temp,@f_itemsets_table,@copies_table, $ctr,$i,$j,$k,
$chr_seq,$tot_recs,$tot_copies,$conf,$strand,$fam1,$fam2,$rel,$flag,
$count,$user_t,$system_t,$cuser_t,$csystem_t);

$conf=$ARGV[3];
chomp $conf;
$strand=$ARGV[4];
chomp $strand;
unless(($strand eq '+') || ($strand eq 'C') || ($strand eq 'B') ){print "Strand not valid!!\n\n";exit;}

$ifname=$ARGV[0];
chomp $ifname;
unless(open(INFILE,$ifname)){print "not able to open ".$ifname."\n\n";exit;}

# slurping in the whole merged freq itemsets file for selected strand
$ctr=0;
while($rec=<INFILE>){
	if($rec =~ /^#/){next;}
	if(length ($rec) < 10){next;}#for avoiding last line
	@temp=split(' ',$rec);
	
	if($temp[8] >= $conf){# get the ones >= conf
		push @f_itemsets_table,[split(' ',$rec)];
		$ctr++;
	}	
}
# record tot recs
$tot_recs = $ctr;
close (INFILE);
print STDERR "Read in merged freq itemsets file, nof records: $ctr\n";

$ifname=$ARGV[1];
chomp $ifname;
unless(open(INFILE,$ifname)){print "not able to open ".$ifname."\n\n";exit;}

# slurping in the records from copies file for given strand
$ctr=0;
while($rec=<INFILE>){
	if($rec =~ /^#/){next;}
	if(length ($rec) < 10){next;}#for avoiding last line
	@temp=split(' ',$rec);
	
	# this may not be needed 
	# we have separate files for each strand
	if($temp[3] eq $strand){# get the ones for given strand only
		push @copies_table,[split(' ',$rec)];
		$ctr++;
	}
}
# record tot recs
$tot_copies = $ctr;
close (INFILE);
print STDERR "Read in copies file, nof copies : $ctr\n";

$ifname=$ARGV[2];
chomp $ifname;
unless(open(INFILE,$ifname)){print "not able to open ".$ifname."\n\n";exit;}

# slurping in the whole chromosome file
while($rec=<INFILE>){
	if($rec =~ /^>/){next;}
	if(length ($rec) < 10){next;}#for avoiding last line
	$rec=~ s/\s*//g;
	$chr_seq=$chr_seq.$rec;
}
close (INFILE);
print STDERR "Read in chromosome file, length : ",length($chr_seq),"\n";

# unless(open(OUTFILECONSEQ,">$ifname.conseq.fas")){print "not able to open ".$ifname."conseq.fas\n\n";exit;}

# PROCESS RELATIONSHIPS
foreach $i (@f_itemsets_table){
	my (@copies);
	$i->[0]=~ s/\s//g; $i->[1]=~ s/\s//g; $i->[2]=~ s/\s//g;
  $i->[7]=~ s/\s//g;

	$fam1=$i->[0]; $rel=$i->[1]; 
	$fam2=$i->[2]; $count=$i->[7];
	
	# Note: this may not apply as the reciprocal relationship might have been screened out
	# as having a conf value < threshold when the merged f_itemset was being read in
	# so the copies for the reci. rel. may not have been read into @f_itemsets_table
# 	if ($rel ne "SKIP"){# check if this is an marked record
		#look for reciprocals for U & D relations, only count will be same
		
# 		$flag=0;
# 		if($rel=~ /u/ || $rel=~ /d/){
# 			$k=$rel;
# 			$k=~ tr/ud/du/;# changing to reciprocal relation
# 			$ctr=0;
# 			foreach $j (@f_itemsets_table){
# 				if($j->[0] eq $fam2 && $j->[1] eq $k && 
# 				$j->[2] eq $fam1 && $j->[7] == $count){
# 					$flag=1;
# 					$f_itemsets_table[$ctr][1]="SKIP"; #mark record for skipping later
# 					last;
# 				}
# 				$ctr++;
# 			}
# 		}
# 		
# 		# if no reciprocal
# 		if($flag==0){
# 			print STDERR "No reciprocal found for $fam1 $rel $fam2 with count $count\n";
# 			next; # DO WE CONTINUE??
# 		}

		# look for annotation TODO!!
		
	# pull out copy clusters
	# pull out extended copy clusters (with 20kb, 50kb, 100kb flanking seqs)
	$ctr=0;
	undef @temp;
	foreach $j (@copies_table){
		# clean up input
		$j->[0]=~ s/\s//g; $j->[1]=~ s/\s//g; $j->[2]=~ s/\s//g;
		$j->[0]=~ s/\s//g; $j->[1]=~ s/\s//g; $j->[2]=~ s/\s//g;

		if($j->[0] eq $fam1 && $j->[1] eq $rel && $j->[2] eq $fam2){
			$temp[0]=$j->[4]; $temp[1]=$j->[5]; 
			$temp[2]=$j->[6]; $temp[3]=$j->[7];
			# sort the list of coordinates
			@temp= sort {$a <=> $b} @temp;
			# get the begg and the end to
			# get the whole chunk out
			$copies[$ctr][0]=$temp[0];
			$copies[$ctr++][1]=$temp[3];
		}
	}
	if ($ctr==0){ 
		print STDERR "No copies/clusters found for $fam1 $rel $fam2\n";
		next;
	}
	
	$tot_copies=$ctr;
	
	# write out the copies
	$ctr=0;
	unless(open(OUTFILECOPIES,">".$fam1."_".$rel."_".$fam2.".fas")){print "not able to open  ".$fam1."_".$rel."_".$fam2.".fas\n\n";exit;}
	unless(open(OUTFILECOPIES5K,">".$fam1."_".$rel."_".$fam2.".5k.fas")){print "not able to open  ".$fam1."_".$rel."_".$fam2.".5k.fas\n\n";exit;}
	unless(open(OUTFILECOPIES10K,">".$fam1."_".$rel."_".$fam2.".10k.fas")){print "not able to open  ".$fam1."_".$rel."_".$fam2.".10k.fas\n\n";exit;}
	unless(open(OUTFILECOPIES20K,">".$fam1."_".$rel."_".$fam2.".20k.fas")){print "not able to open  ".$fam1."_".$rel."_".$fam2.".20k.fas\n\n";exit;}
	unless(open(OUTFILECOPIES50K,">".$fam1."_".$rel."_".$fam2.".50k.fas")){print "not able to open  ".$fam1."_".$rel."_".$fam2.".50k.fas\n\n";exit;}
	unless(open(OUTFILECOPIES100K,">".$fam1."_".$rel."_".$fam2.".100k.fas")){print "not able to open  ".$fam1."_".$rel."_".$fam2.".100k.fas\n\n";exit;}

	# TODO : ADD THE LENGTH OF EACH CLUSTER TO THE FASTA HEADER
	foreach $k (@copies){# reusing $k
		# with 0 flanking
		print OUTFILECOPIES ">Cluster: ",$ctr," Start: ",$k->[0]," End: ",$k->[1]," Len: ",($k->[1]-$k->[0]),"\n";
		if($strand eq '+'){
			print OUTFILECOPIES substr($chr_seq,($k->[0]-1),($k->[1] - $k->[0])),"\n";
		}
		elsif ($strand eq 'C'){
			print OUTFILECOPIES comp(substr($chr_seq,($k->[0]-1),($k->[1] - $k->[0]))),"\n";
		}
		else{ # error catcher
			print STDERR "ERROR!! \n";
		}
		
		# with 5k flanking
		print OUTFILECOPIES5K ">Cluster: ",$ctr," Start: ",$k->[0]-5000," End: ",$k->[1]+5000," Len: ",($k->[1]-$k->[0] + 10000),"\n";
		if($strand eq '+'){
			print OUTFILECOPIES5K substr($chr_seq,(($k->[0]-1)-5000),(($k->[1] - $k->[0])+10000)),"\n";
		}
		elsif ($strand eq 'C'){
			print OUTFILECOPIES5K comp(substr($chr_seq,(($k->[0]-1)-5000),(($k->[1] - $k->[0])+10000))),"\n";
		}
		else{ # error catcher
			print STDERR "ERROR!! \n";
		}
		
		# with 10k flanking
		print OUTFILECOPIES10K ">Cluster: ",$ctr," Start: ",$k->[0]-10000," End: ",$k->[1]+10000," Len: ",($k->[1]-$k->[0] + 20000),"\n";
		if($strand eq '+'){
			print OUTFILECOPIES10K substr($chr_seq,(($k->[0]-1)-10000),(($k->[1] - $k->[0])+20000)),"\n";
		}
		elsif ($strand eq 'C'){
			print OUTFILECOPIES10K comp(substr($chr_seq,(($k->[0]-1)-10000),(($k->[1] - $k->[0])+20000))),"\n";
		}
		else{ # error catcher
			print STDERR "ERROR!! \n";
		}
		
		# with 20k flanking
		print OUTFILECOPIES20K ">Cluster: ",$ctr," Start: ",$k->[0]-20000," End: ",$k->[1]+20000," Len: ",($k->[1]-$k->[0] + 40000),"\n";
		if($strand eq '+'){
			print OUTFILECOPIES20K substr($chr_seq,(($k->[0]-1)-20000),(($k->[1] - $k->[0])+40000)),"\n";
		}
		elsif ($strand eq 'C'){
			print OUTFILECOPIES20K comp(substr($chr_seq,(($k->[0]-1)-20000),(($k->[1] - $k->[0])+40000))),"\n";
		}
		else{ # error catcher
			print STDERR "ERROR!! \n";
		}
		
		# with 50k flanking
		print OUTFILECOPIES50K ">Cluster: ",$ctr," Start: ",$k->[0]-50000," End: ",$k->[1]+50000," Len: ",($k->[1]-$k->[0] + 100000),"\n";
		if($strand eq '+'){
			print OUTFILECOPIES50K substr($chr_seq,(($k->[0]-1)-50000),(($k->[1] - $k->[0])+100000)),"\n";
		}
		elsif ($strand eq 'C'){
			print OUTFILECOPIES50K comp(substr($chr_seq,(($k->[0]-1)-50000),(($k->[1] - $k->[0])+100000))),"\n";
		}
		else{ # error catcher
			print STDERR "ERROR!! \n";
		}
		
		# with 100k flanking
		print OUTFILECOPIES100K ">Cluster: ",$ctr," Start: ",$k->[0]-100000," End: ",$k->[1]+100000," Len: ",($k->[1]-$k->[0] + 200000),"\n";
		if($strand eq '+'){
			print OUTFILECOPIES100K substr($chr_seq,(($k->[0]-1)-100000),(($k->[1] - $k->[0])+200000)),"\n";
		}
		elsif ($strand eq 'C'){
			print OUTFILECOPIES100K comp(substr($chr_seq,(($k->[0]-1)-100000),(($k->[1] - $k->[0])+200000))),"\n";
		}
		else{ # error catcher
			print STDERR "ERROR!! \n";
		}
		$ctr++;
	}
	undef @copies;
	close (OUTFILECOPIES);
	close (OUTFILECOPIES5K);
	close (OUTFILECOPIES10K);
	close (OUTFILECOPIES20K);
	close (OUTFILECOPIES50K);
	close (OUTFILECOPIES100K);
	print STDERR "Wrote file for $fam1 $rel $fam2 relationship...\n";
# 	}# if end
}# foreach end

# close (OUTFILECONSEQ);

# calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print "\# Runtime details : \n";
print "\# System time for process: $system_t\n";
# print OUTFILEDATA "\# System time for children: $csystem_t\n";
print "\# User time for process: $user_t\n";
# print OUTFILEDATA "\# User time for children: $cuser_t\n";



exit;
