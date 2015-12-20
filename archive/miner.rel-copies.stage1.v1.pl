#!/usr/bin/perl -w
# MGEL
# Surya Saha 07/25/07
# reading the copy list from out2F_itemsets for all relationships with confidence 
# > threshold and getting the extended subseq (flanking 20kb, 50kb, 100kb) from the chromosome seq
# Writing out the copy sequences for each such relationship in separate files


# NOTE: ADD IN CODE FROM DEMO.MA-CLUSTALW.PL. IT HAS THE THRESHOLD FUNCTIONALITY WORKED IN.. IDIOT!!

# merged_f_itemsets
#fam1, Category, fam2, Occ, Occ conf, Avg Rand Coin, Avg Rand Coin conf, Diff, Diff conf, Avg Dist, Std Dev, fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen
# R=139	U	R=347	23	1	0	0	23	1	2633	3833.08	R=139	62	1492	R=347	23	68
# R=347	U	R=248	23	1	0	0	23	1	86	2.44	R=347	23	68	R=248	34	1514
# R=441	U	R=347	23	1	0	0	23	1	858	2798.19	R=441	34	698	R=347	23	68
# R=1026	U	R=918	21	1	0	0	21	1	-25	2.64	R=1026	24	68	R=918	21	206

# copies
# #fam1	rel	fam2	Strand	fam1-st	f1-end	f2-st	f2-end
# R=291	d2	R=39	B	8936	9225	9987	10255
# R=291	d3	R=112	+	8936	9225	10808	10865
# R=291	d3	R=112	C	8936	9225	10836	10892
# R=291	d3	R=401	+	8936	9225	10880	11115
# R=291	d3	R=688	C	8936	9225	10966	11115
# R=291	d3	R=688	+	8936	9225	10966	11115


# v1: 

# Stage 1
# v1: reduced nof cmdline parameters
#   : adapted to disregard strand info for relationships
# v2: added in option not to print extended subseqs
# 
use strict;
use warnings;
use POSIX;

my $ver=2.0;

unless (@ARGV == 5){
	print "USAGE: $0 < merged .f_itemsets.tab file> < .copies.tab file> <chr seq file> <conf limit> <-extend/core(need flanking regions??)\n";
	exit;
}

# get the complement
sub comp{
	my $DNA;
	$DNA=$_[0];
	$DNA=~ s/\s*//g;# clean it
	$DNA=~ tr/ACGTacgt/TGCAtgca/;
	return $DNA;
}

my ($ifname,$rec,@temp,@mf_itemsets_table,@copies_table, $ctr,$i,$j,$k,
$chr_seq,$tot_recs,$tot_copies,$conf,$fam1,$fam2,$rel,$flag,$extend,
$count,@copies,$user_t,$system_t,$cuser_t,$csystem_t);

$conf=$ARGV[3];
chomp $conf;

$extend=$ARGV[4];
chomp $extend;
if(!($extend eq "-extend" || $extend eq "-core")){ print STDERR "flag can be only -extend or -core\nExiting..\n"; exit;}
# convert to 0 or 1
if($extend eq "-extend"){ $extend=1;}
elsif ($extend eq "-core"){ $extend=0;}

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
		push @mf_itemsets_table,[split(' ',$rec)];
		$ctr++;
	}	
}
# record tot recs
$tot_recs = $ctr;
close (INFILE);
print STDERR "Version: $ver\n";
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
	
	push @copies_table,[split(' ',$rec)];
	$ctr++;
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
foreach $i (@mf_itemsets_table){
	@copies=();
	$i->[0]=~ s/\s//g; $i->[1]=~ s/\s//g; $i->[2]=~ s/\s//g; $i->[7]=~ s/\s//g;

	$fam1=$i->[0]; $rel=$i->[1]; $fam2=$i->[2]; $count=$i->[7];
	
	# pull out copy subseq
	# pull out extended copy subseq(with 20kb, 50kb, 100kb flanking seqs)
	$ctr=0;
	undef @temp;
	foreach $j (@copies_table){
		# clean up input
		$j->[0]=~ s/\s//g; $j->[1]=~ s/\s//g; $j->[2]=~ s/\s//g; $j->[3]=~ s/\s//g;

		if($j->[0] eq $fam1 && $j->[1] eq $rel && $j->[2] eq $fam2){
			$temp[0]=$j->[4]; $temp[1]=$j->[5]; 
			$temp[2]=$j->[6]; $temp[3]=$j->[7];
			# sort the list of coordinates
			@temp= sort {$a <=> $b} @temp;
			# get the begg and the end to
			# get the whole chunk out
			$copies[$ctr][0]=$temp[0];
			$copies[$ctr][1]=$temp[3];
			$copies[$ctr++][2]=$j->[3];# record strand for subseq
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
	if ($extend){
		unless(open(OUTFILECOPIES5K,">".$fam1."_".$rel."_".$fam2.".5k.fas")){print "not able to open  ".$fam1."_".$rel."_".$fam2.".5k.fas\n\n";exit;}
		unless(open(OUTFILECOPIES1K,">".$fam1."_".$rel."_".$fam2.".10k.fas")){print "not able to open  ".$fam1."_".$rel."_".$fam2.".10k.fas\n\n";exit;}
# 		unless(open(OUTFILECOPIES20K,">".$fam1."_".$rel."_".$fam2.".20k.fas")){print "not able to open  ".$fam1."_".$rel."_".$fam2.".20k.fas\n\n";exit;}
# 		unless(open(OUTFILECOPIES50K,">".$fam1."_".$rel."_".$fam2.".50k.fas")){print "not able to open  ".$fam1."_".$rel."_".$fam2.".50k.fas\n\n";exit;}
# 		unless(open(OUTFILECOPIES100K,">".$fam1."_".$rel."_".$fam2.".100k.fas")){print "not able to open  ".$fam1."_".$rel."_".$fam2.".100k.fas\n\n";exit;}
	}
	
	foreach $k (@copies){# reusing $k
		# with 0 flanking
		print OUTFILECOPIES ">Subseq: ",$ctr," Start: ",$k->[0]," End: ",$k->[1]," Len: ",($k->[1]-$k->[0]),"\n";
		if($k->[2] eq '+'){
			print OUTFILECOPIES substr($chr_seq,($k->[0]-1),($k->[1] - $k->[0])),"\n";
		}
		elsif ($k->[2] eq 'C'){
			print OUTFILECOPIES comp(substr($chr_seq,($k->[0]-1),($k->[1] - $k->[0]))),"\n";
		}
		else{ # error catcher
			print STDERR "ERROR!! \n";
		}
		
		if ($extend){
			# with 5k flanking
			print OUTFILECOPIES5K ">Subseq: ",$ctr," Start: ",$k->[0]-5000," End: ",$k->[1]+5000," Len: ",($k->[1]-$k->[0] + 10000),"\n";
			if($k->[2] eq '+'){
				print OUTFILECOPIES5K substr($chr_seq,(($k->[0]-1)-5000),(($k->[1] - $k->[0])+5000)),"\n";
			}
			elsif ($k->[2] eq 'C'){
				print OUTFILECOPIES5K comp(substr($chr_seq,(($k->[0]-1)-5000),(($k->[1] - $k->[0])+5000))),"\n";
			}
			else{ # error catcher
				print STDERR "ERROR!! \n";
			}
			
			# with 1k flanking
			print OUTFILECOPIES1K ">Subseq: ",$ctr," Start: ",$k->[0]-1000," End: ",$k->[1]+1000," Len: ",($k->[1]-$k->[0] + 20000),"\n";
			if($k->[2] eq '+'){
				print OUTFILECOPIES1K substr($chr_seq,(($k->[0]-1)-1000),(($k->[1] - $k->[0])+1000)),"\n";
			}
			elsif ($k->[2] eq 'C'){
				print OUTFILECOPIES1K comp(substr($chr_seq,(($k->[0]-1)-1000),(($k->[1] - $k->[0])+1000))),"\n";
			}
			else{ # error catcher
				print STDERR "ERROR!! \n";
			}
			
# 			# with 20k flanking
# 			print OUTFILECOPIES20K ">Subseq: ",$ctr," Start: ",$k->[0]-20000," End: ",$k->[1]+20000," Len: ",($k->[1]-$k->[0] + 40000),"\n";
# 			if($k->[2] eq '+'){
# 				print OUTFILECOPIES20K substr($chr_seq,(($k->[0]-1)-20000),(($k->[1] - $k->[0])+20000)),"\n";
# 			}
# 			elsif ($k->[2] eq 'C'){
# 				print OUTFILECOPIES20K comp(substr($chr_seq,(($k->[0]-1)-20000),(($k->[1] - $k->[0])+20000))),"\n";
# 			}
# 			else{ # error catcher
# 				print STDERR "ERROR!! \n";
# 			}
# 			
# 			# with 50k flanking
# 			print OUTFILECOPIES50K ">Subseq: ",$ctr," Start: ",$k->[0]-50000," End: ",$k->[1]+50000," Len: ",($k->[1]-$k->[0] + 100000),"\n";
# 			if($k->[2] eq '+'){
# 				print OUTFILECOPIES50K substr($chr_seq,(($k->[0]-1)-50000),(($k->[1] - $k->[0])+50000)),"\n";
# 			}
# 			elsif ($k->[2] eq 'C'){
# 				print OUTFILECOPIES50K comp(substr($chr_seq,(($k->[0]-1)-50000),(($k->[1] - $k->[0])+50000))),"\n";
# 			}
# 			else{ # error catcher
# 				print STDERR "ERROR!! \n";
# 			}
# 			
# 			# with 100k flanking
# 			print OUTFILECOPIES100K ">Subseq: ",$ctr," Start: ",$k->[0]-100000," End: ",$k->[1]+100000," Len: ",($k->[1]-$k->[0] + 200000),"\n";
# 			if($k->[2] eq '+'){
# 				print OUTFILECOPIES100K substr($chr_seq,(($k->[0]-1)-100000),(($k->[1] - $k->[0])+100000)),"\n";
# 			}
# 			elsif ($k->[2] eq 'C'){
# 				print OUTFILECOPIES100K comp(substr($chr_seq,(($k->[0]-1)-100000),(($k->[1] - $k->[0])+100000))),"\n";
# 			}
# 			else{ # error catcher
# 				print STDERR "ERROR!! \n";
# 			}
		}
		$ctr++;
	}
	undef @copies;
	close (OUTFILECOPIES);
	if ($extend){
		close (OUTFILECOPIES5K);
		close (OUTFILECOPIES1K);
# 		close (OUTFILECOPIES20K);
# 		close (OUTFILECOPIES50K);
# 		close (OUTFILECOPIES100K);
	}
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
