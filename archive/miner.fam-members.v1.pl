#!/usr/bin/perl -w
# MGEL
# Surya Saha 09/28/07
# reading cmd line input fasta file and .out file which is sorted on the start position
# writing out a mfasta file for each family

#   1935 10.6  0.0  2.8 chr12|12012               8936      9225 27487989 C  R=291           Unknown             (0)     283       2
#   1009 22.1  2.2  4.2 chr12|12012               9987     10255 27486959 +  R=39            Unknown               1     264     (0)
#    226 29.3  0.0  0.0 chr12|12012              10808     10865 27486349 +  R=112           Unknown               1      58     (0)
#   1558 11.8  1.2  0.8 chr12|12012              11272     10850 27485686 C  R=39            Unknown             (6)     258       1
#   1558 11.8  1.2  0.8 chr12|12012              11272     10840 27485686 C  R=39            Unknown             (6)     258       1

# v1: 
# IS THERE A PROBLEM WITH THE AVG LEN CALCULATION??

use strict;
use warnings;
use POSIX;
use Graph;

unless (@ARGV == 2){
	print "USAGE: $0 <input .out file> <input .fas file>\n";
	exit;
}

my ($ifname,$rec,$chr_seq,@temp,@fams,@fam_info,@copy_info,%temphash,$ctr,$i,$j,
$tot_recs,$user_t,$system_t,$cuser_t,$csystem_t);

$ifname=$ARGV[0];
chomp $ifname;
unless(open(INFILE_OUT,$ifname)){print "not able to open ".$ifname."\n\n";exit;}
$ifname=$ARGV[1];
chomp $ifname;
unless(open(INFILE_FAS,$ifname)){print "not able to open ".$ifname."\n\n";exit;}

# get the complement
sub comp{
	my $DNA;
	$DNA=$_[0];
	$DNA=~ s/\s*//g;# clean it
	$DNA=~ tr/ACGTacgt/TGCAtgca/;
	return $DNA;
}

#to get the index position of a family in the @fam_info array
#it might be faster to just get info from the @counts array 
#once we have the index pos
#params: $fam
sub get_index{
	my ($fam,$ctr);
	$fam=$_[0];
	$fam=~ s/\s*//g;
	$ctr=0;
	foreach (@fam_info){
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

# slurping in the whole freq itemsets file
$ctr=0;
while($rec=<INFILE_OUT>){
	if($rec =~ /^#/){next;}
	if(length ($rec) < 10){next;}#for avoiding last line
	
	push @copy_info,[split(' ',$rec)];
	$fams[$ctr]=$copy_info[$ctr][9];# get the family
	$ctr++;
}
close (INFILE_OUT);
# record tot recs
$tot_recs = $ctr;
$chr_seq='';
# get the chr sequence
while($rec=<INFILE_FAS>){
	if($rec =~ /^>/){next;}
	if(length ($rec) < 10){next;}#for avoiding last line
	$rec=~ s/\s*//g;
	$chr_seq=$chr_seq.$rec;
}
close (INFILE_FAS);

# sort and get the uniq list of the family
@temp=sort @fams;
@fams=@temp;
%temphash = map { $_, 1 } @fams;
@fams = keys %temphash;
#sorting agn
@temp=sort @fams;
@fams=@temp;

# set the copy count, avg copy length of each family to 0
# @fa_info = fam name, copy count, avg copy length 
for($ctr=0;$ctr<@fams;$ctr++){
	$fam_info[$ctr][0]=$fams[$ctr];
	$fam_info[$ctr][1]=0;
	$fam_info[$ctr][2]=0;
}

# 0..13
#   1935 10.6  0.0  2.8 chr12|12012               8936      9225 27487989 C  R=291           Unknown             (0)     283       2
#   1009 22.1  2.2  4.2 chr12|12012               9987     10255 27486959 +  R=39            Unknown               1     264     (0)
#    226 29.3  0.0  0.0 chr12|12012              10808     10865 27486349 +  R=112           Unknown               1      58     (0)
#    241 26.3  0.0  0.0 chr12|12012              10836     10892 27486322 C  R=112           Unknown             (0)      58       2
#    239 29.8  4.2  8.8 chr12|12012              10880     11115 27486099 +  R=401           Unknown               1     226     (0)
#   1108  9.5  0.0  1.4 chr12|12012              10966     11115 27486099 C  R=688           Unknown             (0)     148       1
#   1558 11.8  1.2  0.8 chr12|12012              11272     11528 27485686 C  R=39            Unknown             (6)     258       1
#   1688  8.6  0.0  0.8 chr12|12012              13726     13983 27483231 C  R=443           Unknown            (31)     256       1
#    944 17.1  0.0  0.0 chr12|12012              13873     14059 27483155 C  R=499           Unknown             (1)     187       1

# now sort it by fam 
@temp = sort {($a->[9] cmp $b->[9])} @copy_info;
@copy_info=@temp;

# write out the family copy sequences
$ctr=0;
foreach(@copy_info){
	if($ctr==0) {
		$i=$copy_info[$ctr][9];# init the family name
		$j=1;# init the copy num
		
		unless(open(OUTFILE,">>$i.fas")){print "not able to open $i.fas\n\n";exit;}
	}
	
	if($copy_info[$ctr][9] eq $i){ # if current family
		print OUTFILE ">Copy:",$j++," sub%:",$copy_info[$ctr][1]," del%:",$copy_info[$ctr][2]," ins%:",$copy_info[$ctr][3]," len:",($copy_info[$ctr][6]-$copy_info[$ctr][5])," Strand:",$copy_info[$ctr][8]," Start:",$copy_info[$ctr][5]," End:",$copy_info[$ctr][6],"\n";
		# print the sequence
		if($copy_info[$ctr][8] eq '+'){
			print OUTFILE substr($chr_seq,$copy_info[$ctr][5]-1,($copy_info[$ctr][6] - $copy_info[$ctr][5])),"\n";
		}
		elsif($copy_info[$ctr][8] eq 'C'){
			# add seq info
			print OUTFILE  comp(substr($chr_seq,$copy_info[$ctr][5]-1,($copy_info[$ctr][6] - $copy_info[$ctr][5]))),"\n";
		}
		
		# increment count and length
		$fam_info[get_index($copy_info[$ctr][9])][1]++;
		$fam_info[get_index($copy_info[$ctr][9])][2]=$fam_info[get_index($copy_info[$ctr][9])][2] + ($copy_info[$ctr][6]-$copy_info[$ctr][5]);

	}
	elsif($copy_info[$ctr][9] ne $i){ # if new family
		print OUTFILE "\n";
		close OUTFILE;
		$i=$copy_info[$ctr][9];# init the family name
		$j=1;# init the copy num
		
		unless(open(OUTFILE,">>$i.fas")){print "not able to open $i.fas\n\n";exit;}

		print OUTFILE ">Copy:",$j++," sub%:",$copy_info[$ctr][1]," del%:",$copy_info[$ctr][2]," ins%:",$copy_info[$ctr][3]," len:",($copy_info[$ctr][6]-$copy_info[$ctr][5])," Strand:",$copy_info[$ctr][8]," Start:",$copy_info[$ctr][5]," End:",$copy_info[$ctr][6],"\n";
		# print the sequence
		if($copy_info[$ctr][8] eq '+'){
			print OUTFILE substr($chr_seq,$copy_info[$ctr][5]-1,($copy_info[$ctr][6] - $copy_info[$ctr][5])),"\n";
		}
		elsif($copy_info[$ctr][8] eq 'C'){
			# add seq info
			print OUTFILE  comp(substr($chr_seq,$copy_info[$ctr][5]-1,($copy_info[$ctr][6] - $copy_info[$ctr][5]))),"\n";
		}
		
		# increment count and length
		$fam_info[get_index($copy_info[$ctr][9])][1]++;
		$fam_info[get_index($copy_info[$ctr][9])][2]= $fam_info[get_index($copy_info[$ctr][9])][2]+ ($copy_info[$ctr][6]-$copy_info[$ctr][5]);
	}
	else{ # error catcher
		print STDERR "ERROR!! \n";
	}
	$ctr++;
}

close OUTFILE;

unless(open(OUTFILE,">$ifname.fam.info")){print "not able to open $ifname.fam.info\n\n";exit;}

# set the avg length and write it out
# @fa_info = fam name, copy count, avg copy length 
for($ctr=0;$ctr<@fams;$ctr++){
	$fam_info[$ctr][2]=int($fam_info[$ctr][2]/$fam_info[$ctr][1]);
}

# now sort it by count and avg length
@temp = sort {($b->[1] <=> $a->[1])} @fam_info;
@fam_info=@temp;

# write out the family info
print OUTFILE "Version : 1.0\n\n";
print OUTFILE "RepeatScout Version : 1.0.1\n";
print OUTFILE "Total families : ",scalar @fams,"\n";
print OUTFILE "Total copies : ",$tot_recs,"\n\n";
print OUTFILE "#Fam\tCount\tAvglen\n";
for($ctr=0;$ctr<@fams;$ctr++){
	print OUTFILE "$fam_info[$ctr][0]\t$fam_info[$ctr][1]\t$fam_info[$ctr][2]\n";
}
close OUTFILE;


# calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print "\# Runtime details : \n";
print "\# System time for process: $system_t\n";
# print OUTFILEDATA "\# System time for children: $csystem_t\n";
print "\# User time for process: $user_t\n";
# print OUTFILEDATA "\# User time for children: $cuser_t\n";


exit;
