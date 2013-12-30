#!/usr/bin/perl -w
#RepeatMasker .OUT Summary
#Using the .OUT file, the related FASTA file, and an ANNOT file,
#summarize .OUT file data and provide annotations involved along with it.
#Summary of total data near the EOF
#
#Jonathan Harper 6/24/2008
# Use of arrays increases the runtime


use strict;
use warnings;

#Round off to 2 decimal places
sub round2{
	my ($num);
	$num=$_[0];
	$num=$num*100;
	$num=int($num);
	$num=$num/100;
	return $num;
}

#Calculate the standard deviation of the entire population given.
sub stdev(\@) {
	my @data = @{$_[0]};
	my $sum;my $count = $#data+1;
	foreach(@data) {
		$sum+=$_;
	}
	#Find the average value of elements
	my $avg = $sum/$count;
	$sum=0;
	
	#Find the difference between each from the average, square it
	foreach(@data) {
		$_ = $_-$avg;
		$_ = $_ * $_;
		$sum += $_;
	}
	
	#And return the square root of the average of those values, the Std Dev.
	return sqrt($sum/$count);
}

my($temp,$outfilename,$fastafn,$annotfn,%data_table,%count,$user_t,$system_t,$cuser_t,$csystem_t);

if($#ARGV+1 != 3) {
	print "Usage: $0 <CHR OUT file> <REP LIB FASTA file> <REP LIB ANNOT file>\n";
	exit;
}

$outfilename = $ARGV[0];
$fastafn = $ARGV[1];
$annotfn = $ARGV[2];

open(OUTFILE,$outfilename) or die ("Could not open $outfilename: $!.");

#Finds all families, putting them in a hash to their number of elements.
#This way, we have a unique list of repeat families.
while(<OUTFILE>) {
	my $rec = [split(' ',$_)];
	my($fam); 
	foreach $fam (keys(%count)) {
		if($fam eq $rec->[9]) {
			$count{$fam}++;
		}
	}
	if(!$count{$rec->[9]}) {
		$count{$rec->[9]}=1;
	}
}
close(OUTFILE);

open(SUMMARY,">$fastafn.summ") or die ("Could not open $fastafn.summ for write: $!");
#Print the form of the document at the top.  All doc. comments will start with a pound symbol.
$temp=localtime();
print SUMMARY "\# Time: ",$temp,"\n\n";
print SUMMARY  "# Fam\tAvgSWScore\tStdDevScore\tAvgEleLength\tStdDevLength\tNumElements\tSubs%\tDel%\tIns%\tConseqLength\n";
print SUMMARY "# Ref Fam\tConseqLength\tAnnotCount\tAnnotRatio\n";
print SUMMARY "# Annotations tabbed\n\n";

#Initiate variables for the next block and the array values for the sum data.
my(@totaldata,$key,$totalcount);
for (0..8) {$totaldata[$_]=0;}

open(OUTFILE,$outfilename) or die ("Could not open $outfilename: $!.");
open(FASTAFILE,$fastafn) or die ("Could not open $fastafn for read: $!");
open(ANNOTFILE,$annotfn) or die ("Could not open $annotfn for read: $!");

#Slurping in files to avoid access time.
my(@fastafile,@outfile,@annotfile);
while(<OUTFILE>) {
	push(@outfile,$_);
	}
close(OUTFILE);
while(<FASTAFILE>) {
	push(@fastafile,$_);
}
close(FASTAFILE);
while(<ANNOTFILE>) {
	push(@annotfile,$_);
	}
close(ANNOTFILE);

#Each key in %count is a repeat family.  Iterate through those
foreach $key (sort keys %count) {
	my(@sw_scores,@elem_lengths,$rec,$score,$numelements,$cseqlength,$inreplength,$elelength,$pctsub,$pctdel,$pctins);
	($score,$cseqlength,$elelength,$pctsub,$pctdel,$pctins,$inreplength) = (0,0,0,0,0,0,0);
	$numelements = $count{$key};
	$totalcount++;
	
	#Open the RepeatMasker .OUT file to retrieve data about the repeat
	foreach (@outfile) {
		$rec = [split(' ',$_)];
		if($rec->[9] eq $key) {
			push(@sw_scores,$score);
			push(@elem_lengths,$elelength);
			$score+=$rec->[0];
			$elelength += abs($rec->[6] - $rec->[5]);
			$inreplength += $rec->[12];
			$pctsub += $rec->[1];
			$pctdel += $rec->[2];
			$pctins += $rec->[3];
		}
	}
	
	#Open the FASTA file related to it to find the length of the sequence.
	my $lastkey = 0;
	my $fastalength=0;
	my $k;
	foreach (@fastafile) {
		if($_ =~ /^>/) {
			if($fastalength != 0 && $lastkey eq $key) {
				$cseqlength=$fastalength;last;
			}
			else {
				$fastalength=0;
				$_ =~ /(R=[\d]+)/;
				$lastkey = $1;
				next;
			}
		}
		else {
			chomp($_);
			$fastalength += length($_);
		}
	}
	# for last repeat family in fasta file
	$cseqlength=$fastalength;
	#Calculate standard deviations of the SW Score and Element lengths
	my $scorestdev = stdev(@sw_scores);
	my $lengthstdev = stdev(@elem_lengths);
	
	$score /= $numelements;$elelength /= $numelements;$pctsub /= $numelements;$pctdel /= $numelements;$pctins /= $numelements;
	
	$score = int($score);$elelength = int($elelength);$pctsub = int($pctsub);$pctdel = int($pctdel);$pctins = int($pctins);
	
		print SUMMARY "$key\t$score\t",&round2($scorestdev)."\t$elelength\t",&round2($lengthstdev),"\t$numelements\t$pctsub\t$pctdel\t$pctins\t$cseqlength\n";
	
	$totaldata[0]+=$score;$totaldata[1]+=$scorestdev;$totaldata[2]+=$elelength;$totaldata[3]+=$lengthstdev;$totaldata[4]+=$numelements;$totaldata[5]+=$pctsub;$totaldata[6]+=$pctdel;$totaldata[7]+=$pctins;$totaldata[8]+=$cseqlength;
	
	my $annotcount=0;
	
	#Retrieve annotations to print afterwards from the ANNOT file.  Tab them over except the summary first line. 
	foreach $_ (@annotfile) {
		if($_ =~ m/$key\s/) {
			if($annotcount>0) {
				print SUMMARY "\t$_";
			}
			else {
				print SUMMARY "$_";
			}
			$annotcount++;
		}
		
		elsif ($annotcount>0) {print SUMMARY "\n";last;}
	}
	
	if(!$annotcount) {print SUMMARY "\n";}
		
}

#Printing out an overall summary below.
print SUMMARY "\n# Data Summary:";
print SUMMARY  "\n# Fam\tAvgSWScore\tStdDevScore\tAvgEleLength\tStdDevLength\tNumElements\tSubs%\tDel%\tIns%\tConseqLength\n";

#The only one that doesn't need to be averaged is the number of elements.  So we exclude it..
for (0..8) {unless($_==4){$totaldata[$_] /= $totalcount;}}
my $i;
print SUMMARY "\# ";
foreach $i (@totaldata) {
	$i = &round2($i);
	print SUMMARY "$i\t";
}

print SUMMARY "\n";
close(SUMMARY);

# calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print "\# Runtime details : \n";
print "\# System time for process: $system_t\n";
print "\# User time for process: $user_t\n";

exit;