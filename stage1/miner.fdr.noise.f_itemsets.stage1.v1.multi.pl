#!/usr/bin/perl -w
# MGEL
# Surya Saha 12/25/08
# FOR NOISE DATA
# reading multiple noise f_itemset files and writing out FDR statistics
# writing out number of rules with confidence above certain thresholds

my $ver=1.0;
use strict;
use warnings;
use POSIX;

# fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, Occurence, Category, Avg. dist., Std. dev.
# R=424	72	120	R=437	67	346	2	U	4512	2699.49
# 0	1	2	3	4	5	6	7	8	9

unless (@ARGV == 2){
	print "USAGE: $0 <core input noise f_itemsets file name> <replicates>\n";
	print "NOTE: Make sure core name is only extended by .<#>.f_itemsets.stg1.tab for all replicates\n";
	exit;
}

my ($icorefname,$rec,$reps,@temp,%rel_conf_ctrs,%rel_count_ctrs,%noise_hash,@values,$conf,
$ctr,$i,$j,$k,$l,$key,$value,$user_t,$system_t,$cuser_t,$csystem_t);

$icorefname=$ARGV[0];
chomp $icorefname;
$reps=$ARGV[1];
chomp $reps;

# declaring all hashes for replicates
# for($i=1;$i<=$reps;$i++){ %{noise_hash_."$i"}={1};}

unless(open(OUTFILEDATA,">$icorefname.FDR.tab")){print "not able to open ".$icorefname."FDR.tab\n\n";exit;}

sub round2{
	my ($num);
	$num=$_[0];
	$num=$num*100;
	$num=int($num);
	$num=$num/100;
	return $num;
}

# fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, Occurence, Category, Avg. dist., Std. dev.
# R=424	72	120	R=437	67	346	2	U	4512	2699.49
# 0		1	2	3		4	5	6	7	8	9

# slurping in the all noise itemset files
# put the noise data into a hash
for($i=1;$i<=$reps;$i++){
	unless(open(INFILENOISEDATA,$icorefname.".".$i.".f_itemsets.stg1.tab")){print "not able to open ".$icorefname.".".$i.".f_itemsets.stg1.tab"."\n\n";exit;}
	$ctr=0;
	
	while($rec=<INFILENOISEDATA>){
		if($rec =~ /#/){next;}
		if(length ($rec) < 10){next;}#for avoiding last line
		@temp = split(' ',$rec);
		# for families with elements within a range
		# COMMENT OUT IF FOR ALL FAMS
# 		if ((($temp[1]>=50) && ($temp[1]<=100)) && (($temp[4]>=50) && ($temp[4]<=100))) 
# 		if (($temp[1]<=50) && ($temp[4]<=50))
# 		if (($temp[1]>=100) && ($temp[4]>=100))
		{
			# calculating conf
			if($temp[1] > $temp[4]){
				$values[0]= &round2($temp[6]/$temp[4]);
			}
			elsif($temp[4] > $temp[1]){
				$values[0]= &round2($temp[6]/$temp[1]);
			}
			else{#both are equal, same family relation??
				$values[0]= &round2($temp[6]/$temp[1]);
			}
			$values[1]=$temp[6];
			$noise_hash{"$i $ctr"}=[@values];
			@values=();
			$ctr++;
		}
	}
	close (INFILENOISEDATA);
}

# initialize counters
for ($i=0.05;$i<=1.05;$i+=0.05){
	$rel_conf_ctrs{"$i"}=0;
	$rel_count_ctrs{"$i"}=0;
}


$i=localtime();
# calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;

print STDERR "\# Runtime details after preparing all hashes:\n";
print STDERR "\# System time for process: $system_t\n";
print STDERR "\# User time for process: $user_t\n\n";

#filling up the @merged_table array the merged file
# @real_table
# fam1, fam1-count, fam1-avglen, fam2, fam2-count, fam2-avglen, Occurence, Category, Avg. dist., Std. dev.
# R=424	72	120	R=437	67	346	2	+	4512	2699.49
# 0	1	2	3	4	5	6	7	8	9

# updating the counters
while(($key, $value)= each %rel_conf_ctrs){
	
	@values= values %noise_hash;
	foreach $j (@values){
		if(@$j[0] >= $key){
			$rel_conf_ctrs{$key}++;#increment ctr if valid
			$rel_count_ctrs{$key}=$rel_count_ctrs{$key}+@$j[1];#increment ctr if valid
		}
	}
}

# print the connected components
$i=localtime();
# calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;

print OUTFILEDATA "\# Version: $ver\n";
print OUTFILEDATA "\# Time: $i\n\n";
print OUTFILEDATA "\n**********************************************************************\n\n";
print OUTFILEDATA "\# Core noise f_itemsets file used : $icorefname\n";
print OUTFILEDATA "\# Replicates used : $reps\n";

# getting the ctrs in a sorted array
print OUTFILEDATA "Confidence\tNumber \tAvg \n";
print OUTFILEDATA "threshold\t\tof rules\tcount\n";
foreach $key (sort keys %rel_conf_ctrs){
	if ($rel_conf_ctrs{$key}>0 && $rel_conf_ctrs{$key}<=$reps){
		print OUTFILEDATA $key,"\t\t\t",$rel_conf_ctrs{$key},'*';
	}
	else{
		print OUTFILEDATA $key,"\t\t\t",int($rel_conf_ctrs{$key}/$reps);
	}
	if ($rel_conf_ctrs{$key} > 0){
# 		print OUTFILEDATA "\t\t",int(int($rel_count_ctrs{$key}/$rel_conf_ctrs{$key})/$reps);
		print OUTFILEDATA "\t\t",int($rel_count_ctrs{$key}/$rel_conf_ctrs{$key});
	}
	print OUTFILEDATA "\n";
}

print OUTFILEDATA "\n\n**********************************************************************\n\n";
# calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print OUTFILEDATA "\# Runtime details : \n";
print OUTFILEDATA "\# System time for process: $system_t\n";
print OUTFILEDATA "\# User time for process: $user_t\n";

close (OUTFILEDATA);


# calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print "\# Runtime details : \n";
print "\# System time for process: $system_t\n";
print "\# User time for process: $user_t\n";


exit;
