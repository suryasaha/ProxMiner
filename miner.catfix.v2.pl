#!/usr/bin/perl -w
#
#Jonathan Harper 06/11/2008
#
#.CAT Name Fix
#
#This file reads in a .cat file and a FASTA file in order to restore the .CAT
#file's original name from the source (FASTA file).
#Usage: perl catfix.pl <path of .CAT file> <path of FASTA file>
#
#v2: more comments added, basics fixed
#  : total rewrite, some optimization, split up into parts for easier reading, added meta-data

use strict;
use warnings;
use Tie::File;

# declare all variables here
my($catfile,$fastafile,$outfile,@cat_table,@meta_data,$user_t,$system_t,$cuser_t,$csystem_t,$counter,$i);

unless(@ARGV == 2) {
	print "Usage: perl catfix.pl <path of .CAT file> <path of FASTA file>\n";
	exit;
}

$catfile = $ARGV[0];
$fastafile = $ARGV[1];
chomp $catfile;chomp $fastafile;

open(CATFILE,$catfile) or die ("Could not open $catfile for read: $!");

$counter = 0;

#Slurping in the CAT file here.  There is a set of 'meta-data' at the end of every CAT file
#which needs to be preserved.  The first lines of it start with '#' and are easy to track but
#the others I skip over and add to the end of the next one by using a counter since the number
#of them is constant between CAT files.
while(<CATFILE>) {
	if($_ !~ /^#/ && $counter==0) {
		if($_ !~ /RepeatMasker/) {
			push(@cat_table,[split(' ',$_)]);
		}
		
		else {
			$counter=3;
			push(@meta_data,$_);
			$counter--;
			next;
		}
	}
	else {
		if($counter>0) {$counter--;}	
		push(@meta_data,$_);
	}
}

#Sorting by repeat name ASCII-betically.  R=1,R=104,R=111, et cetera.
@cat_table = sort {$a->[4] cmp $b->[4]} @cat_table;

close(CATFILE);

open(FASTAFILE,$fastafile) or die ("Could not open $fastafile for read: $!");

#This loops through the FASTA file, matching each header line (which starts with a >)
#with an entry in the CAT table.  All whitespace must be removed (preserved by replacing
#with dash ("-") characters) because of the format of the CAT file.
$counter = 0;
while(<FASTAFILE>) {
	if($_ =~ /^>/) {
		$_ =~ s/\s+$//g; # remove trailing whitespace
		$_ =~ s/^>//g; # remove >
		$_ =~ s/\s/-/g; # replace ' ' wt -
		my $rec = $_;
		foreach $i (@cat_table) {
			my $match;
			if($i->[8] eq "C") {$match=$i->[9];} else {$match=$i->[8];}
			if($rec =~ m/$match/) {
				if($i->[8] eq "C") {$i->[9]=$rec;} else {$i->[8]=$rec;}
			}
		}
		$counter++;
		if ($counter%50 == 0) {print STDERR '.';}
	}
}
print STDERR "\n";

open(CATFILE,">$catfile") or die ("Could not open $catfile for read: $!");
#Finally, this prints the cat_table for all fields and then adds the meta-data on the end.
foreach $i (@cat_table) {
	my $string = "";
	$counter = 0;
	do{
		$string = $string.$i->[$counter].' ';
		$counter++;
	} until (!$i->[$counter]);
	$string = $string."\n";
	print CATFILE $string;
}


foreach $i (@meta_data) {
	print CATFILE $i;
}

close(CATFILE);

#calculating time taken
($user_t,$system_t,$cuser_t,$csystem_t) = times;
print STDERR "\# Runtime details : \n";
print STDERR "\# System time for process: $system_t\n";
print STDERR "\# User time for process: $user_t\n";

exit;
