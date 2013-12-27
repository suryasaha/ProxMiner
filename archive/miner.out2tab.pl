#!/usr/bin/perl -w
# MGEL
# Surya Saha 3/15/07
# reading cmd line input RM .out file,
# and writing hit details to tabbed output file
# which is sorted on the start position



use strict;
use warnings;

unless (@ARGV == 1){
	print "USAGE: $0 <input file> \n";
	exit;
}


my ($ifname,$div_lim,$del_lim,$ins_lim,$rec,@temp,@table,@famnames,$ctr,$i,$j);

$ifname=$ARGV[0];
chomp $ifname;
unless(open(INFILEDATA,$ifname)){print "not able to open ".$ifname."\n\n";exit;}

unless(open(OUTFILEDATA,">$ifname.tab")){print "not able to open ".$ifname.".tab \n\n";exit;}

#slurping in the whole report file
while($rec=<INFILEDATA>){
	if($rec =~ /#/){next;}
	if(length ($rec) < 2){next;}
	push @table, [split(' ',$rec)];
}


#for $i (0 .. $#table){
#	for $j (0 .. $#{$table[$i]}){
#		print "record $i $j is $table[$i][$j]\n";
#	}
#}

#sorting @table on start position
@temp = sort {$a->[5] <=> $b->[5]} @table;
@table=@temp;

print OUTFILEDATA "\#SW\tsub%\tdel%\tins%\tstart\tend\tfam\tq_st\tq_nd\n";
#printing the table array in tabbed format
for $i (0 .. $#table){
	print OUTFILEDATA $table[$i][0],"\t",$table[$i][1],"\t",$table[$i][2],"\t",
	$table[$i][3],"\t",$table[$i][5],"\t",$table[$i][6],"\t",$table[$i][9],"\t",
	$table[$i][12],"\t",$table[$i][13],"\t";
	print OUTFILEDATA "\n";
}

close (INFILEDATA);
close (OUTFILEDATA);

exit;
