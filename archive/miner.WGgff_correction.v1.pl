#!/usr/bin/perl -w
# MGEL
# Surya Saha 6/1/09
# reading cmd line input .gff file with whole genome coordinates for RScout families,
# and writing a new WG GFF output file with correct chromosome numbers and coordinates
# Output is sorted on the start position



use strict;
use warnings;
use POSIX;

unless (@ARGV == 1){
	print "USAGE: $0 <input WG .gff file>\n";
	exit;
}

# find the correct coordinate on the correct chromosome
sub resolv_coord{
	my (@temp_coord,$ii,$jj,$vcoord,@chrs,$rcoord,$rchr,$chrlen);
	$vcoord=$_[0];
	# hard coding lengths of 12 rice chrs
	$chrs[0]=43596771; $chrs[1]=35925388; $chrs[2]=36345490;
	$chrs[3]=35244269; $chrs[4]=29874162; $chrs[5]=31246789;
	$chrs[6]=29688601; $chrs[7]=28309179; $chrs[8]=23011239;
	$chrs[9]=22876596; $chrs[10]=28462103; $chrs[11]=27497160;
	
	for($ii=0;$ii<12;$ii++){
		$chrlen=0;
		# get the actual chr lengths
		for ($jj=0;$jj<=$ii;$jj++){
			$chrlen=$chrlen+$chrs[$jj];
		}
		# add in 15k padding and find chr, get real coord
		if ($vcoord <= ($chrlen + 15000*$ii)) {
			$rchr=$ii+1;# real chr
			# get real coord
			$chrlen=$chrlen-$chrs[$ii];
			$rcoord=$vcoord-((15000*$ii)+$chrlen);
			last;
		}
	}
	
	$temp_coord[0]=$rchr;
	$temp_coord[1]=$rcoord;
	
	return @temp_coord;
}


#my ($ifname,$len,$rec,@temp,%temphash,@counts,@table,@famnames,@outarr,$temp,$ctr,$i,$j);
my ($ifname,$rec,@table,@temp,$start,$end,$chr);

$ifname=$ARGV[0];
chomp $ifname;

unless(open(INFILEDATA,$ifname)){print "not able to open ".$ifname."\n\n";exit;}

unless(open(OUTFILEDATA,">updated.$ifname")){print "not able to open updated.".$ifname." \n\n";exit;}
print OUTFILEDATA "[Repeat]\nfeature = repeat_unit:ProxMiner\nglyph = segments\n";

#@table
# chr01-12	RptSct-1.0.2	repeat_region	2484	2552	229	+	.	R=4368
# chr01-12	RptSct-1.0.2	repeat_region	2493	2565	251	-	.	R=11696
# chr01-12	RptSct-1.0.2	repeat_region	2510	2617	707	+	.	R=4606
#  0              1               2           3    4    5   6    7     8  

#slurping in the whole report file
#adding correct chr, start and end positions 
#and writing it out
while($rec=<INFILEDATA>){
	if($rec =~ /#/){next;}
	if(length ($rec) < 2){next;}
	
	push (@table,[split(' ',$rec)]);
	
	@temp=&resolv_coord($table[0][3]);
	$start=$temp[1];
	@temp=&resolv_coord($table[0][4]);
	$end=$temp[1];
	$chr=$temp[0];
	
	print OUTFILEDATA "chr$chr\t$table[0][1]\t$table[0][2]\t$start\t$end\t$table[0][5]\t$table[0][6]\t$table[0][7]\tFam $table[0][8]\n";
	
	
	@table=();
}

close (INFILEDATA);
close (OUTFILEDATA);

exit;
