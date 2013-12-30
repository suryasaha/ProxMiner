#!/usr/bin/perl
# MGEL
# Surya Saha
# 06/26/08

# Reading the tbalstn report (conserved domain mfasta sequences against rice chr 12) 
# Reporting all hit in GFF3 format to feed to the spatial relation
# finding pipeline.

# # <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
# chr12	RptSct-1.0.2	repeat_unit	8936	9225	1935	-	.	R=291
# chr12	RptSct-1.0.2	repeat_unit	9987	10255	1009	+	.	R=39
# chr12	RptSct-1.0.2	repeat_unit	10808	10865	226	+	.	R=112
# chr12	RptSct-1.0.2	repeat_unit	10836	10892	241	-	.	R=112
# chr12	RptSct-1.0.2	repeat_unit	10880	11115	239	+	.	R=401

use lib "/home/ssaha/tools/bioperl_modules";
use strict;
use warnings;
use Bio::Perl;
use Bio::SearchIO;


eval {require Bio::SearchIO; };
if ( $@ ) {
print STDERR "Cannot find Bio::SearchIO\n";
print STDERR "You will need the bioperl-run library to parse a tblastn report\n";
return 0;
}

unless (@ARGV == 5){
	print "USAGE: $0 <tblastn report file>  <min match length(bp)> <report type(blastn,tblastn..) <seqname> <feature>\n";
	exit;
}

my ($seqname,$feature,$ver,$i,$j,$ipfname,$minlen,@params,$factory,$input,$blast_report,$rtype,$sbjct,$result,$hsp);

$ver=1.0;

$ipfname=$ARGV[0];
chomp $ipfname;
unless(open(OUTFILE,">$ipfname.gff")){print "not able to open ".$ipfname.".gff\n\n";exit;}
$minlen=$ARGV[1];
chomp $minlen;
$rtype=$ARGV[2];
chomp $rtype;
$seqname=$ARGV[3];
chomp $seqname;
$feature=$ARGV[4];
chomp $feature;

print OUTFILE "\# Version: $ver\n";
$i=localtime();
print OUTFILE "\# Time: $i\n\n";
print OUTFILE "\# $rtype report file: $ipfname\n";
print OUTFILE "\# Minimum length of DNA hit: $minlen\n\n";
print OUTFILE "\# <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]\n";

$blast_report=new Bio::SearchIO(-format => 'blast', -file => $ipfname, -report_type => $rtype);


while($result= $blast_report->next_result){# result is Bio::Search::Result::ResultI
	
	print "Query:",$result->query_name(); #print STDERR "Query:",$result->query_name();
	print "\tLength:",$result->query_length(); #print STDERR "\tLength:",$result->query_length();

	while ($sbjct = $result->next_hit){#sbjct is Bio::Search::Hit::HitI
		print "\tAln length is ", $sbjct->length_aln(),'(AA), ',$sbjct->length_aln()*3,"(bp)\n";
		#print STDERR "\tAln length is ", $sbjct->length_aln(),'(AA), ',$sbjct->length_aln()*3,"(bp)\n";
		print "\tHSP lengths(bp): "; #print STDERR "\tHSP lengths(bp): ";
		if ($sbjct->num_hsps > 0){	
			while ($hsp=$sbjct->next_hsp){# $hsp is Bio::Search::HSP::HSPI
				# minlen should be a % of the query or hit length
				if ($hsp->length() >= $minlen){
					print $hsp->length(),' '; #print STDERR $hsp->length(),' ';
					if ($hsp->length() > ($sbjct->length_aln()*3)){ 
						print "\nERR: HSP is longer than alignment length.\n";
						#print STDERR "\nERR: HSP is longer than alignment length.\n";
					}
					# print to gff file
					print OUTFILE "$seqname\t$ipfname\t$feature\t",$hsp->start('hit')
					,"\t",$hsp->end('hit'),"\t";
					
					if ($hsp->strand('hit') == 1) {print OUTFILE '+';}
					elsif ($hsp->strand('hit') == -1) {print OUTFILE '-';}
					
					print OUTFILE "\t",$hsp->score,"\t",$sbjct->frame(),"\t",$result->query_name(),"\n";
				}
			}
		}
		print "\n"; #print STDERR "\n";
	}
}


close (OUTFILE);

exit
