#!/usr/bin/perl -w
#This script is written to parse through ciselements out put file sent 
#by jimmy to mingli. it reads through the file line by line and break 
#apart the elements into the respective columns and make it is easy to
#sort by excel.
#Version:1.0
#Author: Ratnesh Singh
use strict;

my $file= $ARGV[0];
open FILE,"$file" or die ("Cant open $file");
open OUT,">$file.parsed" or die ("cant open output file");
print OUT"supercontig\tstrand\torf\tUStr_start\tUStr_end\tcis_start\tciselement\n";
while (<FILE>){
	s/^\s//;
	s/\s+/ /g;
#	print $_;
	my ($super,$at,$cis_start,@ciselement)=();
	   ($super,$at,$cis_start,@ciselement)=split(/ /,$_);
	my  $ciselement=join("_",@ciselement);
		$cis_start=~s/\[//;
		$cis_start=~s/\]//;
		
#	print "super: $super \t start: $cis_start \t cis: $ciselement \n";
	my(@super)=split(/_/,$super);
	my $supercontig=$super[0].'_'.$super[1];
	my $strand= 'MINUS' if($super[2] eq "-");
	   $strand='PLUS' if($super[2] eq "+");
#	print $supercontig, "\t", $strand,"\n";	
		$super[3]=~s/\[//;
		$super[3]=~s/\]//;
	my ($UStr_start,$UStr_end)= split(/,/,$super[3]);

#	print "start: $UStr_start end: $UStr_end \n";
	my $orf=$super[4];
#	print "ORF is :: $orf \n"
print "$supercontig\t$strand\t$orf\t$UStr_start\t$UStr_end\t$cis_start\t$ciselement\n";
print OUT"$supercontig\t$strand\t$orf\t$UStr_start\t$UStr_end\t$cis_start\t$ciselement\n";

}
exit;