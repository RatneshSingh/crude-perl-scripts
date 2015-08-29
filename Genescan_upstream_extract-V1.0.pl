#!/usr/bin/perl -w


##################################################################################
#   This script will read through genscan output file and collect                #
#   information related to ORF startpoint and print them in a file               #
#   and then use this information to extract upstream sequences from             #
#   the supercontigs                                                             #
#   Version:1.0                                                                  #
#   Author: Ratnesh Singh  email bugs at ratnesh@hawaii.edu                      #
##################################################################################


use strict;
use warnings;

#Declair global variables
my ($header,$seq,%sequence,$len);


#open the directory containing files(current directory".").
opendir(DIR, ".");

#read filesnames (readdir(DIR)  having .txt pattern in an array
my @files = grep(/genscan/,readdir(DIR));
closedir(DIR);
print "files to be read @files\n";


#read all the sequence in memory as hash.
open (SEQ,$ARGV[0]) or die "Cant open $ARGV[0]";
$/="\n>";
while(<SEQ>){
	($header,my@seq)=split(/\n/,$_);
	$header=~s/>//;
	$seq=join("",@seq);
	$seq=~s/\s//g;
	$sequence{$header}=$seq;
	
	}
$len=$ARGV[1];
chomp($len);
#mkdir(./Upstream_ORF);
open(OUT2,">Upstrem_ORFs.txt");
#processing files and gathering information to cut upstream region.
$/="\n";
foreach my $file (@files) {

#null sequence name before reading new file
	my $sequence=();
	open(FILE,"$file") or die "Can't open file $file";
#	open(OUT2,">./Upstream_ORF/$")
#reading each file one a a time and parsing through for ORF start site.
	while (<FILE>) {
		if(/^\s*$/){next;}

		if(/^Sequence\s+(\w+)/){$sequence=$1;print "\nsequence name is:$sequence\n";next;}
		if(/Init/){ 
			$_=~s/(\s)+/ /g;
			$_=~s/^\s+//g;
#			$_=~s/\+/PLUS/g;
#			$_=~s/\-/MINUS/g;
			my	($GnEx,$Type,$Strand,$Begin,$End,$Len,$Fr,$Ph,$I_Ac,$Do_T,$CodRg,$P,$Tscr)=split(/ /,$_);	
	
			my $u_str= extract_sequence($sequence,$Begin,$Strand,$len);
			if($u_str){
				
			print OUT2 (">$sequence".'_'."$Strand".'_'.'ORF'.'_'."$GnEx".'('."$Begin".')',"\n$u_str\n");
#			print "\n\n$_\n";
#			print OUT"Sequence:$sequence\tGenEx:$GnEx\tType:$Type\tStrand:$Strand\tBegin:$Begin\tEnd:$End\tLen:$Len\tFr:$Fr\tPh:$Ph\tI_Ac:$I_Ac\tDo_T:$Do_T\tCodRg:$CodRg\tP:$P\tTscr:$Tscr\n";
			}
			next;
			}
	
		else{
			next;
			}
  		}
  	}
#close(OUT);
close(OUT2);
exit;
	
	
	
sub extract_sequence{
	
	my ($seq_header,$begin,$strand,$len)=@_;
	$strand=~s/\-/MINUS/;
	$strand=~s/\+/PLUS/;
	my$upstream_start=0;
#calculation for upstream regions start and end.
#
# (+ strand)--(Upstream start = $begin-1501)---------------------------ORF-start-->(Upstream end=$begin)-------------------------------------------------------------
# (- strand)--------------------------------------------------------<--ORF-start-(upstream start=$begin)-----------(Upstream end=$begin+1501)------------------------
#	
	if($strand eq 'PLUS'){
		$upstream_start=$begin-$len-1;
#		$upstream_end=$begin;
	}
	if($strand eq 'MINUS'){
		$upstream_start=$begin;
#		$upstream_end=$begin+1501;

	}
	my $seq_upstream=substr($sequence{$seq_header},$upstream_start,$len);
	return($seq_upstream);
	
	
	
	
}