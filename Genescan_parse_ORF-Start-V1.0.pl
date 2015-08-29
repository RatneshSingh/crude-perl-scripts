#!/usr/bin/perl -w
#This script will read through genscan output file and collect
#information related to ORF startpoint and print them in a file
#which can be used by any other script to cut upstream seqs 
#from the bigger sequences
#Version:1.0
#Author: Ratnesh Singh
use strict;

#open the directory containing files(current directory".").
opendir(DIR, ".");
#read filesnames (readdir(DIR)  having .txt pattern in an array
my @files = grep(/genscan/,readdir(DIR));
closedir(DIR);
print "files to be read @files\n";
open(OUT,">output.txt");
#processing files and printing information.
foreach my $file (@files) {
my $sequence=();
open(FILE,"$file") or die "Can't open file $file";

#reading each file one a a time and parsing through.
while (<FILE>) {
	if(/^\s*$/){next}

	if(/^Sequence\s+(\w+)/){$sequence=$1;print "\nsequence name is:$sequence\n";next;}
	if(/Init/){ 
		$_=~s/(\s)+/ /g;
		$_=~s/^\s+//g;
#		$_=~s/\+/PLUS/g;
#		$_=~s/\-/MINUS/g;
	my	($GnEx,$Type,$Strand,$Begin,$End,$Len,$Fr,$Ph,$I_Ac,$Do_T,$CodRg,$P,$Tscr)=split(/ /,$_);	
			
#print "\n\n$_\n";
	print OUT"Sequence:$sequence\tGenEx:$GnEx\tType:$Type\tStrand:$Strand\tBegin:$Begin\tEnd:$End\tLen:$Len\tFr:$Fr\tPh:$Ph\tI_Ac:$I_Ac\tDo_T:$Do_T\tCodRg:$CodRg\tP:$P\tTscr:$Tscr\n";
	next;}
	
	else{next;}
  }
  }
close(OUT);
	exit;