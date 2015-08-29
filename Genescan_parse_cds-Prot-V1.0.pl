#!/usr/bin/perl -w
#This script will read through genscan output file and collect
#coding sequences from it into different file.keep this script
#in the folder containing genscan output and run from their.
#Version:1.0
#Author: Ratnesh Singh


use strict;

#open the directory containing files(current directory".").
opendir(DIR, ".");




#read filesnames (readdir(DIR)  having .txt pattern in an array
my @files = grep(/genscan/,readdir(DIR));
closedir(DIR);
#print "files to be read @files\n";
open(PEPT,">papaya_pept_genscan.fasta");
open(CDS,">papaya_CDS_genscan.fasta");



#processing files and printing information one after one.
foreach my $file (@files) {
	my (@file)=();
	open(FILE,"$file") or die "Can't open file $file";

	$/="\n>";
	while(<FILE>){
			my $line=$_;
			if($line=~/>/){
				$line=~s/>//g;
				$line=~s/\s*$//g;
				push(@file,$line);
			}
	}
	# remove the first block from the array as it contains table of coordinates
	shift @file;

	# process rest of the array.
	my $size=@file; print "number of coding sequences in $file ::  $size\n";

	for(my $i=0;$i<$size;$i++){
								my($header,@sequence)=split(/\n/,$file[$i]);
								my$sequence= join("",@sequence);
								$sequence=~ s/\s//g;

								my $bypass=$header;
								$bypass=~ s/\|/\t/g;
								my $pattern1="GENSCAN_predicted_peptide";
								my $pattern2="GENSCAN_predicted_CDS";

								print PEPT ">$header\n$sequence\n" if ($bypass=~ /$pattern1/i);


								print CDS ">$header\n$sequence\n" if ($bypass=~ /$pattern2/i);

	}
}

exit;