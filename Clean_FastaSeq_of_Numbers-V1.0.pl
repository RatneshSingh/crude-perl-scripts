#!/usr/bin/perl -w
use strict;

my($inputfile,$header,@sequence,$sequence,$outputfile,$patternfile,$pattern,$count,@pattern,$ta,%seq_hash);

print " \n****Welcome!!! This script will read pattern to search from the given file and pick out sequence containing pattern from another input file. can be used AS MANUAL MODE***  \n\n\nusage : perl script_name sequence_file outputfile\n\n";

#opening inputfile containing sequences in fasta format.
if ($ARGV[0]){$inputfile=$ARGV[0];}else{print "Enter input file containing sequences\n"; $inputfile=<STDIN>;}
chomp ($inputfile);
open INPUT,"<$inputfile" or die "Cannot open $inputfile.....\n\n";


#opening outputfile
if ($ARGV[1]){$outputfile=$ARGV[1];}else{print "Enter output file name\n"; $outputfile=<STDIN>;}
chomp($outputfile);
open OUT,">$outputfile" or die "cannot create $outputfile.....\n\n";


#Read database sequence in to a hash %seq_hash{}, remove ">" from header. remove any white spaces and newline from sequence.
print "reading Sequences from input file.....Plz wait...\n";
$/="\n>";
$count=0;

while(<INPUT>){#
    chomp;
    ($header,@sequence)=split("\n",$_);
    $header=~s/>//;
    $sequence= join("",@sequence);
    $sequence=~ s/\s//g;
    $sequence=~s/\n//g;
    $sequence=~s/\d//g;


	print OUT">$header\n$sequence\n";
	$count++;
}#

#print "Done....\n$count sequences cleaned and written in the file $outputfile\n\n";
#@seq_count=();
