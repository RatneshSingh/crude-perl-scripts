#!/usr/bin/perl -w
use strict;

#this script will concatenate trailing fragments of nucleotide under one header into one line.
#In other words it will remove "\n" character from the graments of sequence under one header.

my($inputfile,$header,@sequence,$sequence,$outputfile,$patternfile,$pattern,$count,@pattern,$ta);

print " \nWelcome!!! This script will Concatenate several fragment of one sequence under one Header in to One line\n";

#opening inputfile containing sequences in fasta format.
if ($ARGV[0]){$inputfile=$ARGV[0];}else{print "Enter input file containing sequences\n"; $inputfile=<STDIN>;}
chomp ($inputfile);
open INPUT,"<$inputfile" or die "Cannot open $inputfile.....\n\n";

#opening outputfile
if ($ARGV[1]){$outputfile=$ARGV[1];}else{print "Enter output file name\n"; $outputfile=<STDIN>;}
chomp($outputfile);
open OUT,">$outputfile" or die "cannot create $outputfile.....\n\n";

#opening file containing TA or pattern to search seperated by new line.
if ($ARGV[2]){$patternfile=$ARGV[2];}else{print "Enter pattern file containing pattern seperated by newline\n"; $patternfile=<STDIN>;}
chomp ($patternfile);
open PATTERN,"<$patternfile" or die "Cannot open $inputfile.....\n\n";


#make an array from the input patterns t search
while(<PATTERN>){
$pattern=$_ if($_ ne /\s*$/);
chomp($_); 
push (@pattern,$pattern);
}
#$count=@pattern; #just to check if array is formed or not
#print "\n@pattern\n";
$/="\n>";

foreach $ta(@pattern){
pickseq($ta);
if($header=~/^>/){print OUT"$header\n$sequence\n\n";}
else{print OUT">$header\n$sequence\n\n";
}
}

#******************************************************************
sub pickseq {
my ($pattern)=@_;

while(<INPUT>){
chomp;

if(/$pattern/){
my($header,@sequence)=split("\n",$_);
my ($sequence)=join("",@sequence);
$sequence=~ s/\s//;
return ($header,$sequence);
next;}
}
}
}




