#!/usr/bin/perl -w
use strict;


my($inputfile,$header,@sequence,$sequence,$outputfile,$patternfile,$pattern,$count,@pattern,$ta);

print " \n****Welcome!!! This script will read pattern to search from the given file and pick out sequence containing pattern from another input file***  \n";



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




#make an array from the input patterns search
while(<PATTERN>){
chomp;
$pattern=$_ if($_ ne /\s*$/);
chomp($_); 

#my($a,$b,$c,$d)= split(/\|/,$_);
#chomp ($d);
$pattern=$_;
push (@pattern,$pattern);
}
#$count=@pattern; #just to check if array is formed or not
#print "\n@pattern\n";
$/="HMMER2.0  [2.3.2]";



#read fasta file and search for the pattern...

while(<INPUT>){#
chomp;

#print "This is line of input file $_\n this is pattern $ta\n";

foreach $ta(@pattern){##
chomp($ta);

#print "this is $ta after chomping $ta";

if($_=~ /$ta/){###
#print "Pattern:$ta found";
print OUT">$_ \n";
#($header,$sequence)=split("\n",$_);
#$sequence=join("",@sequence);# just to save time. remove # sign to return back to normal.
#$sequence=~ s/\s//; # remove # to return normal curative script. 

#print "\nThis is for:$ta \nThis  is header\n$header\nThis is sequence \n$sequence";

#if($header=~/^>/){print OUT"$header\n$sequence\n\n";}
#else{print OUT">$header\n$sequence\n\n";}
}###
#else {print "$ta not found";}
}##
}#
