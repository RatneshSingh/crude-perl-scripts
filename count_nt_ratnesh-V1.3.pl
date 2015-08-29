#!/usr/bin/perl -w

#This scripts calculates the DNA sequence length and output as table or adds in seq name
#
# usage: perl ScriptName -f SeqFile
#
# other options (optional):
# -l\t ThreshHold length\tFilter out sequences of length.
# -t\t TRUE/YES \t Prints Headers and selengths as a table [Default].
# -h\t TRUE/YES \t Prints seqlengths attached in header.
# -g\t TRUE/YES \t calculates and prints GC percents of seqs in table.
# -c\t TRUE/YES \t trim the header till first white space.
#	updated on : 03/02/2011
#	created by Ratnesh Singh

 use warnings;
 use strict;
 use Getopt::Std;

 our($opt_f,$opt_l,$opt_t,$opt_h,$opt_g,$opt_c,$opt_d);
 getopt('flthgcd');

 my $help="This scripts calculates the DNA sequence length and output as table or adds in seq name

 usage: perl ScriptName -f SeqFile

 other options:
 -l\t FilterThreshHold\tFilter out sequences of length.
 -t\t TRUE/YES \t Prints Headers and selengths as a table [default].
 -h\t TRUE/YES \t Prints seqlengths attached in header.
 -d\t Char     \t Use this delimiter to join SeqLength to Header [_]
 -g\t TRUE/YES \t calculates and prints GC percents of seqs in table.
 -c\t TRUE/YES \t trim the header till first white space.

 ";

print "$help\n\nProcessing Sequences. Please wait.....\n";
#read file for the sequences and load it in hash
if(!defined $opt_f){print "Please provide sequence file after '-f' flag \n\n$help\n"; die;}
else{open(FASTAFILE,$opt_f) or die"Cannot find the file $opt_f \n";}

if(defined $opt_l||defined $opt_t||defined $opt_h||defined $opt_g){}
else{$opt_t='TRUE';}

#open out files#


if(defined $opt_l){
		open(OUT2,">largerThan.$opt_l.$opt_f");
		open(OUT3,">smallerThan.$opt_l.$opt_f");
}

if(defined $opt_t){open(OUT4,">$opt_f._SeqLengths.table_");}
if(defined $opt_h){open(OUT5,">SeqLengthsinName_.$opt_f");}
if(defined $opt_g){open(OUT6,">$opt_f._GCcontent.table");}

if(!$opt_d){$opt_d="_";}
chomp($opt_d);

my ($header,$seq,$length,$A,$C,$T,$G,$N,$GC);
#define record seperator ($/) as "\>"

$/="\n>";

#seperate header and sequence in a hash called $header and @seqarray
my$removed_count=0;
my $counter=0;
while (<FASTAFILE>){
	chomp;
   #split record into header and sequence and feed to $header and @seq
	($header,my@seq)= split (/\n/,$_);
	$header=~s/>//;
	if(defined $opt_c){($header,my@rest)=split(/\s+/,$header); chomp($header);}
	$seq=join ("",@seq); #join elements of @seq array to remove any spaces.
	$seq=~ tr/a-z/A-Z/;
	$A=$seq=~tr/A//;  #count the number of A T G C in $seq
	$G=$seq=~tr/G//;
	$T=$seq=~tr/T//;
	$C=$seq=~tr/C//;
	$N=$seq=~tr/KMRYSWVBHDNX//;

	$length=$A+$G+$C+$T+$N;
	chomp ($length,$header,$seq);
	$GC=($G+$C)*100/$length if $length >0;





	if(defined $opt_l){

		if($length>$opt_l){

			print OUT2 ">$header\n$seq\n";
    	}
		else {
			#print "removing\n$header\n$seq\n";
			print OUT3">$header\n$seq\n";$removed_count++;
		}

	}


	# Create table of headers and sequence length in two colums;
	if(defined $opt_t){
		print OUT4 "$header\t$length\n";
	}

	if(defined $opt_h){
		print OUT5 ">$header$opt_d$length\n$seq\n";
	}

	if(defined $opt_g){
		print OUT6 "$header\t$GC\n";

	}

$counter++;
}

if(defined $opt_l){
print "\ntotal number of sequences removed as they are smaller than $opt_l nt:\t$removed_count\n\n";
}

print "Total number of sequences processed:$counter\n\n.";

close FASTAFILE;
close OUT2;
close OUT3;
close OUT4;
close OUT5;
exit;
