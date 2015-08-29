#!/usr/bin/perl -w

              							#read file for the sequences and load it in hash

open(FASTAFILE,$ARGV[0]) or die"Cannot find the file \n";
open(OUT,">nuc_number.txt");
								#define record seperator ($/) as "\>"

$/="\n>";

								#seperate header and sequence in a hash called $header and @seqarray

while (<FASTAFILE>)
{
chomp;
							        #split record into header and sequence and feed to $header and @seq
($header,@seq)= split (/\n/,$_);
								#join elements of @seq array to remove any spaces.
$seq=join ("",@seq);
	
$P=$seq=~tr/P//;                                                #count the number of A T G C in $seq
$length=length$seq;
$perP=($P/$length)*100;
$NonP= $length-$P;

#calculate %GC content in $seq
		
print "In sequence\n>$header\nnumber of P:$P,  number of NonP:$NonP, Sequence length :$length, percentP : $perP % \n";
			   


}

close FASTAFILE;
close OUT;		
exit;		
