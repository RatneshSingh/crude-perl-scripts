#!/usr/bin/perl -w
use strict;
use Getopt::Std;

our($opt_g,$opt_l,$seq_length,$GC_value,$seq_number,$opt_n,$opt_v);
getopt('glnv');


my$usage="\nUsage:\nperl script options...

-l	length of random sequence[1000]
-g	GC content of random sequence[0.5]
-n	number of random sequence requested[1]
-v	verbose. Show seq on screen also
\n\n";

print "$usage";
# set GC value and sequence length
if(defined $opt_l){$seq_length=$opt_l;} else {$seq_length=1000;}
if(defined $opt_g){$GC_value=$opt_g;} else {$GC_value=0.5;}
if(defined $opt_n){$seq_number=$opt_n;} else {$seq_number=1;}

open OUT,">Random_sequence.fasta";

for(my$j=1;$j<=$seq_number;$j++){
	# generate random sequence
	
	my$random_seed= int(rand(100))^time ^ $$;
	srand($random_seed);
	my$newseq=();
	for (my$i=0;$i<=$seq_length;$i++) {
		my$letter=();
		my$range=99;
		my$random_number=int(rand($range)+1);
		
		if($GC_value<1){$GC_value=$GC_value*100;}
		
		if ($random_number <= $GC_value) {
		
			if(($random_number%2)==0){$letter="G";}
			else{$letter="C";} 
		}
	
		else {
			if(($random_number%2)==0){$letter="A";}
			else{$letter="T";}
			
		}

		$newseq.= $letter;

	}
	print OUT">Random_seq $j\n$newseq\n";
if(defined $opt_v){	print ">Random_seq $j\n$newseq\n";}

}