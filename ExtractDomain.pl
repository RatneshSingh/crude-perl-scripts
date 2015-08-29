use strict;
use warnings;
open FASTA,"$ARGV[0]";

my %seq;
my @fulllist;

$/="\n>";
while (<FASTA>){

	chomp;
	#split record into header and sequence and feed to $header and @seq
	(my$header,my@seq)= split (/\n/,$_);
	$header=~s/>//;	
	#$header=~s/\\/_/g;
	$header=~s/\s*$//g;
								
	#join elements of @seq array to remove any spaces.
	my$seq=join ("",@seq);
	chomp ($header,$seq);

	push (@fulllist,$header);

	$seq{$header}=$seq;

	}


foreach my$seq2(keys %seq){
my$count=0;
$count++ while $seq{$seq2}=~/[VILMF][\s-]*[VILMF][\s-]*[VILMF][\s-]*[VILMF][\s-]*D/gi;

print "\nThere are $count domains in sequence $seq2";
my@domains=$seq{$seq2}=~/([VILMF][\s-]*[VILMF][\s-]*[VILMF][\s-]*[VILMF][\s-]*D)/gi;
 print "\t@domains";


}



