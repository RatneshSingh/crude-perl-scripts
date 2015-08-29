#! /usr/bin/perl -w
use strict;
use Getopt::Std;

our($opt_s);
getopt('s');


print "Sequence file : $opt_s\n";
my $out='cleaned_'."$opt_s";
open FASTAOUT,">$out";


ReadFasta($opt_s);

close FASTA;
close FASTAOUT;

print "Sequences were saved in $out\n ";







#-------------------------------------Start ReadFasta---------------------------------------+

sub ReadFasta{ # to read fasta format files into hash. returns hash.
	
	my $seqfile=shift(@_);
	
	
	my ($header,@sequence);
	chomp $seqfile;
	open FASTA,"$seqfile";
	print "reading Sequences from input file.....Plz wait...\n";
	my%seq_hash=();
	$/="\n>";    # Change record seperator to read Fasta
	my$count=1;
	while(<FASTA>){
    	chomp;
    	($header,@sequence)=split("\n",$_);
    	
    	$header=~s/>//;						# Remove Leading > from Header
    	$header=~s/\s*$//;					# Remove trailing spaces from header
    	$header=~s/^\s*//;					# Remove Leading spaces from Header
    	
    	my$sequence= join("",@sequence);
    	$sequence=~ s/\s//g;				#remove spaces from sequences
    	$sequence=~s/\n//g;					# remove new lines from sequences
    	$sequence=~s/[^ATGC]//gi;		# Remove evrything non ATGC
    	if($header=~/^\s*$/){next;}			# nex if header is empty line
    
 		
 #   	$seq_hash{$header}=$sequence;     #feed headers and sequences in hash.
 
		print FASTAOUT">$header\n$sequence\n";
		$count++;
	}
	
#	my @seq_count=keys (%seq_hash);
	my $seq_count=$count;

	print "Done....\nNumber of sequences read form input file = $seq_count\n\n";
	#@seq_count=();
	$/="\n";    # Record seperator set back to default newline.

#	return(%seq_hash);
return();

}
#-------------------------------------End ReadFasta---------------------------------------+