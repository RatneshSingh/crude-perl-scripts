#!/usr/bin/perl -w
use strict;
use Getopt::Std;

our($opt_s,$opt_l);
getopt('sl');

my$usage="Usage: perl script -s sequence -l length_of_fragment\n";
die "Cannot find Seqenuce file\n$usage" if !defined $opt_s;
$opt_l=500 if !defined $opt_l;


my $output=$opt_l."_nt"."_$opt_s";
open OUT,">$output";
my$sequence=ReadFasta($opt_s);

#foreach (keys%sequence){
	my$seq_length=length($sequence);
	my$count=0;
	my$seq_name='papaya_all_concat';
	for(my$i=0;$i<=$seq_length;$i=$i+$opt_l){
		
		my$end_site=$i+$opt_l;
		my$start_site=$i+1;
		
		if(($seq_length-$i)>$opt_l){#for fragments larger than chopping size
			my$dna_fragment=substr($sequence,$i,$opt_l);
			print OUT">$seq_name ($start_site - $end_site)\n$dna_fragment\n";
		}
		else{#for fragments smaller than chopping size.

			my$dna_fragment=substr($sequence,$i);
			print OUT">$seq_name ($start_site - $seq_length)\n$dna_fragment\n";
			#last;
		}
	
		$count++;
	}
	print "$count fragments of $opt_l bases each from $_ .\n";

#}

#-------------------------------------Start ReadFasta---------------------------------------+

sub ReadFasta{ # to read fasta format files into hash. returns hash.
	
	my $seqfile=shift(@_);
	
	
	my ($header,@sequence);
	chomp $seqfile;
	open FASTA,"$seqfile";
	print "reading Sequences from input file.....Plz wait...\n";
	
	$/="\n>";    # Change record seperator to read Fasta

	while(<FASTA>){
    	chomp;
    	($header,$sequence)=split("\n",$_);
    	
    	$header=~s/>//;						# Remove Leading > from Header
    	$header=~s/\s*$//;					# Remove trailing spaces from header
    	$header=~s/^\s*//;					# Remove Leading spaces from Header
    	
#    	my$sequence= join("",@sequence);
#    	@sequence=();
    	$sequence=~ s/\s//g;
    	$sequence=~s/\n//g;
#    	if($header=~/^\s*$/){next;}
#    	$seq_hash{$header}=$sequence;     #feed headers and sequences in hash.
		
	}
	
#	my @seq_count=keys (%seq_hash);
#	my $seq_count=@seq_count;

#	print "Done....\nNumber of sequences read form input file\n\n";
#	@seq_count=();
	$/="\n";    # Record seperator set back to default newline.

	return($sequence);
#	%seq_hash=();
$sequence=();
}
#-------------------------------------End ReadFasta---------------------------------------+