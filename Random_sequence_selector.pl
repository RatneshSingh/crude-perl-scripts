#!/usr/bin/perl -w
use strict;
use Getopt::Std;

our($opt_s,$opt_l,$opt_n,$opt_r);
getopt('slnr');

my$usage="Usage: perl script options.....
-s	sequence 
-l	length_of_fragment 
-n	Number of fragments to sample[100] 
-r	remove the Ns[False]\n";

print "\n$usage\n";
die "Cannot find Seqenuce file\n$usage" if !defined $opt_s;
$opt_l=500 if !defined $opt_l;


my $output=$opt_l."_nt_".$opt_n."_".$opt_s;
open OUT,">$output";
my%sequence=ReadFasta($opt_s);

foreach (keys%sequence){ # for every sequence
	my$seq_length=length($sequence{$_});
	my$count=0;
	my$seq_name=$_;
	
	
	for(my$i=0;$i<$opt_n;$i++){ # start choppping sequence from start
		my$random_number=int(rand($seq_length-$opt_l-1)); # generate random start site
		my$end_site=$random_number+$opt_l;
		my$start_site=$random_number+1;
		
		if(($seq_length-$random_number)>$opt_l){#for fragments larger than chopping size
			my$dna_fragment=substr($sequence{$_},$random_number,$opt_l);
			print OUT">$seq_name($start_site - $end_site)\n$dna_fragment\n";
		}
		else{#for fragments smaller than chopping size.

			my$dna_fragment=substr($sequence{$_},$random_number);
			print OUT">$seq_name($start_site - $seq_length)\n$dna_fragment\n";
			print "Sequence smaller than requested length\n"; 
			}
	
		$count++;
	}
	print "$count fragments of $opt_l bases each from $_ .\n";

}

#-------------------------------------Start ReadFasta---------------------------------------+

sub ReadFasta{ # to read fasta format files into hash. returns hash.
	
	my $seqfile=shift(@_);
	
	
	my ($header,@sequence);
	chomp $seqfile;
	open FASTA,"$seqfile";
	print "reading Sequences from input file.....Plz wait...\n";
	my%seq_hash=();
	$/="\n>";    # Change record seperator to read Fasta
	my$sequence1=();
	while(<FASTA>){
    	chomp;
    	($header,@sequence)=split("\n",$_);
    	
    	$header=~s/>//;						# Remove Leading > from Header
    	$header=~s/\s*$//;					# Remove trailing spaces from header
    	$header=~s/^\s*//;					# Remove Leading spaces from Header
    	if($header=~/^\s*$/){next;}

    	my$sequence= join("",@sequence);
    	$sequence=~ s/\s//g;
    	$sequence=~s/\n//g;
    	if(lc($opt_r) eq lc('true')){$sequence=~s/N//gi;} # remove Ns from the sequence
    	$sequence1.=$sequence; # join the sequences into one sequence
    	$header=();
    	@sequence=();
    	$sequence="";		
	}
	$seq_hash{$seqfile}=$sequence1;     #feed headers and sequences in hash.

	my @seq_count=keys (%seq_hash);
	my $seq_count=@seq_count;
	my $seqlength=length($sequence1);
	print "Done....\nNumber of sequences read form input file = $seq_count\n\n";
	print "Done....\nTotal length of sequence read form input file = $seqlength\n\n";
	@seq_count=();
	$/="\n";    # Record seperator set back to default newline.

	return(%seq_hash);

}
#-------------------------------------End ReadFasta---------------------------------------+