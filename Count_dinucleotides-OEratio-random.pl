#!/usr/bin/perl -w
use strict;
use Getopt::Std;

our($opt_f,$opt_l,$opt_n,$opt_p,$opt_r,$opt_a);
getopt('fnlpra');

# Write usage rules and options
my$usage="
****************************************************************
This script calculates the Observed/Expected 
ratios of all the dinucleotide pair in a given sequence

Optionally this can be used to calculate the OE ratio
of randomally selected fragments of specified length
from the given sequence for statistical test purposes

Usage:

=>To calculate OE ratio for individual files 
(multiple sequences in each file will be merged)

	perl\tscript\tfile1\tfile2\tfile.........


=>To calculate OE ratio for each sequence in one file
(Each sequence in the file will be treated individually)\n 
	
	perl script -f[file containing multiple sequences] 
				-a Remove ambigous bases from sequence [false]


=> to calculate the OE ratio for random selected fragments

	perl script [options...........]
	
-f	file containing multiple sequences
-r	Only test the OE ratio in randomaly selected sequences
	(Exclusion of this option will calculate in random sequence)
-l	Length of randomly selected sequence [10000]
-n	Number of random sequences [100]
-a Remove ambigous bases from sequence [false]
-p	print frequencies in file [FALSE]
*****************************************************************
\n\n"
;

print $usage;


if(!defined @ARGV){print "\nParameters missing!!!\n"; die;}

#my%count;
my@basepairs;
my%basepair;
my%freq;
my@bases=("A","T","G","C");

#create all possible pairs of ATGC and save in hash
foreach my$base1(@bases){
	foreach my$base2(@bases){
		my$basep=$base1.$base2;
		#print "\nbase pair:$basep";
		$basepair{$basep}=0;
				
	}
}

if(defined $opt_f){my$ouputfile="OEratio_".$opt_f.".txt";
open OUT,">$ouputfile";}
else{my$outputfile="OEratio_"$ARGV[0].".txt"; open OUT,">$outputfile";}
#transfer above created pairs in an array
foreach(keys%basepair){push(@basepairs,$_);}

#print headers for output file
print "\tG+C";
print OUT"\tG+C";

for(my$i=0;$i<=scalar@basepairs-1;$i++){print "\t$basepairs[$i]";print OUT"\t$basepairs[$i]";}
my$count;
if(!defined $opt_r){
	# for individuals files with one sequence in each file
	if(!defined $opt_f){
		
		foreach my$file(@ARGV){
		if($file=~/^-/){next}
		my$sequence=ReadFasta($file);
		count_OE($sequence,$file);
		}	
	}
# When one file with multiple sequence is provided

	elsif($opt_f){
		my%sequences= ReadFastaTohash($opt_f);
		foreach my$file(keys %sequences){
			my$sequence=$sequences{$file};
			count_OE($sequence,$file);
		}	
	}

	close OUT;

}

if(defined $opt_r){

	$opt_l=10000 if !defined$opt_l;
	$opt_n=100 if !defined $opt_n;
	die "Cannot find sequence file" if !defined $opt_f;
	
		open OUT,">$opt_f.$opt_n.RandomFrag.$opt_l.Length.OEratio.out";
		print OUT"\t";
		for(my$i=0;$i<=scalar@basepairs-1;$i++){print "\t$basepairs[$i]";print OUT"\t$basepairs[$i]";}

		my$sequence1= ReadFasta($opt_f); # join multiple sequence to one
			
		for(my$i=1;$i<=$opt_n;$i++){
			my$length=length($sequence1);
			my$random=int(rand($length-$opt_l-1));
			my$sequence=substr($sequence1,$random,$opt_l);
			count_OE($sequence,$i);
		}

	
}

print"\n\n";
#printf OUT"\n\nFrequency of Bases:\nA\t%.3f\nT\t%.3f\nG\t%.3f\nC\t%.3f",$freq{'A'},$freq{'T'},$freq{'G'},$freq{'C'};	




#-------------------------------------Start ReadFasta---------------------------------------+

sub ReadFasta{ # to read fasta format files into hash. returns hash.
	
	my $seqfile=shift(@_);
	
	
	my ($header,@sequence);
	my$sequence2=();
	my$sequence=();
	chomp $seqfile;
	open FASTA,"$seqfile";
#	print "reading Sequences from input file.....Plz wait...\n";
	my%seq_hash=();
	$/="\n>";    # Change record seperator to read Fasta

	while(<FASTA>){
    	chomp;
    	($header,@sequence)=split("\n",$_);
    	
    	$header=~s/>//;						# Remove Leading > from Header
    	$header=~s/\s*$//;					# Remove trailing spaces from header
    	$header=~s/^\s*//;					# Remove Leading spaces from Header
    	
    	$sequence= join("",@sequence);
    	$sequence=~ s/\s//g;
    	$sequence=~s/\n//g;
    	if($header=~/^\s*$/){next;}
 #   	$seq_hash{$header}=$sequence;     #feed headers and sequences in hash.
		$sequence2.=$sequence;
	}
	
#	my @seq_count=keys (%seq_hash);
#	my $seq_count=@seq_count;

#	print "Done....\nNumber of sequences read form input file = $seq_count\n\n";
#	@seq_count=();
	$/="\n";    # Record seperator set back to default newline.
	$sequence=~s/\s+//g;
	return($sequence2);

}
#-------------------------------------End ReadFasta---------------------------------------+

sub ReadFastaTohash{
	
	 my$seqfile=shift(@_);
	 my($header,@sequence);
	chomp $seqfile;
	open FASTA,"$seqfile";
	print "reading Sequences from input file.....Plz wait...\n";
	my%seq_hash=();
	$/="\n>";    # Change record seperator to read Fasta

	while(<FASTA>){
    	chomp;
    	($header,@sequence)=split("\n",$_);
    	
    	$header=~s/>//;						# Remove Leading > from Header
    	$header=~s/\s*$//;					# Remove trailing spaces from header
    	$header=~s/^\s*//;					# Remove Leading spaces from Header
    	
    	my$sequence= join("",@sequence);
    	$sequence=~ s/\s//g;
    	$sequence=~s/\n//g;
    	
    	if($header=~/^\s*$/){next;}
    	$seq_hash{$header}=$sequence;     #feed headers and sequences in hash.
		
	}
	
	my @seq_count=keys (%seq_hash);
	my $seq_count=@seq_count;

	print "Done....\nNumber of sequences read form input file = $seq_count\n\n";
	@seq_count=();
	$/="\n";    # Record seperator set back to default newline.

	return(%seq_hash);
}

sub count_OE{

my $sequence=shift(@_);
my $file=shift(@_);	
my$seq_length=$sequence=~tr/ATGCatgc//;
		my$GCnumber=$sequence=~tr/GCgc//;
		my$GCpcent=$GCnumber/$seq_length;
	
		my$Anumber=$sequence=~tr/Aa//;
		$freq{'A'}=$Anumber/$seq_length;
	
		my$Tnumber=$sequence=~tr/Tt//;
		$freq{'T'}=$Tnumber/$seq_length;
	
		my$Gnumber=$sequence=~tr/Gg//;
		$freq{'G'}=$Gnumber/$seq_length;
	
		my$Cnumber=$sequence=~tr/Cc//;
		$freq{'C'}=$Cnumber/$seq_length;
	
		printf "\n%s\t%.3f",$file,$GCpcent;
		if($opt_p){printf OUT"\n%s\t%.3f",$file,$GCpcent;}
	
		#Calculate frequency	
		for(my$i=0;$i<=scalar@basepairs-1;$i++){
	#		print "\tpair:$basepairs[$i]";
			$count++ while $sequence=~/$basepairs[$i]/gi;
			my$basefrequency=$count/$seq_length;
		
			printf"\t%.3f",$basefrequency;
			if($opt_p){printf OUT"\t%.3f",$basefrequency;}

			
			$count=0;
		}

#Calculate Obs/Exp ratios

		printf "\n%s\tObs/Exp ratio",$file;
		printf OUT"\n%s\tObs/Exp ratio",$file;
	
		for(my$i=0;$i<=scalar@basepairs-1;$i++){
#			print "\tpair:$basepairs[$i]";
			$count++ while $sequence=~/$basepairs[$i]/gi;
			(my$b1,my$b2)=split("",$basepairs[$i]);
			$freq{$b1}=~s/[^\d.]//gi;
			$freq{$b2}=~s/[^\d.]//gi;
				
			my$baseOEratio=($count/$seq_length)/($freq{$b1}*$freq{$b2});
		
			printf"\t%.3f",$baseOEratio;
			printf OUT"\t%.3f",$baseOEratio;

		
			$count=0;
		}
	
	
}	

