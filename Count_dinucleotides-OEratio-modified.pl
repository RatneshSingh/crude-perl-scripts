#!/usr/bin/perl -w
use strict;
use Getopt::Std;

our($opt_f);
getopt('f');
my$usage="Usage:\nperl\tscript\tfile1\tfile2\tfile.........\n\n\nor\n\n perl script -f[file containing multiple sequences]\n\n";

print $usage;


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

open OUT,">Base_pair_freq_Summary.txt";
#transfer above created pairs in an array
foreach(keys%basepair){push(@basepairs,$_);}

#print headers for output file
print "\tG+C";
print OUT"\tG+C";

for(my$i=0;$i<=scalar@basepairs-1;$i++){print "\t$basepairs[$i]";print OUT"\t$basepairs[$i]";}
my$count;
if(!defined $opt_f){
# count frequencies forA,T,Gand C.= No.of N/length of seq under consideration

	foreach my$file(@ARGV){

		my$sequence=ReadFasta($file);
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
		printf OUT"\n%s\t%.3f",$file,$GCpcent;
	
		#Calculate frequency	
		for(my$i=0;$i<=scalar@basepairs-1;$i++){
	#		print "\tpair:$basepairs[$i]";
			$count++ while $sequence=~/$basepairs[$i]/gi;
			my$basefrequency=$count/$seq_length;
		
			printf"\t%.3f",$basefrequency;
			printf OUT"\t%.3f",$basefrequency;

			
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
}	
# Process each sequence in a file individually

elsif($opt_f){

	my%sequences= ReadFastaTohash($opt_f);
	
	foreach my$file(keys %sequences){
		

		my$sequence=$sequences{$file};
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
		printf OUT"\n%s\t%.3f",$file,$GCpcent;
	
	#Calculate frequency	
		for(my$i=0;$i<=scalar@basepairs-1;$i++){
#		print "\tpair:$basepairs[$i]";
			$count++ while $sequence=~/$basepairs[$i]/gi;
			my$basefrequency=$count/$seq_length;
		
			printf"\t%.3f",$basefrequency;
			printf OUT"\t%.3f",$basefrequency;

		
			$count=0;
		}

#Calculate Obs/Exp ratios

		printf "\n%s\tObs/Exp ratio",$file;
		printf OUT"\n%s\tObs/Exp ratio",$file;
	
		for(my$i=0;$i<=scalar@basepairs-1;$i++){
#		print "\tpair:$basepairs[$i]";
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

}	
print"\n\n";
printf OUT"\n\nFrequency of Bases:\nA\t%.3f\nT\t%.3f\nG\t%.3f\nC\t%.3f",$freq{'A'},$freq{'T'},$freq{'G'},$freq{'C'};	




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
    	@sequence=();
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
	%seq_hash=();
}
