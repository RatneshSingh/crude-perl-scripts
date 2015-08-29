#!/usr/bin/perl -w
use strict;
use Getopt::Std;
# this script can remove sequences from multi fasta files based on % of masked region and length. 
our($inputfile,$outputfile,$count,%seq_hash);
our($opt_s,$opt_o,$opt_l,$opt_p,$opt_n);

getopt('solpn');


my$usage= "
This script will filter sequences masked by the any masker programs from non masked ones. User can provide filtering criteria.  

\n\nusage : perl script_name [options]

options:
-s	Sequence file containing sequences in fasta format
-l	Min Length of the sequence to keep [1].
-n	Min Number of non-masked nucleotide in sequence [0]. 
-p	Min Percent of sequence mask allowed [100].
-o	Output file to store resulting sequences.

";

if(!$opt_s){print $usage;}

#opening inputfile containing sequences in fasta format.
if ($opt_s){$inputfile=$opt_s}elsif($ARGV[0]){$inputfile=$ARGV[0]}else{die "\nPlease provide sequence file to filter\n"}
chomp ($inputfile);
open INPUT,"<$inputfile" or die "Cannot open $inputfile.....\n\n";


$opt_l = 1 if !defined $opt_l;
$opt_p = 100 if !defined $opt_p;
$opt_n = 0 if !defined $opt_n;


# Read fasta file in hash. Change the names if name s not uniq
%seq_hash=ReadFasta($inputfile);


#opening outputfile
if ($opt_o){$outputfile=$opt_o;} else{$outputfile='maskfiltered_'.'l-'.$opt_l.'-P_'.$opt_p.'_'.'n-'.$opt_n.'_'.$opt_s;}
if($outputfile){
chomp($outputfile) ;
open OUT,">$outputfile" or die "cannot create $outputfile.....\n\n" ;
}

#*************************************************************************************************************************
#read all the pattern or key words into array @pattern. and search for the matching one in abovemade hash.
my $count_found=my $count_notfound=0;
print"Filtering Masked sequences...Plz wait....\n";
foreach my$Seq_header(keys %seq_hash){
	
	# Count Ns in sequence and calculate % of sequence masked
	my$mask=($seq_hash{$Seq_header}=~tr/Nn/Nn/);
	my$non_mask = ($seq_hash{$Seq_header}=~tr/ATGCatgc/ATGCatgc/);
	my$seq_length = length($seq_hash{$Seq_header});
	my$mask_percent=(100*$mask)/$seq_length;
	
	#filter for % mask and length
	if($seq_length >= $opt_l && $mask_percent <= $opt_p && $non_mask >= $opt_n){
		print OUT">$Seq_header\n$seq_hash{$Seq_header}\n";
		$count_found++;
	}
	else{$count_notfound++;}
	
}

print "\nNumber of sequences written to output = $count_found\n";
print "Number of sequences excluded = $count_notfound\n";


#*************************************************************************************************************************
#exit the program and close all files.
close INPUT;
close OUT;
exit;






#######################################################
#-------------------------------------Start ReadFasta---------------------------------------+

sub ReadFasta{ # to read fasta format files into hash. returns hash.
	
	my $seqfile=shift(@_);
	
	
	my ($header,@sequence);
	chomp $seqfile;
	open FASTA,"$seqfile";
	print "reading Sequences from input file.....Plz wait...\n";
	my%seq_hash=();
	#$seq_hash{'RS_Concatenated'}="";
	
	$/="\n>";    # Change record seperator to read Fasta
	my$last_N=1;
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
    	
    	# Make sure no duplicate header exists else they will be overwritten and lost during hashing due to same key.
    	if(!exists $seq_hash{$header}){
    		$seq_hash{$header}=$sequence;     #feed headers and sequences in hash.
			#$seq_hash{'RS_Concatenated'}=$seq_hash{'RS_Concatenated'}.$Ns.$sequence;
		}	
		else {
		
			# find a uniq header name by adding a number at the end. If header still exists, increase the number by one
			while(exists $seq_hash{$header}){$header=$header.$last_N;$last_N++;}
			
			$seq_hash{$header}=$sequence;
			#$seq_hash{'RS_Concatenated'}=$seq_hash{'RS_Concatenated'}.$Ns.$sequence;

		
		}
	}
	
	my @seq_count=keys (%seq_hash);
	my $seq_count=@seq_count;

	print "Done....\nNumber of sequences read form input file = $seq_count\n\n";
	@seq_count=();
	$/="\n";    # Record seperator set back to default newline.

	return(%seq_hash);

}
#-------------------------------------End ReadFasta---------------------------------------+
