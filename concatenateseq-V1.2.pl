#!/usr/bin/perl -w
use strict;
use Getopt::Std;
#this script will concatenate trailing fragments of nucleotide under one header into one line.
#In other words it will remove "\n" character from the fragments of sequence under one header.

our($opt_s,$opt_n,$opt_o,$opt_c);
getopt('snoc');





my($inputfile,$header,@sequence,$sequence,$outputfile,$add_Ns);

my$usage='
usage: perl scriptname -s sequence_file_in_fasta_format   options

options:
-n Number_of_Ns_to_add_While_Joining_multiple_sequences [10]
-o Output_file_to_save_concatenated_sequences [concatenated_inputfileName]
-c all|ind. Concatenate [all] sequences one header or concatenate fragments of [ind]ividual sequence under their own header [all].
';

print " \nWelcome!!! This script will Concatenate several fragment of one sequence or several sequences in a file under one Header\n$usage

";

if ($opt_s){$inputfile=$opt_s;}else{print "\n\n$usage\n\n\nEnter input file containing sequences\n"; $inputfile=<STDIN>;} chomp ($inputfile);

open INPUT,"<$inputfile" or die "Cannot open $inputfile.....\n\n";

if ($opt_o){$outputfile=$opt_o;}else{$outputfile='Concatenated_'.$inputfile;} chomp($outputfile);

if($opt_c){$opt_c=lc($opt_c);} else {$opt_c='all';}



if($opt_n){$add_Ns=$opt_n;}else{$add_Ns=10;} chomp($add_Ns);

$add_Ns=~s/\D//;
if($add_Ns<1 ||$add_Ns eq ""){$add_Ns=10;}

my$Ns='N'x $add_Ns;

open OUT,">$outputfile" or die "cannot create $outputfile.....\n\n";


my%seq_hash=ReadFasta($inputfile);

if($opt_c eq 'all'){
	
	open OUT,">>$outputfile" or die "cannot create $outputfile.....\n\n";
	print OUT ">$inputfile\n";
	#my$catseq="";
	foreach my$header(keys %seq_hash){
		
		print OUT $seq_hash{$header}.$Ns;

	}
	#print OUT ">$inputfile\n$catseq";

}

elsif($opt_c eq 'ind'){
	
	open OUT,">$outputfile" or die "cannot create $outputfile.....\n\n";
	foreach my$header(keys %seq_hash){
		print OUT ">$header\n$seq_hash{$header}\n";
	}
}







 
###################################################################################################################
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
######################################################################################################################
