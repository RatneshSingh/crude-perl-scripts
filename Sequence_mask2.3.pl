#!usr/bin/perl -w
################################################################################################
#This program will read blast (blastn) tabular output file and sequence(subject)file. 
# It will mask query or subject on the basis of given alignment cordinates.
#=>'program_name'  -s 'sequence_file'  -b 'blast_output_file' -o 'output_file' -m q(uery)|s(ubject))
################################################################################################


use strict;
use Getopt::Std;


our($opt_s,$opt_b,$opt_o,$opt_m);
getopt('sbom');



my($fastafile,$blast_table,$qlength,$output, %seq,@seqheaders, $header, $sequence,@sequence,$extra,$key,$revcom);
my$usage= "\n".'*******************************************************************
This script is made to mask query or subject sequences at the co-ordinates 
given by a blast tabular outputfile.

usage: perl script -options

options:
-s	Sequence file to mask in fasta format
-b	blast result in table format
-m	q|s. To mask q(uery) or s(ubject). Respective co-oridinates will be used based on selection.[q]
-o	Output file to save masked sequences.[Masked_sequence.fa]



***********************************************************************'."\n\n";

#*********************************************************************************************************
# open fasta file containing subject sequences used in the blast as database..
if($opt_s){
$fastafile=$opt_s;
}
else{print "$usage\n"; print "Type the name of file containing sequences in fasta format\n";
$fastafile=<STDIN>;
}
chomp($fastafile);

#open (FASTAFILE,$fastafile) or die"Cannot open the file $fastafile";

#********************************************************************************************************
# open blast output file in tabular format. here it is tblastn where one query is blasted against a set 
# of subject sequences. query is  aminoacid file and subject nucleotide.
#

if($opt_b){
$blast_table=$opt_b;
}
else{print "$usage\n"; print "Type the name of file containing tabular blast output\n";
$blast_table=<STDIN>;
}
chomp($blast_table);

open (BLASTFILE,$blast_table) or die"Cannot open the file $blast_table";

###############################################################
#out put file
if(!$opt_o){$opt_o="Masked_sequence.fa";} chomp($opt_o);
open (OUT,">$opt_o") or die"Cannot open the file $opt_o";

$extra=0;


if(!$opt_m){$opt_m='q';};
chomp($opt_m);

%seq=ReadFasta($fastafile);

#********************************************************************************************************
# parse blast tabular output to find out the cordinates to chop the sequence
#

$/="\n";   # define record seperator for blast file. current record seperator is "\n>" as defined for sequence file.
while(<BLASTFILE>){#1
	chomp;
	next if $_ =~ /^\s*$/;
	my @columns=split(/\s+/,$_);
	
	if(scalar(@columns)>12){print "More than 12 columns are found in blast file. \nPlease Make sure that there are no white spaces in the query or subject names.\n\n"};
	#chomp ($query,$target,$identity,$length,$mismatch,$gap,$qstart,$qend,$sstart,$send,$evalue,$bitscore);
	
	#print "$target,$identity,$length,$mismatch,$gap,$qstart,$qend,$sstart,$send,$evalue,$bitscore\n\n";
	my $extend=0;
	#calculate starting and end point for chopping subject sequence.
	our($target,$sstart,$send);


	if($opt_m eq 'q'){
		$target=$columns[0];
		$sstart=$columns[6];
		$send=$columns[7];

	}
	
	else{
		
		$target=$columns[1];
		$sstart=$columns[8];
		$send=$columns[9];

	}
	
	
	if ($sstart<$send){#2
		$sstart=$sstart-$extra-1;
		$send=$send+$extra;
		$extend=$send-$sstart;
	#$revcom=0;
						}#2
	
	if($sstart>$send) {#3
		my	$sstart_N=$send-$extra-1;
		my	$send_N=$sstart+$extra;
			$extend=$send_N-$sstart_N;
			$sstart=$sstart_N;
			$send=$send_N;
			#$revcom=1;
	}#3

	mask_seq($target,$sstart,$extend);
		
	
}#1


########################################################################
# write modified database

foreach $key(keys %seq){#2
#foreach $key(@seqheaders){

print OUT">$key\n$seq{$key}\n";
print "\nWriting Seq $key";		

}

print " \n\nOut put saved in $opt_o file.\n\n";
close FASTAFILE;
close BLASTFILE;
exit;
######################################################################################################
# subroutine to search header matching subject name and mask on the basis of given cordinates. return choped substring.
######################################################################################################
sub mask_seq {#1

my ($target,$sstart,$extend)=@_;

# my @subject=split(/\|/,$target);  ### use this line if sequences have Accession number or genebank naming format

#my @subject=split(/\s+||\|/,$target);  ### Use this line if Header has unique elements seperated by spaces.

my @subject=($target);  ### Use this line if Header has to be used as whole for identification.



my $accession=$subject[0]; ### change the value in bracket here to define which element should be used for identification.




chomp($accession);

foreach $key(keys %seq){#2
	
	my$keyNew=$key.'X';
	my$accessionNew=$accession.'X';
#	if ($keyNew=~/$accessionNew/i){#3
	if ($key eq $accession){

	my $new_sequence1 = substr($seq{$key},$sstart,$extend);
			$new_sequence1 =~ tr/ATGCatgc/XXXXXXXX/;
			substr($seq{$key},$sstart,$extend)= $new_sequence1;
			my $ssstart=$sstart+1;

	print "For Blast hit $accession Masked $extend N's from \n$key\n at $ssstart\n\n\n"; 

#return ($key,$new_seq);
}#3

}#2

}#1
######################################################################################################################

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
##############################################################################################

