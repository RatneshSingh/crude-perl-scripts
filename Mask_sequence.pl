#!/usr/bin/perl -w
use strict;
use Getopt::Std;

our ($opt_s,$opt_i,$opt_l,$opt_o,$opt_r);
getopt('silor');

my$usage="Mask fraction of sequence by N
usage; perl script options
-s	Sequence file in fasta format
-i	start site to mask from (deafult=1)
-l	length of masked sequence
-r	remove the sequence and replace with 100 Ns [default False]
-o	output_file to save masked fasta (default:masked_sequencefile.fasta)
";


print "$usage\n";




$/="\n>";    
open FASTA,"$opt_s" or die "\nCannot open Sequence file\n";
$opt_o= "Masked_".$opt_s if !defined $opt_o;
open OUT,">$opt_o" or die "Cannot open output file\n$usage";


while(<FASTA>){
    	chomp;
    	(my$header,my@sequence)=split("\n",$_);
    	
    	$header=~s/>//;						# Remove Leading > from Header
    	$header=~s/\s*$//;					# Remove trailing spaces from header
    	$header=~s/^\s*//;					# Remove Leading spaces from Header
    	
    	my$sequence= join("",@sequence);
    	$sequence=~ s/\s//g;
    	$sequence=~s/\n//g;

my $N=();
if (defined $opt_r){$N="N"x100;}else{ $N="N"x$opt_l; }

substr($sequence,$opt_i,$opt_l)=$N;
print OUT">$header\n$sequence";
}