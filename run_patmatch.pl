#! /usr/bin/perl -w
use strict;
use Getopt::Std;

getopt('rpsmM');

our($opt_r,$opt_p,$opt_s,$opt_m,$opt_M);

open PATTERN,"$opt_p" or die "Cannot find pattern file\n";
#my $opt_r; # Type of residue n: nucleotide,p: protein, c: complementary strand
#my $opt_s; # name of sequence file in unjustified fasta format. use unjustify.pl to prepare your sequences
#my $opt_m; # Number of mismatch.
#my $opt_M; # type of mismatch. ids for insertions, deletions, and substitutions; i   for insertions only ;id  for insertions and deletions.

$opt_m=0 if !defined $opt_m;
$opt_M=() if !defined $opt_M;
$opt_r='n' if !defined $opt_r;

my$process_output='|grep -c ">" ';
my$process_output_uniq='|grep ">"|cut -d ":" -f 1|sort|uniq|grep -c ">"';
my $seqsinfile1="grep -c '>' $opt_s";
my $seqsinfile=`$seqsinfile1`;
$seqsinfile=~s/\n//;


# process fasta file for length info
my $length1=0;
$/="\n>";
open FASTA,$opt_s or die "Cannot find Sequence file\n";
while(<FASTA>){
my($header,$seq)=split(/\n/,$_);

my$length= length $seq;
$length1=$length1+$length;
$length=0;
$header=$seq=();
}

print "Length of all sequences;$length1\n";

$/="\n";
print "Element\tFreq\tSequences hit\tSequences searched\tUniq_freq_per_kb\tTotal_freq_per_kb\n";
while(<PATTERN>){

	s/\s*//g;
	
	
	#command to run patmatch: perl scan_pipeline.pl [residue] [pattern] [sequence file] [number of mismatches] [mismatch types]
	my @command=('perl scan_pipeline.pl',"-$opt_r",$_,$opt_s,$opt_m,$opt_M);
	pop @command;
	my $command1=join(" ",@command,$process_output);
	
	my $Total_hits=`$command1`;
	
	my $command2=join(" ",@command,$process_output_uniq);
	#print "sending comnd: $command2\n";
	my $uniq_seq_hits=`$command2`;
	$uniq_seq_hits=~s/\n//;
	$Total_hits=~ s/\n//;
	
	#print "sending command: $command2\n";
	#my$freq=`$command2`;
	#$freq=~s/\n//;
	
	
	#
	my$Uniq_freq_per_kb;
	my$Total_freq_per_kb=$Total_hits*1000/$length1;
#	my$Uniq_freq_per_seq=$Total_hits/$uniq_seq_hits;
	if($uniq_seq_hits){$Uniq_freq_per_kb=$Total_hits/$uniq_seq_hits;}
	else{$Uniq_freq_per_kb=0;}
	
	print "$_\t$Total_hits\t$uniq_seq_hits\t$seqsinfile\t$Uniq_freq_per_kb\t$Total_freq_per_kb\n";
	
	


}

exit;