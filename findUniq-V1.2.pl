#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;

#####################################################################################################
# This script is made to mask sequence regions showing hits in blast file. It will read the coordi- #
# nates given in the blast file (tabular output) and mask the respective query sequence region      #
#  with XXXX.                                                                                        #
#                                                                                                   #
# Author : Ratnesh Singh                                                                            #
# version 1.0                                                                                       #
# contact for bugs: ratnesh@hawaii.edu                                                              #
#####################################################################################################


getopts('s:b:o:rln:x:');
our ($opt_s,$opt_b,$opt_o,$opt_l,$opt_r,$opt_n,$opt_x);
our (%seq);
my($new_header,$new_sequence);

# like the shell getopt, "d:" means d takes an argument
print "-sequence file: $opt_s\n" if defined $opt_s;
print "-blast file: $opt_b\n" if defined $opt_b;
print "-print output as: $opt_o\n" if defined $opt_o;
print "Unprocessed by Getopt::Std:\n" if $ARGV[0];
foreach (@ARGV) {
  print "-- $_\n";
}

my $help="This script will read sequence name and cordinates from
blast file and extract sequences from sequence file provided
 usage: perl $0 options
options:
-s sequence file
-b Blast out put file in table format. If not provided blast will be run for sequence file.
-l convert blast hit regions to lower case or soft masked.
-x replace blasthit regions to this character [default is to change to N]
-r Run repeat masker to mask repeats.
-n Species name for which repeats to be masked[Viridiplantae]
-o output file [Default is: output_seq_extract]

";


die "\nThere is no sequence file specified with -s \n $help" if !defined ($opt_s);
print "\nThere is no blast file specified with -b so script will run blastn for $opt_s against itself. \n " if !defined ($opt_b);



##########################################################################
####### run repeat masker to soft mask repeats.
if($opt_r){
  $opt_n=$opt_n?$opt_n:"Viridiplantae";
  my$mask_to=$opt_l?'-xsmall':'-x';
  system("RepeatMasker -e ncbi -pa 10  -s -species $opt_n $mask_to $opt_s");
  $opt_s.=".masked";
}

##########################################################################
######## run blast to find self hits
$opt_o = "$opt_s.selfBlastnXMasked.fasta" if (!defined $opt_o && !defined $opt_l);
$opt_o = "$opt_s.selfBlastnSoftMasked.fasta" if (!defined $opt_o);
print "saving output in:$opt_o\n";
open OUT,">$opt_o" or die "cannot open Output file\n";


$opt_x=$opt_x?$opt_x:"N";



##### run blast againast self. assumes makeblastdb and blastn are in PATH
if (!defined $opt_b){
  my $blasts_command="blastn -task blastn -query $opt_s -subject $opt_s -outfmt 6 -out $opt_s.selfblastn.table";
  system($blasts_command);
  $opt_b="$opt_s.selfblastn.table";
}
####


print "reading sequence file....Plz wait\n";
open FASTA,"$opt_s" or die "cannot open Sequence file\n";

$/="\n>";
while(<FASTA>){
	chomp;

	my($header,@sequence)=split(/\n/,$_);
	  $header=~s/>//;
	  $header=break_at_space($header);

	my$sequence= join("",@sequence);
	$sequence=~s/\s+|>//g;
	$seq{$header}=$sequence;

}
close (FASTA);
print".............Done\n";
$/="\n";



###################################################################
#parse information from blast file and extract seq start and end  #
###################################################################
open BLAST,"$opt_b" or die "cannot read blast file \n";

while(<BLAST>){
	my $line=$_;
#	print "before line : $line\n";
	if ($line=~/^\s*$/){ next ;};
	if ($line=~/^\#/){ next; } ;
#	print "line: $line\n";
	my @line_info= split(/\t/,$line);

	my$isselfhit=0;
	$isselfhit=self_hits($line);
	next if $isselfhit>0;


	my $query=break_at_space($line_info[0]);

#	print "query:$query\n";
	my $subject=$line_info[1];
#	print"subject:$subject\n";
	my $qstart=$line_info[6];
#	print "Extracting --> $subject: sstart: $sstart\t";
	my $qend=$line_info[7];
	my $qstartN=$qstart-1;
	my $len = $qend - $qstartN;
	print "Masking --> $query: qstart: $qstartN\t qend: $qend \t length : $len \n";

	 ($new_header,$new_sequence)=mask_seq($query,$qstartN,$len);
	#print OUT">$new_header\n$new_sequence\n" if defined $new_sequence;
}

foreach (keys %seq){
   print OUT">$_\n$seq{$_}\n";
}

close(OUT);
close(BLAST);
exit;


################################################################################
# subroutine for extraction of sequences
################################################################################

sub mask_seq{
	my($query,$qstartN,$len1)=@_;
	my $new_sequence1 = substr($seq{$query},$qstartN,$len1);

			 #$new_sequence1 =~ tr/ATGCatgc/NNNNnnnn/

			  my$uc_x=uc$opt_x;
			  $new_sequence1 =~ s/[ATGC]/$uc_x/g;

			  my$lc_x=lc$opt_x;
			  $new_sequence1 =~ s/[atgc]/$lc_x/g;

			  if($opt_l){
			  $new_sequence1 =lc$new_sequence1;
			  }

			substr($seq{$query},$qstartN,$len1)= $new_sequence1;
			return ($query,$seq{$query}) if defined $new_sequence1;
	}
#
#
#sub clean_nonspecific{
#  my$string=shift;
#  $string=~s/[^\w\d\s]/_/g;
#  return($string);
#}

sub break_at_space{
  my$string=shift;
  $string=~s/\s+.*//;

  return($string);
}

sub self_hits{
  my$bline=shift;
  $bline=~s/^\s+//g;
  my@belement=split(/\t/,$bline);

  #### sometimes blast adds lcl| at the start of database entries which can cause problems in finding selfhits
  $belement[0]=~s/^lcl\|//g;
  $belement[1]=~s/^lcl\|//g;

  if(($belement[0] eq $belement[1]) && ($belement[6]+$belement[7])==($belement[8]+$belement[9]))
  {
	#print "Selfhit found at $bline\n";
	return 1
  }
  else {return 0}
}
