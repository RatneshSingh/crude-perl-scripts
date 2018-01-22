#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;

#####################################################################################################
# This script is made to mask sequence regions showing hits in blast file. It will read the coordi- #
# nates given in the blast file (tabular output) and mask the respective query sequence region      #
#  with XXXX.                                                                                       #
#                                                                                                   #
# Author : Ratnesh Singh                                                                            #
# version 1.4		                                                                            #
# 3/10/2014 added min max subroutines to identify start and end as small and large numbers          #
# contact for bugs: ratnesh@hawaii.edu                                                              #
# 06/14/2017 added check for the existsence of multiple sequence with same header in fasta file     #
#####################################################################################################


getopts('s:b:o:rln:x:e:c:');
our ($opt_s,$opt_b,$opt_o,$opt_l,$opt_r,$opt_n,$opt_x,$opt_e,$opt_c);
our (%seq);
$opt_c||=20;
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
-e only mask regions with e-value lower than this [10]
-l convert blast hit regions to lower case or soft masked.
-x replace blasthit regions to this character [default is to change to N]
-r Run repeat masker to mask repeats.
-c Num cpu to use for RepeatMasker.
-n Species name for which repeats to be masked[Viridiplantae]
-o output file [Default is: output_seq_extract]

";


die "\nThere is no sequence file specified with -s \n $help" if !defined ($opt_s) ;
die "\nThe sequence file $opt_s does not exists. \n $help" if (! -e $opt_s) ;
print "\nThere is no blast file specified with -b so script will run blastn for $opt_s against itself. \n " if !defined ($opt_b);



##########################################################################
####### run repeat masker to soft mask repeats.
if($opt_r){
  print "\n Repeatmasker run is requested. Running Repeatmasker with following command.\n";
  $opt_n=$opt_n?$opt_n:"Viridiplantae";
  my$mask_to=$opt_l?'-xsmall':'-x';
  print "\nRepeatMasker -e ncbi -pa $opt_c  -s -species $opt_n $mask_to $opt_s\n";
  system("RepeatMasker -e ncbi -pa $opt_c  -s -species $opt_n $mask_to $opt_s");
  if(-e "$opt_s.masked"){$opt_s.=".masked";}
}

##########################################################################
######## run blast to find self hits
$opt_o = "$opt_s.selfBlastnXMasked.fasta" if (!defined $opt_o && !defined $opt_l);
$opt_o = "$opt_s.selfBlastnSoftMasked.fasta" if (!defined $opt_o);
print "saving output in:$opt_o\n";
open OUT,">$opt_o" or die "cannot open Output file\n";
print "*** -x and -l options are incompatible together. -l will take precedence. Masking to lower case.***\n" if ($opt_l && $opt_x);

$opt_x=$opt_x?$opt_x:"N";
#$opt_e=$opt_e?$opt_e:10;


##### run blast againast self. assumes makeblastdb and blastn are in PATH
if (!defined $opt_b){

  print "\n\n\nRunning blastn against self\n\n";
  my $blasts_command="blastn -task blastn -query $opt_s -subject $opt_s -outfmt 6 -out $opt_s.selfblastn.table";
  print "\n\n$blasts_command\n";

  system($blasts_command);
  $opt_b="$opt_s.selfblastn.table";
}
####


print "reading sequence file....Plz wait\n";
open FASTA,"$opt_s" or die "cannot open Sequence file\n";
my%count;
$/="\n>";
while(<FASTA>){
	chomp;

	my($head,@sequence)=split(/\n/,$_);
	  my$header=quotemeta($head);
    $header=$head;  # not using quotemeta
	  $header=~s/>//;
	  $header=break_at_space($header);

	my$sequence= join("",@sequence);
	$sequence=~s/\s+|>//g;
	$seq{$header}=$sequence;
	$count{$header}+=1;
	print "\nERROR: Found more than one sequence with same header $header. Masking will not be reliable \n" if $count{$header} > 1;
	exit if $count{$header} > 1;
}





close (FASTA);
print".............Done\n";
$/="\n";



###################################################################
#parse information from blast file and extract seq start and end  #
###################################################################
open BLAST,"$opt_b" or die "cannot read blast file \n";
print "\n\nparsing blast file to find hits\n ";
while(<BLAST>){
	my $line=$_;
#	print "before line : $line\n";
	if ($line=~/^\s*$/){ next ;};
	if ($line=~/^\#/){ next; } ;
#	print "line: $line\n";
	my @line_info= split(/\s+/,$line);
	my$qer=quotemeta($line_info[0]); $qer=$line_info[0];  ## not using quotemeta
	my$sub=quotemeta($line_info[1]); $sub=$line_info[1]; ## not using quotemeta


	$line_info[0]=$qer;
	$line_info[1]=$sub;
	my$isselfhit=0;
	$isselfhit=self_hits($line);
	next if $isselfhit>0;

	next if($opt_e && $line_info[10]>$opt_e);

	my $query=break_at_space($line_info[0]);

#	print "query:$query\n";
	my $subject=$line_info[1];
#	print"subject:$subject\n";
	my $qstart=min($line_info[6],$line_info[7]);
#	print "Extracting --> $subject: sstart: $sstart\t";
	my $qend=max($line_info[7],$line_info[6]);
	my $qstartN=$qstart-1;
	my $len = $qend - $qstartN;
	print "Masking --> $query: qstart: $qstart\t qend: $qend \t length : $len \n";

	 #($new_header,$new_sequence)=
	 mask_seq($query,$qstartN,$len);
	#print OUT">$new_header\n$new_sequence\n" if defined $new_sequence;
}

foreach (keys %seq){
   print OUT">$_\n$seq{$_}\n";
}

close(OUT);
close(BLAST);


print "\n\n";
exit;


################################################################################
# subroutine for extraction of sequences
################################################################################

sub mask_seq{
	my($query,$qstartN,$len1)=@_;
	#my$max_len=length($seq{$query})-$qstartN - 1;
	#$len1=$len1<$max_len?$len1:$max_len;
	if (! exists $seq{$query} ) {
        print "\nUnable to find $query in sequence hash table.";
		return;
    }

	print "\n**Masking $query with total length".length($seq{$query})." ans start:$qstartN  extLen:$len1";
	my $new_sequence1 = substr($seq{$query},$qstartN,$len1);

			  if($opt_l){
				#$new_sequence1 =~ tr/ATGCatgc/atgcatgc/
				$new_sequence1 =lc$new_sequence1;
			  }
			  else{
				  my$uc_x=uc$opt_x;
				  $new_sequence1 =~ s/[ATGC]/$uc_x/g;

				  my$lc_x=lc$opt_x;
				  $new_sequence1 =~ s/[atgc]/$lc_x/g;
			  }

			substr($seq{$query},$qstartN,$len1)= $new_sequence1;
			#return ($query,$seq{$query}) if defined $new_sequence1;
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
  $string=~s/\s+.*$//;

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


sub min{
  @_= sort { $a <=> $b }@_;
   return $_[0];

}

sub max{

 @_ = sort { $a <=> $b }@_;
   return $_[-1];


}
