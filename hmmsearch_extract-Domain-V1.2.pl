#!/usr/bin/perl -w
#use warnings;
#use strict;
use Getopt::Std;

#####################################################################################################
# This script is made to extract domains  from sequence file provided by the user, based on the     #
# cordinates given in the HMM output file regular output. It will create tabular output and use it  #
# for extraction co-ordinates.                                                                      #
#                                                                                                   #
# Author : Ratnesh Singh                                                                            #
# version 1.2                                                                                       #
# In case of Bugs contact:ratnesh@hawaii.edu                                                        #
#####################################################################################################


getopt('sho');
our ($opt_s,$opt_h,$opt_o);
our (%seq);


# like the shell getopt, "d:" means d takes an argument
print "-sequence file: $opt_s\n" if defined $opt_s;
print "-hmm coordinates file/sequence: $opt_h\n" if defined $opt_h;
print "-print output as: $opt_o\n" if defined $opt_o;
print "Unprocessed by Getopt::Std:\n" if $ARGV[0];
foreach (@ARGV) {
  print "-- $_\n";
}

my $help="\n\nThis script will read sequence name and cordinates from \n HMM output file and extract domains from sequence file provided\n usage:\n use following options\n -s sequence file \n -h HMM output file in table format \n -o output file [Default is: output_seq_extract]\n\n\n";
die "\nThere is no sequence file specified with -s \n $help" if !defined ($opt_s);



open FASTA,"$opt_s" or die "cannot open Sequence file\n";

print"Reading Fasta file";


$/="\n>";
while(<FASTA>){
	chomp;
	
	my($header,@sequence)=split(/\n/,$_);
	  $header=~s/>//;
#	  $header=~s/\|/_/g;
	  my @header= split(/\s+/,$header);
	  $header=$header[0];

	my$sequence= join("",@sequence);
	$sequence=~s/\s//g;
	$sequence=~s/\n//g;
	
	$seq{$header}=$sequence;
	#print "$header\n";
	#print"$sequence\n\n";
		
}
close (FASTA);
print".............Done\n";


open OUT,">$opt_o" or die "cannot open Sequence file\n";


#---------------------------------------------------------------
# Creating blast file from Hmm output
#---------------------------------------------------------------

#foreach my $file (@files) {

#print "file being read: $opt_h\n";

open BLAST,">$opt_h.Table" or die "Cannot open $opt_h.Table\n";
	{
    local( $/, *HMM ) ; # Slurping whole file
    open( HMM, $opt_h ) or die "sudden flaming death\n";
    my $text = <HMM>;
	#print $text;
	my $Pfam=();
	my @fragment= split(/Alignments(\s)*of(\s)*top-scoring(\s)*domains:/,$text);
	
	if($fragment[0]=~/HMM file:[\s]+(PF[\d]+)/){$Pfam=$1};
	
	$fragment[0]=~ s/hmmsearch [\s\S\w\W\d\D]*?Parsed[\s\S\w\W\d\D]*?E-value[-\s]+//;
	print BLAST"$fragment[0]";
			
	}
	





#---------------------------------------------------------------
# Reading Blast file
#---------------------------------------------------------------

$/="\n";
open BLAST,"<$opt_h.Table" or die "cannot read blast file \n";
	
		while(<BLAST>){
			my $line=$_;
			#	print "before line : $line\n";
			if ($line=~/^\s+$/){ next ;};
			#	print "line: $line\n";
			my @line_info= split(/\s+/,$line);
			my $query=$line_info[0];
#			$query=~s/\|/-/g;
	  		my @query= split(/\s+/,$query);
	  		$query=$query[0];

						
			#	print "query:$query\n";
			#my $subject=$line_info[1];
			#	print"subject:$subject\n";
			my $sstart=$line_info[2];
			#	print "Extracting --> $subject: sstart: $sstart\t";
			my $subend=$line_info[3];
			#	print " subend:$subend\n";
	
			#adjustment for increase in length at both ends
	
			my($start_s,$end_s);
	
			if($sstart > $subend){ $end_s=$sstart - 1; $start_s=$subend - 1;}
			else{$start_s = $sstart - 1; $end_s = $subend - 1 ; }
			my $len = $end_s - $start_s + 1;
			#print "Extracting --> $query: sstart:".($start_s+1)."\t end:". ($end_s+1). "\t length : $len \n";

	
			my ($new_header,$new_sequence)=extract_seq($query,$start_s,$end_s,$len);
			
#			$new_header=~s/_/\|/g;
			
			print "\tStart:$start_s\tEnd:$end_s\tDomain length:$len\n";
			print OUT">$new_header\n$new_sequence\n" if defined $new_sequence;
		}


##################################################################################
sub extract_seq{
	my($query1,$start_s1,$end_s1,$len1)=@_;
	
	my $new_sequence1= substr($seq{$query1},$start_s1,$len1) if defined $seq{$query1};
	my $length_seq = length($seq{$query1}) if defined $seq{$query1};
	print "$query1 length : $length_seq" if defined $length_seq;
	#my $new_header1=$query1.'-'.'start-'.($start_s1+1).'_'.'end-'.($end_s1+1); # use this to add start, stop and length with header
	my $new_header1=$query1;  													# Use this to use simple header.				

	return ($new_header1,$new_sequence1) if defined $new_sequence1;
	print "Cannot find $query1\n" if !defined $seq{$query1}
	}

#################################################################################