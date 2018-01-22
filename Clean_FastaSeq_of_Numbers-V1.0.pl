#!/usr/bin/perl -w
use strict;


my($inputfile,$header,@sequence,$sequence,$outputfile,$patternfile,$pattern,$count,@pattern,$ta,%seq_hash);

print " \n****Welcome!!! This script will read pattern to search from the given file and pick out sequence containing pattern from another input file. can be used AS MANUAL MODE***  \n\n\nusage : perl script_name sequence_file outputfile\n\n";

#opening inputfile containing sequences in fasta format.
if ($ARGV[0]){$inputfile=$ARGV[0];}else{print "Enter input file containing sequences\n"; $inputfile=<STDIN>;}
chomp ($inputfile);
#open INPUT,"<$inputfile" or die "Cannot open $inputfile.....\n\n";


#opening outputfile
if ($ARGV[1]){$outputfile=$ARGV[1];}else{print "Enter output file name\n"; $outputfile=<STDIN>;}
chomp($outputfile);
open OUT,">$outputfile" or die "cannot create $outputfile.....\n\n";


#Read database sequence in to a hash %seq_hash{}, remove ">" from header. remove any white spaces and newline from sequence.


print "reading Sequences from input file.....Plz wait...\n";
#use RS_Subs ();
(my$rseq,my$rconhash,my$revconhash)=ReadFasta($inputfile);
my%seq=%$rseq;
my%conhash=%$rconhash;
my%revconhash=%$revconhash;
foreach my $head (keys %seq) {
    print OUT">$head $revconhash{$head}\n$seq{$head}\n";
	$count++;
}#
close OUT;

open OUT2,">$outputfile.Name_conversion.table" or die "cannot create $outputfile.Name_conversion.table.....\n\n";
foreach my $chead (sort keys %conhash) {
    print OUT2"$conhash{$chead}\t$chead\n";
}#
close OUT2;


#print "Done....\n$count sequences cleaned and written in the file $outputfile\n\n";
#@seq_count=();


###############################
sub ReadFasta{ # to read fasta format files into hash. returns hash.

	my $seqfile=shift(@_);


	#my ($header,@sequence,$oheader);
	chomp $seqfile;
	open FASTA,"$seqfile";
	print "reading Sequences from input file.....Plz wait...\n";
	my%seq_hash;
	#$seq_hash{'RS_Concatenated'}="";
  my%con_hash; ## to collect seq name conversion info.
  my%rev_con_hash; ## to collect reverse of conhash
	$/="\n>";    # Change record seperator to read Fasta

	while(<FASTA>){
      my($oheader,@sequence);
      my$last_N=1;
    	chomp;
    	($oheader,@sequence)=split("\n",$_);
      next if ! $oheader;
    	$oheader=~s/>//;						# Remove Leading > from Header
    	$oheader=~s/\s*$//;					# Remove trailing spaces from header
    	$oheader=~s/^\s*//;					# Remove Leading spaces from Header

      my@header=split(/\s+/,$oheader);
      my$header=$header[0];
    	my$sequence= join("",@sequence);
    	$sequence=~ s/\s//g;
    	$sequence=~s/\n//g;
      $sequence=~s/\*//g;


    	if($header=~/^\s*$/){next;}

    	# Make sure no duplicate header exists else they will be overwritten and lost during hashing due to same key.
    	if(!exists $seq_hash{$header}){
    		$seq_hash{$header}=$sequence;     #feed headers and sequences in hash.
			#$seq_hash{'RS_Concatenated'}=$seq_hash{'RS_Concatenated'}.$Ns.$sequence;
		  }
		  else {

  			# find a uniq header name by adding a number at the end. If header still exists, increase the number by one
        my$theader=$header;
  			while(exists $seq_hash{$theader}){$theader="$header\_$last_N";$last_N++;}
        $header=$theader;
  			$seq_hash{$header}=$sequence;
     		#$seq_hash{'RS_Concatenated'}=$seq_hash{'RS_Concatenated'}.$Ns.$sequence;
  		}

    $con_hash{$oheader}=$header;
    $rev_con_hash{$header}=$oheader;

	}

	my @seq_count=keys (%seq_hash);
	my $seq_count=@seq_count;

	print "Done....\nNumber of sequences read form input file = $seq_count\n\n";
	@seq_count=();
	$/="\n";    # Record seperator set back to default newline.

	return(\%seq_hash,\%con_hash,\%rev_con_hash);

}
