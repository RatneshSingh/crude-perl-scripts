#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;


getopt('sto');
our ($opt_s,$opt_o,$opt_t);

# like the shell getopt, "d:" means d takes an argument
print "\n\nParameter provided\n-sequence file: $opt_s\n" if defined $opt_s;
print "-Conversion table: $opt_t\n" if defined $opt_t;
print "-print output as: $opt_o\n" if defined $opt_o;
print "Unprocessed by Getopt::Std:\n" if $ARGV[0];
foreach (@ARGV) {
  print "-- $_\n";
}



my $help="\n\nUsage options:- \n\nperl\tscript_name\t-s [file]\t-t [file]\t-o [filename]\n\n-s\tSequence file\n-t\tTable containing information for conversion in two columns\n\t(Col1:Tag to add\tCol2: Name to search for)\n-o\tOutputfile\n\n";


print"\n\n\n\nThis script uses user provided conversion table and add tcorresponding tags to the header of the sequences in aspecified database\n\n\n";
#############################################################################################
my (@wholefile,%new,%seq_hash,$list,$seqfile,$outfile);

$list=$opt_t if defined($opt_t) or die"\nThere is no file defined as Conversion File/table\n$help";
chomp ($list)if defined($list);

$seqfile=$opt_s if defined($opt_s)or die"\nThere is no file defined as sequence file\n$help";;
#else{print "Enter Sequence file containing sequences in FASTA format\n"; $seqfile=<STDIN>;}
chomp ($seqfile) if defined ($seqfile);

$outfile=$opt_o if defined($opt_o)or die"\nThere is no file defined as output file\n$help";;
#else{print "Enter the name of out file to save selected sequences\n"; $outfile=<STDIN>;}
chomp ($outfile) if defined ($outfile);



############################################################################################
if(defined($opt_t)){
open (CODE,"$list") or die "Cannot find $list\n";

print "-----Generating conversion table-----\nOriginal Name\tConverted Name\n";
while(<CODE>){
	chomp $_;
	if(/^\s*$/){next;}
my	@line=split(/\s+/,$_);
my	$code=$line[0];
	$code=~s/=/\=/g;
my	$realname=$line[1];
	$new{$realname}= $code;
	
	print "$new{$realname}\t$realname\n";
	

			}
}



##############################################################################################
open SEQ,"$seqfile" or die "Cannot find $seqfile\n";
print "\n\nReading Fasta Sequences from:: $seqfile.....Plz wait...\n";
$/="\n>";

while(<SEQ>){#
    chomp;
   my ($header,@sequence)=split("\n",$_);
    
    $header=~s/>//;
#    $header=~s/-start[\d]*.*$//;
    

 #S	print "$header\n";
    
    my $sequence= join("",@sequence);
    $sequence=~ s/\s//g;
    $sequence=~s/\n//g;
    
  			
	foreach my $key(keys(%new)){
			if($header=~m/$key/i){
			my $newkey= $new{$key}.'-'.$key;
			$header=~s/$key/$newkey/i; 
								}
							}
$seq_hash{$header}=$sequence;
			}




my @seq_count=keys (%seq_hash);
my $seq_count=@seq_count;


print "Done....\nNumber of sequences read form input file = $seq_count\n\n";

@seq_count=();



open OUT,">$outfile";
foreach my $seq_key(keys(%seq_hash)){
		print "$seq_key\n";
		print OUT">$seq_key\n$seq_hash{$seq_key}\n";
}


						
###############################################################################################						
close OUT;
close SEQ;
close CODE;