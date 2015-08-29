#!/usr/bin/perl -w
use strict;
use Getopt::Std;

our ($opt_s,$opt_i,$opt_l,$opt_o,$opt_c,$opt_g);
getopt('silocg');

my$usage="Trims fasta names to user defined length. user can choose a fraction of name as sequence name\n
usage; perl script options
-s	Sequence file in fasta format
-i	start site to trim from (deafult=0)
-l	length of trimmed name (default=30)
-o	output_file to save trimmed fasta (default:NTrim_sequencefile.fasta)
-c	convert ambigous bases to 'N' (default: False)\n
-g	remove 'gi|' from names
";

print "$usage\n";
die "Cannot find sequence file\n" if !defined $opt_s;

$opt_i=0 if !defined $opt_i;
$opt_l=30 if !defined $opt_l;
$opt_o= "NTrim_".$opt_s if !defined $opt_o;

my%fasta=ReadFasta($opt_s); 


open OUT,">$opt_o" or die "Cannot open output file\n$usage";

foreach(keys%fasta){

	my$short_fasta2=$_;
	if($opt_g){$short_fasta2=$_;$short_fasta2=~ s/_gi\|/_/g;$short_fasta2=~ s/\s+gi\|/_/g;$short_fasta2=~ s/\|/_/g;}

	my$short_fasta=substr($short_fasta2,$opt_i,$opt_l);

	if($opt_c){	$fasta{$_}=~s/[^ATGCnatgcn]/n/g;
				#$fasta{$_}=~s/[^atgcn]/n/g;
	}

	print OUT">$short_fasta\n$fasta{$_}\n";
}

close OUT;




###############################################################################
sub ReadFasta{ # to read fasta format files into hash. returns hash.
	
	my $seqfile=shift(@_);
	
	
	my ($header,@sequence);
	chomp $seqfile;
	open FASTA,"$seqfile" or die "\nCannot open Sequence file\n";
	print "reading Sequences from input file.....Plz wait...\n";
	my%seq_hash=();
	$/="\n>";    # Change record seperator to read Fasta

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
    	$seq_hash{$header}=$sequence;     #feed headers and sequences in hash.
		
	}
	
	my @seq_count=keys (%seq_hash);
	my $seq_count=@seq_count;

	print "Done....\nNumber of sequences read form input file = $seq_count\n\n";
	@seq_count=();
	$/="\n";    # Record seperator set back to default newline.
	close FASTA;
	return(%seq_hash);

}
#-------------------------------------End ReadFasta---------------------------------------+

#############################################################################################################
