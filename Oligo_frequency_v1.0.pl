#!/usr/bin/perl -w
use strict;
use Getopt::Std;

our($opt_f,$opt_l,%count,$opt_s,%seq_hash,$opt_o,$opt_c,$out);
getopt('slfco');



my$usage="

This scripts counts the number of all possible combination of n-mer oligo in a given sequence.

usage: perl script -f/-s [options]

-f	File containing sequences in fasta format or
-s	Sequence provided on the commandline.
-l	length of n-mer to be searched [4].
-o	Outputfile to save results [STDOUT].
-c	Concatenate multiple sequences and then count [FALSE].

";





if (!$opt_f){ if(@ARGV){$opt_f=$ARGV[0]}}
if (!$opt_l){$opt_l=4;}
if(!$opt_c){$opt_c='FALSE'}
elsif($opt_c =~ m/TRUE/i){$opt_c='TRUE'};


#Reading sequences in a hash from file or commandline.
if (defined $opt_f){%seq_hash=ReadFasta($opt_f);}
elsif(defined $opt_s){$seq_hash{'New Sequence'}=$opt_s;}
else{ die "Please provide sequence or the name of file containing sequences";}




# concatenate the sequences if asked to do so.
if($opt_c eq 'TRUE'){
	print "\nUser has asked to concatenate all the sequences to one before counting $opt_l mers.\n";
	my$newseq="";
	
	foreach my$name(keys %seq_hash){$newseq=$newseq.$seq_hash{$name};}
	$newseq=~s/\s+//;
	undef %seq_hash;
	$seq_hash{'NewSeq'}=$newseq;
	
	#%seq_hash=%seq_hash1;
}



# Create all possible oligos of a specified length using kmer subroutine copied from web.
# http://www.bioperl.org/wiki/Getting_all_k-mer_combinations_of_residues


die "positive integer needed as the argument!" unless $opt_l > 0 and $opt_l =~ /^\d$/;

print "Creating $opt_l mer table......."; 
my@words=create_mer_table($opt_l);
print "..Done\n";


if($opt_o){ open($out,'>',$opt_o);}
else{ $out=\*STDOUT;}




# count the occurences of each mer in sequences
print "Searching for matches to the $opt_l mers in the sequences\n";
foreach my$seq(keys %seq_hash){
	print "\tSearching in Sequence $seq\n---------------------------------------------------------------\n";
	foreach my$mer(@words){
		$count{$mer} = () = $seq_hash{$seq} =~ /$mer/ig;
		if(!$count{$mer}){$count{$mer}=0;}
	print $out "$seq\t$mer\t$count{$mer}\n" if $count{$mer} !=0;
	if($opt_o){	print "$seq\t$mer\t$count{$mer}\n" if $count{$mer} !=0;} 
	}
	print "....Done\n"
}

#foreach my$mer(keys %count){print "$seq\t$mer\t$count{$mer}\n";}


#print join("\n",@words);
print "\n";


##########################################################################################
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

##------------------------------------------------------------------------------------------------
sub kmers ($;$) {
     my $k = shift;
     my $alphabet = shift || [ qw( A T G C ) ];
     my @bases = @$alphabet;
     my @words = @bases;
     for ( 1 .. --$k ) {
         my @newwords;
         foreach my $w (@words) {
             foreach my $b (@bases) {
                 push (@newwords, $w.$b);
             }
         }
         @words = @newwords;
     }
     return @words;
 }

sub create_mer_table{
my$k=shift;
$k=~s/\D+//;
my@words=kmers($k);
return @words

}