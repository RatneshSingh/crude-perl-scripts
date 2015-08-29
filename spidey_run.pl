use strict;
use warnings;


print "USAGE:\nscript geneseqfile mRNAsequence mRNAtable spidey_outputfile\n";
my (%seq_hash,%mRNA_hash,%genpair);

my$seqfile=$ARGV[0];
my$mRNA_seq=$ARGV[1];
my$mRNA_table=$ARGV[2];
my$spidey_output=$ARGV[3];


open FASTA,"$seqfile" or die "Cannot open $seqfile";
#open GENE,">$gene_output" or die "Cannot open $gene_output\n ";
open MRNA,"$mRNA_seq" or die "Cannot open $mRNA_seq\n ";
open TABLE,"$mRNA_table" or die "Cannot open $mRNA_table\n ";
#open SUMMARY,">$summary_output" or die "Cannot open $summary_output\n ";


####################################################################################################
# Read database sequence in to a hash %seq_hash{}, remove ">" from header. 
# remove any white spaces and newline from sequence.

print "reading Sequences from input file.....Plz wait...\n";
$/="\n>";

while(<FASTA>){#
    chomp;
    my($header,@sequence)=split("\n",$_);
    $header=~s/>//;
    $header=~s/\s//g;
    my$sequence= join("",@sequence);
    $sequence=~ s/\s//g;
    $sequence=~s/\n//g;
    #print ">$header\n";
    #feed headers and sequences in hash.
    $seq_hash{$header}=$sequence;

    #print "$seq_hash{$header}\n\n";
}#
my @seq_count=keys (%seq_hash);
my $seq_count=@seq_count;

print "Done....\nNumber of sequences read form input file = $seq_count\n\n";
@seq_count=();
############################################################################################

$/="\n>";
while(<MRNA>){
		chomp;
    my($mheader,@msequence)=split("\n",$_);
    $mheader=~s/>//;
    $mheader=~s/\s//g;

    my$msequence= join("",@msequence);
    $msequence=~ s/\s//g;
    $msequence=~s/\n//g;
    #print ">$mheader\n";
    #feed headers and sequences in hash.
    chomp($mheader,$msequence);
    $mRNA_hash{$mheader}=$msequence;

    #print "$mRNA_hash{$mheader}\n\n";
}#
my @mRNA_count=keys(%mRNA_hash);
my $mRNA_count=@mRNA_count;

print "Done....\nNumber of mRNA sequences read form input file = $mRNA_count\n\n";
@mRNA_count=0;

###########################################################################################################
# Prepare to run spidey. Arguments needs mRNA file name containing all the mRNA and  a table file with mRNA
# name (column1) and related genomic DNA name (column2). this script will read the table and 
# and group the mRNAs from same genomic DNA in one and run spidey over them.
# 

$/="\n";

while(<TABLE>){
	my($mRNAname,$genomic_name)= split(/\t+/,$_);
	chomp($mRNAname,$genomic_name);
	$mRNAname=~ s/\s//g;
	$genomic_name=~ s/\s//g;
	$genpair{$mRNAname}=$genomic_name;
	
}

foreach my$name(keys(%genpair)){

	chomp (my$mRNA_name=$name);
	open MRNAFILE,">$mRNA_name.tmp";
	print MRNAFILE">$mRNA_name\n$mRNA_hash{$mRNA_name}";
	#print ">$mRNA_name\n$mRNA_hash{$mRNA_name}";

	open GENENAME,">$genpair{$mRNA_name}.tmp";
	my $newGeneName=$genpair{$mRNA_name};
	print GENENAME">$genpair{$mRNA_name}\n$seq_hash{$newGeneName}";
	#print ">$genpair{$mRNA_name}\n$seq_hash{$newGeneName}";

	my$spidy_command= 'spidey.linux -i '."$genpair{$mRNA_name}.tmp".' -m '."$mRNA_name.tmp".' -r p '.">>$spidey_output";
	
	print "\nsending command\n$spidy_command\n";
	`$spidy_command`;
	
	my$deletecommand='rm *.tmp';
	`$deletecommand`;
	
	
	
	}

