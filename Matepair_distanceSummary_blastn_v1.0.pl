#usr/bin/perl -w


#=============================================================================
# This Script Reads blast results in table format run with -v 1 -b 1 options
# Use All the reads as query and assembled sequences as Database to generate 
# reliable results. This Script will read blast result and calculate distances
# between forward and reverse pairs of Sequencing reads. Forward and Reverse 
# are identified by suffix .f1 and .r1 at the ends. If using any other suffix
# need to change it in the script. 
#
# Author: Ratnesh Singh
# Version: 1.0  
# Report any Bug to: ratnesh@hawaii.edu
#=============================================================================


use strict;

my ($blastfile,$queryfile,%reads,@reads);

if (defined $ARGV[0]){$blastfile=$ARGV[0];} else {print "\nEnter file containing blast results in table format,\nwhere reads were used as query against database made of Assembled sequences\nwith options '-v 1 -b 1'\n";$blastfile=<STDIN>;}
if (defined $ARGV[1]){$queryfile=$ARGV[1];} else {print "\nEnter file containing query sequences in FASTA format,\n";$queryfile=<STDIN>;}


chomp $blastfile,$queryfile;

if (-e $blastfile){open BLASTN,"<$blastfile";} else {die "\n$blastfile does not exist. Please check the file name\n"; }
if (-e $queryfile){open QUERY,"<$queryfile";} else {die "\n$queryfile does not exist. Please check the file name\n"; }

open OUT,">$queryfile.$blastfile._MatePairSummary.txt";



#--------------Read blastn results and assign start and orientation of reads on assembly.------------------------------------

while(<BLASTN>){
	
	
	
	my (@line)= split(/\s+/,$_); 
	
	#------Each element of @line will contain following--------------------
	# Query= 0,			Subject=1,	
	# %Identity= 2,		Alnlength=3,	
	# Mismatch=4 ,		Gap=5,	
	# Qstart= 6,		Q_end=7,
	# S_start= 8,		S_end=9,
	# Evalue= 11,		Bitscore=12;
	#-----------------------------------------------------------------------
	
	if($reads{$line[0]}{'alnlnth'}>$line[3]){next;}
	
	else{	
	
	$line[0]=~s/\s*$//;
	$reads{$line[0]}{'contig'}=$line[1];
	$reads{$line[0]}{'qstart'}=$line[6];
	$reads{$line[0]}{'qend'}=$line[7];
	$reads{$line[0]}{'sstart'}=$line[8];
	$reads{$line[0]}{'send'}=$line[9];
	$reads{$line[0]}{'alnlnth'}=$line[3];
	if($line[9]-$line[8]>=0){$reads{$line[0]}{'orientation'}='plus';} else {$reads{$line[0]}{'orientation'}='minus';}
	}

}

#-------------------read names of reads and record them in hash---------------------------------

$/ = "\n>";

while(<QUERY>){
	my($header,@sequence)=split(/\n/,$_);
	$header=~ s/>//;
	$header=~s/\.r1$|\.f1$//;

	push(@reads,$header);

}
# ---------------set default record seperator back to normal.---------------------



$/="\n";



#---------collect uniques from reads list to array @reads2----------------------
my %hsh;
undef @hsh{@reads};
my @reads2 = keys %hsh;


#-----------------find mate pairs--------------------------------------------------------------

print OUT"Query\tF_on_contig\tR_on_contig\tF_start\tR_start\tF_orientation\tR_orientation\tDistance\n";


foreach my $read2(@reads2){

	my $F_read2=$read2.'.f1';
	my $R_read2=$read2.'.r1';
	my $distance=();
	if($reads{$F_read2}{'contig'} eq $reads{$R_read2}{'contig'} ){ $distance= $reads{$F_read2}{'sstart'}-$reads{$R_read2}{'sstart'};}
	elsif (!defined $reads{$F_read2}{'contig'}) {$distance= 'F_read missing';}
	elsif (!defined $reads{$R_read2}{'contig'}) {$distance= 'R_read missing';}
	else{$distance= 'Not on same contig';}

	print OUT"$read2\t$reads{$F_read2}{'contig'}\t$reads{$R_read2}{'contig'}\t$reads{$F_read2}{'sstart'}\t$reads{$R_read2}{'sstart'}\t$reads{$F_read2}{'orientation'}\t$reads{$R_read2}{'orientation'}\t$distance\n";
	print "$read2\t$reads{$F_read2}{'contig'}\t$reads{$R_read2}{'contig'}\t$reads{$F_read2}{'sstart'}\t$reads{$R_read2}{'sstart'}\t$reads{$F_read2}{'orientation'}\t$reads{$R_read2}{'orientation'}\t$distance\n";


}


