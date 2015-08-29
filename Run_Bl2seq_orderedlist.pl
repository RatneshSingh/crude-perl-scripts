#!/usr/bin/perl -w
use strict;

#-----------------------------------------------------------------+
my$usage="Please use following options with script\n\n\nUsage: \nperl	Script	BAC.list	BAC.fasta -D 1 -p blastn -e 1.0e-180\n\n\n";



#------------------------------------------------------+
# Read list of names of sequences in an array in order	
#------------------------------------------------------+
my$listfile;
if($ARGV[0]){$listfile=$ARGV[0];} else{print"$usage\nProvide the name of list file\n";$listfile=<STDIN>;}

my$Log_seq2blast='Log_seq2blast.txt';
open OUT2,">$Log_seq2blast";
	
chomp $listfile;
open LIST,"$listfile" or die "$usage\nNo list file found\nStopping script\n";
my@list;
while(<LIST>){
	chomp $_;	
	$_=~s/>//;
	$_=~s/\s*$//;
	$_=~s/^\s*//;
	if($_ eq /^\s*$/){next;}
	push(@list,$_);

	}

my$list_length=@list;
print "\nTotal number of Names in list are :$list_length\n";



#---------------------------------------------------------------------+
# Read fasta sequences in hash to use them later for blast2n
#---------------------------------------------------------------------+
my$seqfile;
if($ARGV[1]){$seqfile=$ARGV[1];} else{print"$usage\nProvide the name of sequence file\n";$seqfile=<STDIN>;}
	
my %seq_hash= ReadFasta($seqfile);


#my@seq_names= (keys(%seq_hash));				# to check the names of sequences read.
#my$seq_names=join("\n",@seq_names);
#print"$seq_names\n";



#---------------------------------------------------------------------+
# Read list and send sequences for bl2seq in pairs
#---------------------------------------------------------------------+
shift@ARGV;
shift@ARGV;
my$options=join(" ",@ARGV);
print "\noptions obtained:$options\n";



# -----use this section if names in the list and sequences are same----------------------------+

#for(my$i=0;$i<$list_length;$i++){
#	
#	print"\nChecking if $list[$i] and $list[$i+1] exist\n";
#	if (!defined $seq_hash{$list[$i]}){print"\nSequence:$list[$i] not found"; next;}
#	if (!defined $seq_hash{$list[$i+1]}){print"\nSequence:$list[$i+1] not found"; next;}
#
#	Run_bl2seq($list[$i],$seq_hash{$list[$i]},$list[$i+1],$seq_hash{$list[$i+1]},$options);
#
#}
#-----------------------------------------------------------------------------------------------+


#--------Use this section if Names in the list sequences are not exactly same-------------------+
for(my$i=0;$i-1<$list_length;$i++){
	
	if(defined $list[$i+1]){
		my$head1=();
		my$head2=();
	
		foreach my$key(keys%seq_hash){
			if($key=~/$list[$i]/){$head1=$key;}
			if($key=~/$list[$i+1]/){$head2=$key;}
		}
	
		print"\nChecking if $list[$i] and $list[$i+1] exist\t";
		#print OUT2"\nDoes $list[$i] and $list[$i+1] exist\t";
		if (!defined $head1){print"----->No\nSequence:$list[$i] not found";print OUT2"\n----->No:$list[$i] not found\n"; next;}
		elsif (!defined $head2){print"----->No\nSequence:$list[$i+1] not found"; print OUT2"\n----->No:$list[$i+1] not found\n";next;}
	
		else{
			print "----->Yes\nRunning bl2seq for $head1 and $head2\n";
			print OUT2"$list[$i] _ $list[$i+1]\.\.";
			Run_bl2seq($list[$i],$seq_hash{$head1},$list[$i+1],$seq_hash{$head2},$options);
		}
	}
	
#------------------------------------------------------------------------------------------------+

}

print"\n";

exit;






####################################################################################
# Subroutine Area
###################################################################################

#-------------------------------------Start ReadFasta---------------------------------------+

sub ReadFasta{ # to read fasta format files into hash. returns hash.
	
	my $seqfile=shift(@_);
	
	
	my ($header,@sequence);
	chomp $seqfile;
	open FASTA,"$seqfile";
	print "reading Sequences from input file.....Plz wait...\n";
	
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

	return(%seq_hash);

}
#-------------------------------------End ReadFasta---------------------------------------+


#------------------------------------Start Run_bl2seq-------------------------------------+
sub Run_bl2seq{

my$head1=();
my$seq1=();
my$head2=();
my$seq2=();
my$options=();

$head1=shift@_;
$seq1=shift@_;
$head2=shift@_;
$seq2=shift@_;
$options=shift@_;

chomp($head1,$head2,$seq1,$seq2);

#if($head1 or $head2 or $seq1 or $seq2 eq (^\s*$)){next;}


open SEQ1,">$head1.tmp";
print SEQ1">$head1\n$seq1";
open SEQ2,">$head2.tmp";
print SEQ2">$head2\n$seq2";

my$output=$head1.'-Vs-'.$head2.'.bl2seq';

my$command=join(" ",'bl2seq','-i',"$head1.tmp",'-j',"$head2.tmp",$options,'-o',$output);

system($command)==0 or die "system $command failed: $?\n";

my$comb_seqfile=$head1.'-'.$head2.'.fasta';

open SEQCOMB,">$comb_seqfile";					#adding two seq in one file for ease of comparision.
print SEQCOMB">$head1\n$seq1\n>$head2\n$seq2"; #adding two seq in one file for ease of comparision.

#system("rm $head1.tmp");
#system("rm $head2.tmp");

$head1=();
$seq1=();
$head2=();
$seq2=();
$options=();

return();

}
#-----------------------------------End Run_bl2seq----------------------------------------+