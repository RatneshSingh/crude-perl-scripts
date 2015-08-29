#!/usr/bin/perl -w
use strict;
print '
####################################################################################
# THIS SCRIPT IS WRITEEN FOR RUNNING GENEKARKHMM ON SET OF MANY SEQUENCES. IT WILL #
# READ SEQUENCES FROM A GIVEN FASTA FILE AND RUN EACH SEQUENCE THROUGH GENSCAN     #
# ONE BY ONE. THE OUT PUT WILL BE SENT TO A FILE SUGGESTED BY THE USER.            #
#                                                                                  #
#                                                                                  #
# AUTHOR: RATNESH SINGH                                                            # 
# VERSION:1.1                                                                      #
# CONTACT: ratnesh@hawaii.edu                                                      #
# Syntax : perl Script <sequence file>                                             #
####################################################################################
';
my $fasta=();
if(!$ARGV[0]){ print "give file name contatining sequences\n"; $fasta=<STDIN>;}
else{$fasta=$ARGV[0];}
open FASTA,"$fasta" or die "Cant open $ARGV[0]";

# change location of genscan reference file

#my $location= "/usr/local/lib/GENSCAN/Arabidopsis.smat"; 
my $location= "/home/carica/softwares/genemark_hmm_euk.linux_64/o_sativa.mod"; 
print "This is location of ref file : \n $location\n";	


#read fasta file and run through genscan
	
	$/="\n>";
while(<FASTA>){
	chomp;
	my($name,@sequence)=split(/\n/,$_);
	$name=~s/>//;
#	print "Nameis : $name \n";
#	print "Sequence is :@sequence \n";
	my $sequence=join("",@sequence);
#	print "Joined sequence is :$sequence";
	print "Running GenemarkHMM for : $name\n";
	genscan($name,$sequence);

}
exit;



#################################################################################################
#subroutines for handling genscan with perl

sub genscan{
    #subroutine for handling 'genscan' from with a perl script.
    #requires genscan to be installed
    my ($name,$seq) = @_;
    open(TMP,">SeqGenemarkHMM.tmp.$name");
    print TMP ">$name\n$seq\n";
    close(TMP);
    #gmhmme3 -m ~/softwares/genemark_hmm_euk.linux_64/o_sativa.mod -p -b F153_GenemarkHMM.out -f gff3 ../Test_F153_genome.fasta
my @command=("gmhmme3","-m $location","-p -f gff3 SeqGenemarkHMM.tmp.$name  -s $name");
    system("@command >> GenemarkHMM_run.log");# or die "Cannot execute command";
    system ("rm SeqGenemarkHMM.tmp.$name");
}    
