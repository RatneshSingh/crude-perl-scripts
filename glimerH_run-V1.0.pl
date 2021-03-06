#!/usr/bin/perl -w
use strict;
print '
###################################################################################
# THIS SCRIPT IS WRITEEN FOR RUNNING GlimmerHMM ON SET OF MANY SEQUENCES. IT WILL #
# READ SEQUENCES FROM A GIVEN FASTA FILE AND RUN EACH SEQUENCE THROUGH GENSCAN    #
# ONE BY ONE. THE OUT PUT WILL BE SENT TO A FILE SUGGESTED BY THE USER.           #
#                                                                                 #
#                                                                                 #
# AUTHOR: RATNESH SINGH                                                           # 
# VERSION:1.0                                                                     #
# CONTACT: ratnesh@hawaii.edu                                                     #
# Syntax : perl Script <sequence file>                                            #
###################################################################################
';

#Define fasta file containing sequence in fasta format to be used in GlimmerHMM
my $fasta=();
if(!$ARGV[0]){ print "give file name contatining sequences\n"; $fasta=<STDIN>;}
else{$fasta=$ARGV[0];}
open FASTA,"$fasta" or die "Cant open $ARGV[0]";

# change location of gen reference file
my $location=();
if(!$ARGV[1]){ print "give location of training data folder to use\n"; $location=<STDIN>;}
else{$location=$ARGV[1];}
print "GlimmerHMM is trained on: \n $location\n";	


#read fasta file and run through genscan
	
	$/="\n>";
while(<FASTA>){
	chomp;
	my($name,@sequence)=split(/\n/,$_);
	$name=~s/>//;
#	print "Nameis : $name \n";
#	print "Sequence is :@sequence \n";
	my $sequence=join("",@sequence);
	$sequence=~ s/\s//g;
#	print "Joined sequence is :$sequence";
	print "Running GlimmerHMM for : $name\n";
	glimmerHMM($name,$sequence);

}
exit;



#################################################################################################
#subroutines for handling genscan with perl

sub glimmerHMM{
    #subroutine for handling 'glimmerHMM' from with a perl script.
    #requires glimmer to be installed
    my ($name,$seq) = @_;
    open(TMP,">SeqGlimmerHMM.tmp.$name");
    print TMP ">$name\n$seq\n";
    close(TMP);
my @command=("glimmerhmm","SeqGlimmerHMM.tmp.$name","-d"," $location","-o","glimmerHMM.$name.out","-f");
    system("@command")==0 or die "Cannot execute command";
    system ("rm SeqGlimmerHMM.tmp.$name");
}    
