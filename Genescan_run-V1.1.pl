#!/usr/bin/perl -w
use strict;
use Getopt::Std;
print '
################################################################################
# THIS SCRIPT IS WRITEEN FOR RUNNING GENSCAN ON SET OF MANY SEQUENCES. IT WILL #
# READ SEQUENCES FROM A GIVEN FASTA FILE AND RUN EACH SEQUENCE THROUGH GENSCAN #
# ONE BY ONE. THE OUT PUT WILL BE SENT TO A FILE SUGGESTED BY THE USER.        #
#                                                                              #
#                                                                              #
# AUTHOR: RATNESH SINGH                                                        # 
# VERSION:1.1                                                                  #
# CONTACT: ratnesh@hawaii.edu                                                  #
# Syntax : perl Script options..                                               #
#  -s	sequence file							       #
#  -t	arabidopsis|maize|Human. use trained model of species for gene search. #
#        default [Arabidopsis]                                                 #
#  -c	yes|no. print cds for each file [yes]. 				       #
#################################################################################
';
our($opt_s,$opt_t,$opt_c);
$opt_c='yes';
$opt_t='arabidopsis';
getopt('stc');

#Define fasta file containing sequence in fasta format to be used in GlimmerHMM
my $fasta=();
if($opt_s){$fasta=$opt_s}
else{ die "Fasta file name contatining sequences is not provided. On what I am going to run Genscan on?\n";}

$fasta=~s/\s+//g;
open FASTA,"$fasta" or die "Cant open $opt_s";


my $location= "/usr/local/genscan/Arabidopsis.smat" if lc $opt_t eq 'arabidopsis';; 
$location= "/usr/local/genscan/Maize.smat" if lc $opt_t eq 'maize';
$location= "/usr/local/genscan/HumanIso.smat" if lc $opt_t eq 'human';; 


$opt_c=" " if lc$opt_c eq 'no';
$opt_c="-cds" if lc$opt_c eq 'yes';

print "\nGenscan will be using training data for: $location\n";	
$location=~s/\s+//g;

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
	print "Running Genscan for : $name\n";
	genscan($name,$sequence,$opt_c);

}
exit;



#################################################################################################
#subroutines for handling genscan with perl

sub genscan{
    #subroutine for handling 'genscan' from with a perl script.
    #requires genscan to be installed
    my ($name,$seq,$cds) = @_;
    open(TMP,">SeqGenscan.tmp.$name");
    print TMP ">$name\n$seq\n";
    close(TMP);
my @command=("genscan","$location","SeqGenscan.tmp.$name","$cds");
print "sending command to the system:\n@command\n";
    system("@command > genscan_$name");# or die "Cannot execute command";
    system ("rm SeqGenscan.tmp.$name");
}    
