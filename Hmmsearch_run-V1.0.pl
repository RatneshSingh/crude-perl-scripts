#!/usr/bin/perl -w
use strict;
print '
################################################################################
# THIS SCRIPT IS WRITEEN FOR RUNNING HMMSEARCH ON SET OF MODELS. IT WILL       #
# READ LIST OF HMM NAMES AND CUTOFF FROM A GIVEN LIST FILE AND PICKS IT        #
# FROM PFAM DATABASE FILE AND RUN IT THROUGH HMMSEARCH ONE BY ONE.             #
# THE OUT PUT WILL BE SENT TO NEW FILE FOR EACH MODEL.                         #
# AUTHOR: RATNESH SINGH                                                        # 
# VERSION:1.0                                                                  #
# CONTACT: ratnesh@hawaii.edu                                                  #
# Syntax : perl <script> <hmmFile> <sequence_path> <pattern_list_with_evalue>  #
################################################################################
';

my($seq_loc,$hmm,$cutof,%ev,@hmm,);
open HMM,"$ARGV[0]" or die "Cant open $ARGV[0]";

#open FASTA,"$ARGV[1]" or die "Cant open $ARGV[1]";
$seq_loc = $ARGV[1];
chomp $seq_loc;
open LIST,"$ARGV[2]" or die "Cant open $ARGV[2]";



#read hmm name and cutof from list	
while(<LIST>){
	
	if(/\t/){($hmm,$cutof)=split(/\s+/,$_);}
	else {$hmm=$_;$cutof=10;}
	
#	if($cutof eq ""){ $cutof=10;}
	
	chomp ($hmm,$cutof);
	$hmm=~s/\s//g;
	push(@hmm,$hmm);
	$ev{$hmm}=$cutof;
#print "@hmm";
}

#read hmm file into hash.

$/="HMMER2.0  [2.2g]";

while(<HMM>){
chomp;
#print "\n$_ \n";			
			foreach $hmm(@hmm){
#				print "Checking for $hmm \n";
				if($_=~/$hmm/){
				hmmsearch($hmm,$ev{$hmm},$_);	
				}
				else {
					next;
				}
			}			
				
}
print "\n";	
exit;


#################################################################################################
#subroutines for handling hmmsearch with perl

sub hmmsearch{
    #subroutine for handling 'hmmsearch' from with a perl script.
    #requires hmmer to be installed
    
    my ($name,$evalue,$model) = @_;
    open(TMP,">$name.hmm");
    print TMP"HMMER2.0  [2.3.2]";
    print TMP"$model";
    close(TMP);
my @command=("hmmsearch","-E $evalue","$name.hmm","$seq_loc");
    print "command being sent to system \n@command";
    system("@command >> hmmsearch.$name.txt")==0 or print "Cannot execute command";
    system ("rm $name.hmm");
	print "......Done\n";
}    
