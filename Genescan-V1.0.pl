#!/usr/bin/perl -w
use strict;
################################################################################
# THIS SCRIPT IS WRITEEN FOR RUNNING GENSCAN ON SET OF MANY SEQUENCES. IT WILL #
# READ SEQUENCES FROM A GIVEN FASTA FILE AND RUN EACH SEQUENCE THROUGH GENSCAN #
# ONE BY ONE. THE OUT PUT WILL BE SENT TO A FILE SUGGESTED BY THE USER.        #
#                                                                              #
#                                                                              #
# AUTHOR: RATNESH SINGH                                                        # 
# VERSION:1.0                                                                  #
# CONTACT: ratnesh@hawaii.edu                                                  #
# genscan parfname seqfname [-v] [-cds] [-subopt cutoff] [-ps psfname scale]   #
################################################################################

open FASTA,"$ARGV[0]" or die "Cant open $ARGV[1]";
	


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
	genscan($name,$sequence);

}
exit;






sub genscan{
    #subroutine for handling 'genscan' from with a perl script.
    #requires genscan to be installed
    #returns stucture in parentheses format and predicted energy of structure
    my ($name,$seq) = @_;
    
    my ($fold,$energy,@lines) = ();
    open(TMP,">SeqGenscan.tmp.$name");
    print TMP ">$name\n$seq\n";
    close(TMP);
# 	system("touch $name.genscan_ORF");
#   print "genscan  SeqGenscan.tmp.$name >> genscan_output";
my @command=("genscan","/usr/local/lib/GENSCAN/Arabidopsis.smat","SeqGenscan.tmp.$name");
    system("@command >> genscan_$name")==0 or die "Cannot execute command";
    
    system ("rm SeqGenscan.tmp.$name");
#    return($abc);
}    
