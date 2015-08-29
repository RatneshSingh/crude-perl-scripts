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
print "\n
		 THIS SCRIPT IS WRITEEN FOR RUNNING HMMFETCH ON LIST OF HMM MODEL NAMES.
		 IT WILL READ NAMES FROM A GIVEN LIST FILE AND RUN EACH NAME THROUGH 
		 HMMFETCH ONE BY ONE. THE OUT PUT WILL BE SENT TO A FILE NAMED AFTER MODEL. 
		 SYNTAX : SCRIPT HMM_MODEL LIST_DATABASE\n\n";

if (!$ARGV[0]){
	print "Give the name of list containing name of hmm models to pick\n";
	$ARGV[0]=<STDIN>;
}
open LIST,"$ARGV[0]" or die "Cant open $ARGV[0]";
	

if (!$ARGV[1]){
	print "Give the full path for hmm database:\n";
	$ARGV[1]=<STDIN>;
}
my $hmm_database=$ARGV[1];
chomp($hmm_database);

#my $hmm_database="/home/databases_misc/Pfam_fs";
print "Database location provided: $hmm_database \n";
if(-e $hmm_database){
	

	#read fasta file and run through genscan
	
	
	while(<LIST>){
		chomp;
		my $name=$_;
	
#	print "Nameis : $name \n";
	print "\n\nFetching hmm model for: $name\n";
	hmmfetch($name,$hmm_database);

	}
}

else{
	print "cant find $hmm_database \n check if it exist in the given location.\n";
}
close(LIST);
exit;



#################################################################################################
#subroutine for handling 'hmmfetch' from with a perl script.
#requires hmmer to be installed


sub hmmfetch{
        my ($name,$hmm_database) = @_;
 
		my @command=("hmmfetch","$hmm_database","$name");
		#print "command given : @command \n";
    	if(system("@command >> $name.hmm")==0){
    	print "Done fecthing for: $name\n"	
    	} 
    	else{
    		if(-z "$name.hmm"){
    			
    		system("rm $name.hmm");
    		}
 		print "cant find any Hmm model for name: $name\n ";
 		}
}    
