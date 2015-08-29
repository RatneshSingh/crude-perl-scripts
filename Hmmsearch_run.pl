#!/usr/bin/perl -w
use strict;
use Getopt::Std;
my$usage= '
###############################################################################################################
# THIS SCRIPT IS WRITEEN FOR RUNNING HMMSEARCH v3 ON SET OF MODELS. IT WILL                                   #
# READ LIST OF HMM NAMES AND CUTOFF FROM A GIVEN LIST FILE AND PICKS IT                                       #
# FROM PFAM DATABASE FILE AND RUN IT THROUGH HMMSEARCH ONE BY ONE.                                            #
# THE OUT PUT WILL BE SENT TO NEW FILE FOR EACH MODEL.                                                        #
# AUTHOR: RATNESH SINGH                                                                                       # 
# VERSION:1.0                                                                                                 #
# CONTACT: ratnesh@hawaii.edu                                                                                 #
# Syntax : perl <script> -h <hmmFile> -s <sequence_path>                                                      #
#other options:                                                                                               #
# -o	Output file name.										      #
# -l	pattern_list_with_evalue									      #
# -d	Yes|no. Dont delete temporary hmm files fetched from Hmm database. [no]                               #
# -m	Hmm acc Number. Manual entry of hmm and evalue.                                                       #
# -e	E value to use [1]										      #
# -f	Yes|no. Only fetch hmm file. dont run the hmmsearch [no]  					      #
###############################################################################################################

';


our($opt_h,$opt_s,$opt_l,$opt_d,$opt_o,$opt_m,$opt_e,$opt_f);

$opt_d='no';
$opt_f='no';
#$opt_e=1;
#$opt_m='auto';
getopt('hsldomef');

 die "Hmmfile missing. \n$usage" if !$opt_h;
 die "Sequence file missing.\n$usage" if !$opt_s && $opt_l;
 $opt_d='yes' if lc$opt_f eq 'yes';
 $opt_e=1 if lc$opt_f eq 'yes';


my($seq_loc,$hmm,$cutof,%ev,@hmm);
open HMM,"$opt_h" or die "Cant open $opt_h";

#open FASTA,"$opt_s" or die "Cant open $opt_s";
$seq_loc = $opt_s if $opt_s;
chomp $seq_loc if $opt_s;

if(defined$opt_l){open LIST,"$opt_l" or die "Cant open $opt_l";}
elsif($opt_m || !$opt_l){
	print "List of hmm model names and evalues has not been provided. \n Program is working in Manual mode\nName of Hmm model:";
	my $hmm1;
	$hmm1 = $opt_m or $hmm1= <STDIN>;
	chomp $hmm1;
	print $hmm1;
	
	push(@hmm,$hmm1);
	print "\nThe evalue cutof:" if lc$opt_f ne 'yes';
	my$evalue;
	$evalue = $opt_e or $evalue = <STDIN>;
	chomp $evalue;
	print $evalue;
	$ev{$hmm1}=$evalue;
}
else {die "\n did not find hmm name and evalues\n";}

print "\n fetch only flag (-f) is turned on. Program will only fetch the hmm model and exit.\n" if lc$opt_f eq 'yes';


#read hmm name and cutof from list	
if(defined$opt_l){
	while(<LIST>){
	
		if(/\t/){($hmm,$cutof)=split(/\t/,$_);}
		else {$hmm=$_;$cutof=10;}
	
#		if($cutof eq ""){ $cutof=10;}
	
		chomp ($hmm,$cutof);
		$hmm=~s/\s//g;
		push(@hmm,$hmm);
		$ev{$hmm}=$cutof;
		#print "@hmm";
	}
}
my$filename = $seq_loc if $opt_s;
my@filename=split(/\//,$filename) if $opt_s;
$filename=$filename[-1] if $opt_s;
$filename=~s/\//_/g if $opt_s;



#read hmm file into hash.

$/="HMMER3/b [3.0 | March 2010]";

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
    my @command=();
    open(TMP,">$name.hmm");
    print TMP"HMMER3/b [3.0 | March 2010]";
    print TMP"$model";
    close(TMP);
	
    print "command being sent to system \n@command" if $opt_f eq 'no';
    if($opt_f eq 'no'){
	if($opt_o){
	@command=("hmmsearch","-E $evalue","--max","--cpu 40","--tblout hmmsearch.$opt_o.tblout","--domtblout hmmsearch.$opt_o.domtblout","$name.hmm","$seq_loc");
	system("@command >> $opt_o.hmmsearch")==0 or print "Cannot execute command";
	system ("rm $name.hmm") if lc$opt_d eq 'no';
	}
	
	else{
	@command=("hmmsearch","-E $evalue","--max","--cpu 40","--tblout hmmsearch.$name.$filename.tblout","--domtblout hmmsearch.$name.$filename.domtblout","$name.hmm","$seq_loc");
	system("@command >> $name.$filename.hmmsearch")==0 or print "Cannot execute command";
	system ("rm $name.hmm") if lc$opt_d eq 'no';
	}
    }
	print "\nSaving hmm Model as: $name.hmm......" if lc$opt_f eq 'yes';
	print "......Done\n";
}    
