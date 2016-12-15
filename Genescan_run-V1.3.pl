#!/usr/bin/perl -w
use strict;
use Getopt::Std;
print '
#################################################################################
# THIS SCRIPT IS WRITEEN FOR RUNNING GENSCAN ON SET OF MANY SEQUENCES. IT WILL  #
# READ SEQUENCES FROM A GIVEN FASTA FILE AND RUN EACH SEQUENCE THROUGH GENSCAN  #
# ONE BY ONE. THE OUT PUT WILL BE SENT TO A FILE SUGGESTED BY THE USER.         #
#                                                                               #
#                                                                               #
# AUTHOR: RATNESH SINGH                                                         # 
# VERSION:1.2                                                                   #
# CONTACT: ratnesh@hawaii.edu                                                   #
# Syntax : perl Script options..                                                #
#  -s	sequence file							                                #
#  -t	arabidopsis|maize|Human|FullPath to .smat file. use trained model of	#
#		species for gene search.  												#
#        default [Arabidopsis]. Expects these to be in /usr/local/genscan.      #
#  -c	yes|no. print cds for each file [yes]. 				       			    #
#  -p   Num CPU to use[1].	will require ParaFly|| parallel to be installed.	#
#  -l   Split fasta if longer than this Num [2000000]							#
#  -o   Sequence overlap while splitting fasta [0]   							#
#################################################################################
';
our($opt_s,$opt_t,$opt_c,$opt_p,$opt_l,$opt_o);
$opt_c='yes';
$opt_t='arabidopsis';
$opt_p=1;
$opt_l=2000000;
$opt_o=0;
getopt('stcplo');

#Define fasta file containing sequence in fasta format to be used in Genscan
my $fasta=();
if($opt_s){$fasta=$opt_s}
else{ die "Fasta file name contatining sequences is not provided. On what I am going to run Genscan on?\n";}

$fasta=~s/\s+//g;
open FASTA,"$fasta" or die "Cant open $opt_s";

my $location=$opt_t;
if (lc $opt_t eq 'arabidopsis'){ $location= "/usr/local/genscan/Arabidopsis.smat";}
elsif(lc $opt_t eq 'maize'){$location= "/usr/local/genscan/Maize.smat";}
elsif(lc $opt_t eq 'human'){$location= "/usr/local/genscan/HumanIso.smat";}

my$split_len=$opt_l;
my $CPU=$opt_p;
$opt_c=" " if lc$opt_c eq 'no';
$opt_c="-cds" if lc$opt_c eq 'yes';

print "\nGenscan will be using training data for: $location\n";	
$location=~s/\s+//g;

print "This is location of ref file : \n $location\n";	

$opt_o < $opt_l or die "Overlap ($opt_o) cannot be equal or larger that splitting length ($opt_l) threshold. Doing so will enter program into infinite loop.\n. exiting.";


#read fasta file and run through genscan
	
	$/="\n>";
	my@seq_list;
	my %seq;
while(<FASTA>){
	chomp;
	my($name,@sequence)=split(/\n/,$_);
	$name=~s/>//;
	my $sequence=join("",@sequence);
	print "Reading Sequence $name\n";
	
	if (length($sequence) > $split_len) {
		print "$name is larger than $split_len . Splitting in to smaller chunks.\n";
        split_fasta(\$name,\$sequence,\$split_len,\%seq,\$opt_o);
    }else{$seq{$name}=$sequence;}
}

open("CMD",">genscan.cmd");
foreach my$seqname(keys %seq){
	if ($CPU > 1) {
    	open(TMP,">SeqGenscan.tmp.$seqname"); print TMP ">$seqname\n$seq{$seqname}\n"; close(TMP);
		my @command=("genscan","$location","SeqGenscan.tmp.$seqname","$opt_c");
		print CMD join " ",@command," > genscan_$seqname 2> genscan_$seqname.log","\n";
	}else{
		print "Running Genscan for: $seqname\n";
		genscan($seqname,$seq{$seqname},$opt_c);
	}
}

if ($CPU > 1) {
	print "\n Running genscan command in $CPU parralel threads\n";
    system("ParaFly -c genscan.cmd -failed_cmds genscan.cmd.failed -v -shuffle -CPU $CPU") || system("parallel -j $CPU < genscan.cmd") || system("sh genscan.cmd");
}



exit;



#################################################################################################
####subroutines for handling genscan with perl
sub split_fasta{
	my($refname,$refsequence,$refsplit_len,$hRef_seq,$refOverlap)=@_;
	my$seq_length=length($$refsequence);
	if ($seq_length <= $$refsplit_len) {
        $$hRef_seq{$$refname}=$$refsequence;
		return;
    }
    
	for(my$i=0;$i<=$seq_length;$i=$i + $$refsplit_len - $$refOverlap){
		my$start=$i;
		$start=0 if $start < 0;
		my$end=$i+$$refsplit_len ;
		$end=$seq_length  if $end > $seq_length ;
		
		my$subseq_len=$end - $start;
		
		my$subseq=substr($$refsequence,$i,$subseq_len);
		my$new_name=join("",$$refname,"_",$i+1,"-",$end);
		$$hRef_seq{$new_name}=$subseq;
	}
}


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


