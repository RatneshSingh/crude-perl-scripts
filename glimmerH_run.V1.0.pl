#!/usr/bin/perl -w
use strict;
use Getopt::Std;
print '
#################################################################################
# THIS SCRIPT IS WRITEEN FOR RUNNING GLIMMERHMM ON SET OF MANY SEQUENCES. IT WILL  #
# READ SEQUENCES FROM A GIVEN FASTA FILE AND RUN EACH SEQUENCE THROUGH GLIMMERHMM  #
# ONE BY ONE. THE OUT PUT WILL BE SENT TO A FILE SUGGESTED BY THE USER.         #
#                                                                               #
#                                                                               #
# AUTHOR: RATNESH SINGH                                                         # 
# VERSION:1.2                                                                   #
# CONTACT: ratnesh@hawaii.edu                                                   #
# Syntax : perl Script options..                                                #
#  -s	sequence file							                                #
#  -t	arabidopsis|rice|Human|zebrafish|celegans|FullPath to .smat file. use trained model of	#
#		species for gene search.  												#
#       default [Arabidopsis]. Expects these to be in /usr/local/GlimmerHMM/trained_dir/.   #
#  -c   options to add to glimmerhmm.
#  -p   Num CPU to use[1].	will require ParaFly|| parallel to be installed.	#
#  -l   Split fasta if longer than this Num [2000000]							#
#  -o   Sequence overlap while splitting fasta [0]   							#
#################################################################################
';
our($opt_s,$opt_t,$opt_c,$opt_p,$opt_l,$opt_o);
$opt_c=' -g ';
$opt_t='arabidopsis';
$opt_p=1;
$opt_l=2000000;
$opt_o=0;
getopt('stpclo');

#Define fasta file containing sequence in fasta format to be used in Glimmerhmm
my $fasta=();
if($opt_s){$fasta=$opt_s}
else{ die "Fasta file name contatining sequences is not provided. On what I am going to run Glimmerhmm on?\n";}

$fasta=~s/\s+//g;
open FASTA,"$fasta" or die "Cant open $opt_s";

my $location=$opt_t;
if (lc $opt_t eq 'arabidopsis'){ $location= "/usr/local/GlimmerHMM/trained_dir/arabidopsis";}
elsif(lc $opt_t eq 'human'){$location= "/usr/local/GlimmerHMM/trained_dir/human";}
elsif(lc $opt_t eq 'celegans'){$location= "/usr/local/GlimmerHMM/trained_dir/Celegans";}
elsif(lc $opt_t eq 'rice'){$location= "/usr/local/GlimmerHMM/trained_dir/rice";}
elsif(lc $opt_t eq 'zebrafish'){$location= "/usr/local/GlimmerHMM/trained_dir/zebrafish";}

my$split_len=$opt_l;
my $CPU=$opt_p;
my $options=$opt_c;


print "\nGlimmerhmm will be using training data for: $location\n";	
$location=~s/\s+//g;

print "This is location of ref file : \n $location\n";	

$opt_o < $opt_l or die "Overlap ($opt_o) cannot be equal or larger that splitting length ($opt_l) threshold. Doing so will enter program into infinite loop.\n. exiting.";


#read fasta file and run through glimmerhmm
	
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

open("CMD",">glimmerhmm.cmd");
foreach my$seqname(keys %seq){
	my $seq_file=$seqname;
		$seq_file=~s/\s+\S+$//g;
	if ($CPU > 1) {
    	open(TMP,">SeqGlimmerhmm.tmp.$seq_file"); print TMP ">$seqname\n$seq{$seqname}\n"; close(TMP);
		my @command=("glimmerhmm","SeqGlimmerhmm.tmp.$seq_file","$location","$options");
		print CMD join " ",@command," > $seq_file.glimmer 2> $seq_file.glimmerhmm.log","\n";
	}else{
		print "Running Glimmerhmm for: $seq_file\n";
		glimmerHMM($seqname,$seq{$seqname},$location,$options);
	}
}

if ($CPU > 1) {
	print "\n Running glimmerhmm command in $CPU parralel threads\n";
    system("ParaFly -c glimmerhmm.cmd -failed_cmds glimmerhmm.cmd.failed -v -shuffle -CPU $CPU") || system("parallel -j $CPU < glimmerhmm.cmd") || system("sh glimmerhmm.cmd");
}



exit;



#################################################################################################
####subroutines for handling glimmerhmm with perl
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
		my$new_name=join("",$$refname,"_",$i+1,"-",$end,"   ",length($subseq),"nt");
		$$hRef_seq{$new_name}=$subseq;
	}
}

#################################################################################################
#subroutines for handling glimmerhmm with perl

sub glimmerHMM{
    #subroutine for handling 'glimmerHMM' from within a perl script.
    #requires glimmer to be installed
    my ($name,$seq,$location,$options) = @_;
	my $tmp_file=">SeqGlimmerHMM.tmp.$name";
    open(TMP,">$tmp_file");
    print TMP ">$name\n$seq\n";
    close(TMP);
my @command=("glimmerhmm "," $tmp_file "," $location ","$options", " > "," $tmp_file.glimmer ");
    system("@command")==0 or die "Cannot execute command";
    system ("rm $tmp_file");
}    
