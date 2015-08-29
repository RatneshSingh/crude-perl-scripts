#!/usr/bin/perl

#------------------------------------------------------------------------
# This script takes cis-elements sequences(-p) and count its frequency in 
# sequences(-s) provided for test. It creates 1000 or specified number of
# random sequences from provided sequences(-d) and calculates z-score which
# can be further used to calculate p-value with excel using formula 
# p-value=1-NORMSDIST(Z-score))
# This script uses 'awk' command. Runs faster but no mismatches can be used.
# to use mismatches use one depending on patmatch.
# created by : Ratnesh Singh
# In case of any bugs or suggestion contact: ratnesh@hawaii.edu
# usage: perl script -s seq -d seq -p cis-elemsnt -m number_of_datasets -n number of seq in each datasets
#------------------------------------------------------------------------




use warnings;
use strict;
use Getopt::Std;
use List::Util qw(sum);
use Statistics::Descriptive;



getopt('sdpmnrM');
my$usage="perl script {options.....}
-s 'real sequences to test in fasta'
-d 'surrogate sequences in fasta' 
-p 'pattern file containing ciselements'
-m 'number of datasets to test' [1000]
-n 'number of sequences in each dataset'
-r 'residue type (-n nucleotide,-p protein, -c complementry strand)[n]'
-M 'Number of mismatches allowed'[0]
-t 'type of mismatch: i (insertions);ids (insertion,deletions,substitutions); id (insertions,deletions)'
";


our($opt_m,$opt_n,$opt_s,$opt_p,$opt_r,$opt_M,$opt_t,$opt_d,%rseq,%seq);

print"\n\n=============================================================================\n";
print "-s:\t$opt_s\n" if defined $opt_s;
print "-d:\t$opt_d\n" if defined $opt_d;
print "-p:\t$opt_p\n" if defined $opt_p;
print "-m:\t$opt_m\n" if defined $opt_m;
print "-n:\t$opt_n\n" if defined $opt_n;
print "-r:\t$opt_r\n" if defined $opt_r;
print "-M:\t$opt_M\n" if defined $opt_M;
print "-t:\t$opt_t\n" if defined $opt_t;
print"\n=============================================================================\n\n";
die "$usage" if !defined $opt_s;
die "$usage" if !defined $opt_d;


#addition in command to count frequency of cis element
my$process_output='|grep -c ">" ';


#Assign default value to $opt_m if not provided by the user. It cannot be less than 3
$opt_m=1000 if !defined $opt_m;
die "\nError:\nNumber of datasets(-m) cannot be < 3\n" if($opt_m < 3);

$/="\n>";
open RSEQ,"$opt_s" or die "\n\nCannot find file containing sequences to test\n\n$usage";
# counting the number of sequences in real seq file. Assign this number to $opt_n if it is not defined in options
while(<RSEQ>){
      (my$header,my@seq)=split(/\n/,$_);
      $header=~s/>//;
      my$rseq_1line= join("",@seq);
      $rseq_1line=~ s/\s//g;
      $rseq{$header}=$rseq_1line;
      
}
close RSEQ;
$opt_n=scalar keys %seq if !defined $opt_n;


#print "Value of opt_n :$opt_n\n";

#open output file to collect values
open OUT,">Cis_elements_stats.txt";
#print OUT"Cis_Site\tTrue_Hits\tExpected_Hit\tVariance\tStDev\tTrimmed_mean\tZ-score\n";
printf OUT"%-15s %-10s %-10s %-10s %-10s %-10s %-10s\n","Cis_Site","True_Hits","Expec_Hit","Var","Std","Trimmed mean","Z-score";

#print "Cis_Element\tReal_Hit\tExpected_Hit\tVar\tStd Dev\tTrimmed mean\tZ-score\n";


# open sequence file, convert names to numbers.
open SEQ,"$opt_d" or die "Cannot find sequence file\n$usage";
my $count=1;
my @file_list=(); # initiate an array @file-list to capture the names of files generated. use this to delete it later
while(<SEQ>){
      (my$header,my@seq)=split(/\n/,$_);
      $header=$count;
      
      my$seq_1line= join("",@seq);
      $seq_1line=~ s/\s//g;
      $seq_1line=~ s/>//g;

     # $header=~ s/\s*//g;
      $header=~ s/[^\d]//g;

      $seq{$header}=$seq_1line;
      $count++;
}
 close SEQ;
#count the number of sequences in %seq hash use it as range for random number generator
my $range= scalar keys %seq;

#print "Range is : $range\n";
      
# create m datasets of n sequences
#print "\n******************Creating Surrogate datasets*******************\n";
for(my$i=1;$i<=$opt_m;$i++){
       open DATASET,">ds_$i.temp";
       push(@file_list,"ds_$i.temp");
       for(my$n=1;$n<=$opt_n;$n++){
                my $random=int(rand($range));
                print DATASET">$random\n$seq{$random}\n";
       }
	   close DATASET;
	
}
#print "\n******************Done Creating Surrogate datasets*******************\n";
printf "%-15s %-10s %-10s %-10s %-10s %-10s %-10s\n","Cis_Site","True_Hits","Expec_Hit","Var","Std","Trimmed mean","Z-score";
#print "Cis_Site\tTrue_Hits\tExpec_Hit\tVar\tStd\tTrimmed mean\tZ-score\n";


#open pattern file and run patmatch for each pattern. count the number of hits in m datasets individually and in real seq.
#command to run patmatch: perl scan_pipeline.pl [residue] [pattern] [sequence file] [number of mismatches] [mismatch types]
$/="\n";

print "Provide the pattern file containing cis elements" if !defined $opt_p;
$opt_p=<STDIN> if !defined $opt_p;

#print "\nprovide residue type (-n,-p,-c)" if !defined $opt_r;
$opt_r='n' if !defined $opt_r;
$opt_M=0 if !defined $opt_M;
$opt_t="" if !defined $opt_t;

open PATTERN,"$opt_p" or die "cannot open $opt_p \n$usage";
while(<PATTERN>){
	if(/^\s*$/){next;} 
	my$pattern=$_ if $_ ne "";
	$pattern=~s/[^a-zA-Z]//g;
	
	# CONVERT AMBIGOUS CODES TO PATTERN MATCH
	$pattern=~s/K/\[GT\]/gi;
	$pattern=~s/M/\[AC\]/gi;
	$pattern=~s/R/\[AG\]/gi;
	$pattern=~s/Y/\[CT\]/gi;
	$pattern=~s/S/\[CG\]/gi;
	$pattern=~s/W/\[AT\]/gi;
	$pattern=~s/U/T/gi;
	$pattern=~s/B/\[CGT\]/gi;
	$pattern=~s/V/\[ACG\]/gi;
	$pattern=~s/H/\[ACT\]/gi;
	$pattern=~s/D/\[AGT\]/gi;
	$pattern=~s/K/\[GT\]/gi;
	$pattern=~s/N/\[ACGT\]/gi;
	next if $pattern eq "";
	#count frequency in real seq
	#my @command=('perl scan_pipeline.pl',"-$opt_r",$pattern,$opt_s,$opt_M,$opt_t);
	my @command=('awk -F ',"$pattern", "\'{s+=(NF-1)} END {print s}\'"); 
	
	#pop @command;
	my $command1=join(" ",@command,$opt_s);
	#print"Command being sent for real hits:$command1\n";
	my $real_seq_hits=`$command1`;
	#print "real hit:$real_seq_hits\n";	
	my @hits=();
	foreach my $seqfile(@file_list){
		my @command=('awk -F ',"$pattern", "\'{s+=(NF-1)} END {print s}\'"); 
		#pop @command;
		my $command1=join(" ",@command,$seqfile);
		#print "Sending command;$command1\n";
		my $Total_hits=`$command1`;
		push(@hits,$Total_hits);
		$Total_hits=0;
	}
	#my$average = sum(@hits)/@hits;
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@hits); 

	my $mean = $stat->mean();
	my $var  = $stat->variance();
	my $tm   = $stat->trimmed_mean(.25);
	my $std  = $stat->standard_deviation();
	$Statistics::Descriptive::Tolerance = 1e-10;

                       
#    my $Zscore=($real_seq_hits - $mean)/$std if ($std!=0);
    
    my $Zscore=($real_seq_hits - $mean)/($std/sqrt()) if ($std!=0);

    
    #Change the pattern to ambigous codes
    $pattern=~s/\[GT\]/K/gi;
	$pattern=~s/\[AC\]/M/gi;
	$pattern=~s/\[AG\]/R/gi;
	$pattern=~s/\[CT\]/Y/gi;
	$pattern=~s/\[CG\]/S/gi;
	$pattern=~s/\[AT\]/W/gi;
	$pattern=~s/\[CGT\]/B/gi;
	$pattern=~s/\[ACG\]/V/gi;
	$pattern=~s/\[ACT\]/H/gi;
	$pattern=~s/\[AGT\]/D/gi;
	$pattern=~s/\[GT\]/K/gi;
	$pattern=~s/\[ACGT\]/N/gi;

	
	if(defined $Zscore){
	printf OUT"%-15s%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n",$pattern,$real_seq_hits,$mean,$var,$std,$tm,$Zscore;
	printf "%-15s%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n",$pattern,$real_seq_hits,$mean,$var,$std,$tm,$Zscore;
	}
	else{
	$Zscore=' NA';
	printf OUT"%-15s%10.2f%10.2f%10.2f%10.2f%10.2f%-12s\n",$pattern,$real_seq_hits,$mean,$var,$std,$tm,$Zscore;
	printf "%-15s%10.2f%10.2f%10.2f%10.2f%10.2f%10s\n",$pattern,$real_seq_hits,$mean,$var,$std,$tm,$Zscore;
	}
	
	
	#print OUT"$pattern\t$real_seq_hits\t$average\n";
}

system('rm *.temp');
