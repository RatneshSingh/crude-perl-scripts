#!/usr/bin/perl -w
use strict;
use Getopt::Std;

our($opt_b,$opt_m,$opt_l,$opt_o,$opt_f,$opt_e,$opt_n,$opt_p);

$opt_m=4;
$opt_b=0.2;
$opt_e=0.6;
$opt_n=10;
$opt_o='CpG_Size_bin.result';
$opt_p='newcpgreport';

getopt('bfmloenp');

my$usage='
usage: perl script options....

-l	list of files to parse
or
-f	file to parse
-e	Obs_Exp ratio minimum cutoff[0.6]
-m	maximum size of cluster to be considered in Kb[4]
-b	bin size in kb [0.2]
-o	outputfile
-n	number of char in file name to be used
-p	cpgcluster|cpgnewreport. Name of program name used to create table [newcpgreport.
';


# Acquire input to parse
my@files=();
if($opt_l){
	print "\nlist of file names has been provided:$opt_l\nReading file names";
	open LIST,"$opt_l" if defined $opt_l;
	foreach(<LIST>){push(@files,$_);}
}

elsif($opt_f){print "The file to parse:$opt_f\n";push(@files,$opt_f);}
else{die "Provide a file or list of files to parse\n\n$usage\n";}
print "\nfiles to parse: @files\n\n";


open OUT,">$opt_o" or die "Unable to creat output file to save results. Check your permissions";
my(%size_freq);
for(my$i=0;$i<=$opt_m;$i=$i+$opt_b){print OUT"\t$i";print "\t$i";}

foreach(@files){
	open TABLE,"$_" or die "Unable to find file $_\n\n";
	
	my$filename=substr($_,0,$opt_n);
	print "\n$filename";
	print OUT"\n$filename";
	# initialize and set value to zero;
	for(my$i=0;$i<=$opt_m;$i=$i+$opt_b){$size_freq{$i}=0;}

	while (<TABLE>){
			# next if header line or an empty line found.
			if($_=~/Seq/||$_=~/start/||$_=~/end/ ||$_=~/^\s+$/|| $_=~/CGI/||$_=~/From/||$_=~/OEratio/|| $_=~/pvalue/){next}
			my($seq,$start,$end,$size,$percent_GC,$Obs_Exp,$Num_CG,$Mean_Distance,$pvalue,$Sum_GC);
			# split line and catch values.
			($seq,$start,$end,$size,$Sum_GC,$percent_GC,$Obs_Exp)=split(/\s+/,$_) if $opt_p eq 'newcpgreport';
			($seq,$start,$end,$size,$percent_GC,$Obs_Exp,$Num_CG,$Mean_Distance,$pvalue)=split(/\s+/,$_) if $opt_p eq 'cpgcluster';
			chomp($seq,$start,$end,$size,$Sum_GC,$percent_GC,$Obs_Exp) if $opt_p eq 'newcpgreport';
			chomp ($seq,$start,$end,$size,$percent_GC,$Obs_Exp,$Num_CG,$Mean_Distance,$pvalue) if $opt_p eq 'cpgcluster';
			my$size_kb=$size/1000;
			$Obs_Exp=0.001 if (!$Obs_Exp);
			$Obs_Exp=0.001 if ($Obs_Exp eq "");

			for(my$i=0;$i<=$opt_m;$i=$i+$opt_b){
				if($Obs_Exp>=$opt_e){
					if (($size_kb>$i) && ($size_kb <=($i+$opt_b))){$size_freq{$i}++;}
				}
			}
			
	}

	for(my$i=0;$i<=$opt_m;$i=$i+$opt_b){
		print OUT"\t$size_freq{$i}";
		print "\t$size_freq{$i}";
				
	}
}
print "\n";