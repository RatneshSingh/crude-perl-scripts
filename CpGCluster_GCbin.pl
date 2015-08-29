#!/usr/bin/perl -w
use strict;
use Getopt::Std;

our ($opt_m,$opt_d,$opt_t,$opt_o,$opt_p);
getopt('tmdop');


my$usage='
This script will read the output from CpGcluster program and 
bin GC percentage and OE ratios.

usage: perl script options...
-t	table from CpG cluster program
-m	Max bin 
-d	Bin size
-o	Output results
-p	Number of GC binning. User will be asked to input ranges

';

my(@GC_user,%GC);
# check parameters and assign default value if not provided
if(!$opt_t){die "\n\n\nError:Cannot find table from CpG cluster program\n\n$usage\n";}
else{open TABLE,"$opt_t";}


if($opt_p){
print "Please enter the max limit of $opt_p ranges you want to use. Larger first\n";
for(my$i=1;$i<=$opt_p;$i++){print"\n$i:"; my$range=<STDIN>; push(@GC_user,$range);}
}
else{@GC_user=(50,35,20);}

# sort array in reverse numerical order
@GC_user= sort{ $b <=> $a } @GC_user;

my$lengthGCuser=scalar@GC_user;


$opt_d=0.05 if !defined $opt_d;
$opt_m=2 if !defined $opt_m;
$opt_o="Gcbin".$opt_t if !defined $opt_o;


# Ask user to input max range of GC bin. collect in a@GC_user

#set number of bins to use for.
#my$binNum=int($opt_m/$opt_d)+1)
my$count=0;

if(!$opt_t){die "Cannot find table from CpG cluster program\n\n$usage\n";}
else{open TABLE,"$opt_t";}


LOOP: while(<TABLE>){
	if($_=~/CGI/){next LOOP;}
	if($_=~/^\s+$/){next LOOP;}
	$count++;
	my($CGI,$From,$To,$Length,$GCpcent,$OEratio,$CpG,$MeanDist,$pvalue)=split(/\s+/,$_);
	chomp($CGI,$From,$To,$Length,$GCpcent,$OEratio,$CpG,$MeanDist,$pvalue);
	
	# count CGIs for each GC% and in bin size;
	for(my$i=0;$i<=$lengthGCuser-1;$i++){

		if(int$GCpcent>=int$GC_user[$i]){
			for(my$j=0;$j<=$opt_m;$j=$j+$opt_d){
				if(($OEratio>=$j)&&($OEratio<($j+$opt_d))){$GC{$GC_user[$i]}{$j}++;}	
			}
		}
	}
	$CGI=$From=$To=$Length=$GCpcent=$OEratio=$CpG=$MeanDist=$pvalue=()
}


open OUT,">$opt_o";
for(my$j=0;$j<=$opt_m;$j=$j+$opt_d){printf OUT"\t%.2f - %.2f",$j,$j+$opt_d;}

#foreach my$GCpcent(keys%usedGCuser){
for(my$i=0;$i<=$lengthGCuser-1;$i++){
	printf OUT"\n$GC_user[$i]";
	for(my$j=0;$j<=$opt_m;$j=$j+$opt_d){
		if($GC{$GC_user[$i]}{$j}){printf OUT"\t%.2f",$GC{$GC_user[$i]}{$j};	}
		else{printf OUT"\t0.00";}
	}
			
}
print "Total number of CGIs found: $count";

