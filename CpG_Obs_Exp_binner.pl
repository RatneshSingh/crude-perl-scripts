#!/usr/bin/perl -w
use strict;
use Getopt::Std;

our($opt_t,$opt_o,$opt_s,$opt_d,$opt_m,$opt_l);
getopt('tosdml');


my$usage='
Usage: perl Script options....
-t	parsed newcpgreport in table format
	(not used when using multiple files -m )
-o	output file to save results
-s	smallest Obs/Exp value to start with [0]
-d	bin size for Obs/Exp values
-m	max bin size [ Auto setect]
-l	list of files to be parsed
';

if(defined$opt_o){open OUT,">$opt_o";} else {die "$usage";}
my$cycle=0;
# read multiple files to parse
if($opt_l){open FILES,"$opt_l";while(<FILES>){parse($_);}}
elsif($opt_t) {parse($opt_t);}
else{die "Provide table file or list of files\n\n$usage\n";}
$opt_d=0.05 if !defined $opt_d;




########################################################################################



sub parse{

	my$file=shift@_;
	chomp($file);
	$file=~s/\s//g;
	open TABLE,"$file";
	my@ObsExp;
	while (<TABLE>){
	
		# next if header line or an empty line found.
		if($_=~/Seq/||$_=~/start/||$_=~/end/ ||$_=~/^\s+$/){next;}
	
		# split line and catch values.
		my($seq,$start,$end,$size,$Sum_GC,$percent_GC,$Obs_Exp)=split(/\s+/,$_);
		chomp($seq,$start,$end,$size,$Sum_GC,$percent_GC,$Obs_Exp);
	
		#create arrays of all the things you need.
		$Obs_Exp=~s/\s+//g;
		push(@ObsExp,$Obs_Exp) if ($Obs_Exp ne "");
	}

	# set max obs/exp number to next integer of the max Obs values possible

	my$max_Obs=0;
	if(!$opt_m){
		foreach (@ObsExp){next if $_=~/^\s*$/;if($_>$max_Obs){$max_Obs=(int($_/$opt_d)+1)*$opt_d;}}
	}
	else{$max_Obs=$opt_m;}


	# assign complete number for each obs value and count the times they are present. 
	# count the total number of cPG islands

	my$count=0;
	my%modObs_Value=();	
	foreach my$Obs_Value(@ObsExp){
		my$Obs_ValueNew=(int($Obs_Value/$opt_d)+1);
		$modObs_Value{$Obs_ValueNew}++;
		$count++;
	}
	if($count==0){$count=1;}
	my$max_bin=$max_Obs/$opt_d;

	# print headings(bin sizes)
	if($cycle<1){
		for (my$i=1;$i<=$max_bin;$i++){
			my$bins=($i-1)*$opt_d;
			my$binl=$i*$opt_d;
			print OUT"\t$bins\-$binl";
			$cycle++;
		}
	}
	
	
	
	
	
	print OUT"\n$file testline";
	# printing % obs values for each bin.
	for (my$i=1;$i<=$max_bin;$i++){
		$modObs_Value{$i}=0 if !defined $modObs_Value{$i};
		my$percent_obs=$modObs_Value{$i}*100/$count;
		printf OUT"\t%.2f",$percent_obs;
	}
#	print OUT"\n";


}




