#!/usr/bin/perl -w
use strict;
use Getopt::Std;

our($opt_t,$opt_o,$opt_w,$opt_l,$opt_m);
getopt('towlm');


my$usage='
Usage: perl Script options....

-l	List of files (parsed newcpgreport in table format) to be scanned.
-t	parsed newcpgreport in table format
	(not used when using multiple files -l )
-o	output file to save results
-w	Window size (bp) for scanning [1000]
-m	length of the Molecule
';

if(defined$opt_o){open OUT,">$opt_o";};# else {die "$usage";}
$opt_w = 100000 if !defined $opt_w;
die($usage) if !defined $opt_m;
# read multiple files to parse
if($opt_l){	open FILES,"$opt_l"; while(<FILES>){parse($_);} }

elsif($opt_t) {parse($opt_t);}

else{die "Provide table file or list of files\n\n$usage\n";}




########################################################################################



sub parse{

	my$file=shift@_;
	chomp($file);
	#$file=~s/\s//g;
	open TABLE,"$file";
	my(@StartEndPair,%window,@window,%CGIfeatures);

	while (<TABLE>){
	
		# next if header line or an empty line found.
		if($_=~/Seq/||$_=~/start/||$_=~/end/ ||$_=~/^\s+$/){next;}
	
		# split line and catch values.
		my($seq,$start,$end,$size,$Sum_GC,$percent_GC,$Obs_Exp)=split(/\s+/,$_);
		chomp($seq,$start,$end,$size,$Sum_GC,$percent_GC,$Obs_Exp);
	
		#create arrays of all the things you need.
		my$start_end_pair=$start.'/'.$end;
		push(@StartEndPair,$start_end_pair);
		$CGIfeatures{$start_end_pair}{'size'}=$size;
		
		
		
	}

	# Create window bin for scanning.
	
	for(my$i=1;$i<$opt_m;$i=$i+$opt_w){
		push(@window,$i);
		$window{$i}{'CGInumber'}=0;
		$window{$i}{'CGIlenght'}=0;
	}
	my$winlength=@window;
	# compare CpGIs in each window and count the number 

	foreach my$pair(@StartEndPair){
		(my$start,my$end)=split(/\//,$pair);
		$start=~s/\D//g;
		$end=~s/\D//g;	
		for(my$j=0;$j<$winlength-1;$j++){
			if($start>=$window[$j] && $end <= $window[$j+1]){
			$window{$window[$j]}{'CGInumber'}++;
			$window{$window[$j]}{'CGIlenght'}=$window{$window[$j]}{'CGIlenght'}+$CGIfeatures{$pair}{'size'};
			}
			else{next;}
		}
	}
	
	my@sortedwindows=keys %window;
	@sortedwindows=sort { $a <=> $b } @sortedwindows;
	
	if(defined $opt_o){	print OUT "\t".'CpGIs'."\t".'CpGIs_Length'."\n";}
		print  "\t".'CpGIs'."\t".'CpGIs_Length'."\n";

	
	foreach my$key(@sortedwindows){
		if(defined $opt_o){	
			print OUT $key."\t".$window{$key}{'CGInumber'}."\t".$window{$key}{'CGIlenght'}."\n";
		}
		print  $key."\t".$window{$key}{'CGInumber'}."\t".$window{$key}{'CGIlenght'}."\n";
	}	



}




