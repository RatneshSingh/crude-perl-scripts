#!/usr/bin/perl    -w
use strict;
system(cls);
print "\n\n
|**************************************************************************|
|X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X |
|**************************************************************************|
|This script will read blast regular output file and parse through counting| 
|the occurrence of evalues below the threshold level.It starts with the    |
|lowest evalue assigned and reaches to highest limit by increment of user  |
|assigned values and calculates hits and nohits for each evalue between    |
|lower and higher limits.                                                  |
|**************************************************************************|
|X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X |
|**************************************************************************|";
my ($e_count,$score_line,$e_value,$filename,$elimit_L,$elimit_H,$inc,@e_value,$nohit_count,$real_nohits);
#**************************************************************************

#ask for blastfile name and open it. die if cannot.
if ($ARGV[0]){$filename=$ARGV[0];} 
else { print "\nEnter the blast output file name::\n";
$filename=<STDIN>;
}  
chomp $filename;
open BLASTFILE,"$filename" or die "\nCant open $filename\n";
#****************************************************************************

#ask for lower limit of e value to start counting with. 
if($ARGV[1]){$elimit_L=$ARGV[1];			}
else{
print "\n\nLOWER limt of Expect value (decimal numbers or in the format 1.0e-100): \n";
$elimit_L=<STDIN>;
	}
chomp $elimit_L;
#******************************************************************************

#ask for the upper limit till e value counting is done
if($ARGV[2]){$elimit_H=$ARGV[2];}
else{
print "\n\nUPPER limt of Expect value (decimal numbers or in the format 1.0e-100): \n";
$elimit_H=<STDIN>;
	}
chomp $elimit_H;
#*********************************************************************************

#ask for the increment value to increase in multiples till e value reaches higher limits
if($ARGV[3]){$inc=$ARGV[3];	}
else{
print "\n\nWrite the range of increment you want between lower and upper e values \n";
$inc=<STDIN>;
}
chomp $inc;
#*******************************************************************************


#open outout file with the prefix of blastfile name and write the header in the first line.
open OUT,">$filename"."_"."evalue_range.txt";
print OUT"Expect_Value\tHits_found\tNohits_Found\n";
print "\nExpect_value\tHits_found\tNo_Hit_found\n";
#*******************************************************************************

#read blast file and make array of evalues from first hit for each sequences
#$/="Query=";
my $chunk=0;
my $real_nohits=0;
while (<BLASTFILE>){ #1

			if(/Score =/){#2
			$score_line=$_;$score_line =~ /Score.*Expect.*= ([-.e\d]+)/;
			$e_value = $1;
				if($e_value=~ /^e/){$e_value="1".$e_value;}
			push (@e_value,$e_value);
						}#2
			if(/(No hits found)/){$real_nohits++;}
#$chunk++;
#if($chunk==11282){print "$_";}

					}#1
#*******************************************************************************


#count the occurrence of evalues equal or less than e_limit provided and write into file and stdout.
for(my $elimit=$elimit_H;$elimit>=$elimit_L;$elimit=$elimit/$inc){#open for loop

chomp $elimit;
$e_count=0;
$nohit_count=0;

		foreach $e_value(@e_value){ if($e_value<=$elimit){$e_count++;}else {$nohit_count++;}}

$nohit_count=$nohit_count+$real_nohits;
print "elimit= $elimit \t ecount= $e_count \t nohit=$nohit_count \t real nohits=$real_nohits \n"; 
print OUT"$elimit\t$e_count\t$nohit_count\n";
#print "$elimit\t$e_count\t$nohit_count\n";
if($elimit==0){exit;}
}#close for loop.
#**************************************************************************************

#close all open files and exit.
close OUT;
close BLASTFILE;  
exit,

