#!/usr/bin/perl    -w
use strict;
system(cls);
print "\n\n
|**************************************************************************|
|X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X |
|**************************************************************************|
|This script will read blast regular output file and parse through counting| 
|the occurrence of evalues provided by the user.you can give as many evalue|
|you want. Give all the values seperated by spaces. You can use it on the  |
|command line or interactively. for command line use as                    |
|=> program_name  blastfile_name  evalues_seperated by_spaces              |
|**************************************************************************|
|X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X |
|**************************************************************************|\n\n";
#*******************************************************************************

my ($e_limit,$e_count,$score_line,$e_value,$filename,@e_limit,@e_value,$real_nohits,$nohit_count);
#*******************************************************************************
#ask for blast file and open it. die if cannot
if ($ARGV[0]){$filename=$ARGV[0];} 

else { print "\nEnter the blast output file name::\n";
$filename=<STDIN>;
}  
chomp $filename;
open BLASTFILE,"$filename" or die "\nCant open $filename\n";
#*******************************************************************************

#ask for evalue(s). to count hits and non hits for.
if($ARGV[1]){
@e_limit=@ARGV;
shift @e_limit;
}
else{

print "\n\nGive Expect value(s). it should be number or in the format like 1e-10 seperated by spaces if many \n";

@e_limit=split (/\s/,<STDIN>);}
#*******************************************************************************

#open outout file with the prefix of blastfile name and write the header in the first line.
open OUT,">$filename"."_"."multiple_evalues.txt";
print OUT"elimit\tevalue_count\tnohit_found\n";
print "elimit\tevalue_count\tnohit_found\n";
#*******************************************************************************


#read blast file and make array of evalues from first hit for each sequences
$/="Query=";
$real_nohits=0;
while (<BLASTFILE>){ #1

			if(/Score =/){#2
			$score_line=$_;$score_line =~ /Score.*Expect.*= ([-.e\d]+)/;
			$e_value = $1;
				if($e_value=~ /^e/){$e_value="1".$e_value;}
			push (@e_value,$e_value);
						}#2
			if(/(No hits found)/){$real_nohits++;}

					}#1
#*******************************************************************************

#count the occurrence of evalues equal or less than e_limit provided.
foreach $e_limit(@e_limit){##1

chomp $e_limit;
$e_count=0;

		foreach $e_value(@e_value){##2
				if($e_value<=$e_limit){$e_count++;}else {$nohit_count++;}
								}##2
$nohit_count=$nohit_count+$real_nohits;
print OUT"$e_limit\t$e_count\t$nohit_count\n";
print "$e_limit\t$e_count\t$nohit_count\n";
}##1
#*******************************************************************************
close OUT;
close BLASTFILE;  
exit;
