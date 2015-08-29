#!/usr/bin/perl -w
use strict;
use Getopt::Std;


our($opt_f,$opt_o,$opt_d,$opt_l);
getopt('fodl');

my$usage="\n\nThis script is to parse newcpgreport output to the table\n\nusage:\nperl script -f newcpgreport-file -o output_file\n";

open FILE,"$opt_f" or die "$usage\n";
open OUT,">$opt_o" if defined $opt_o;

if(defined $opt_d){open OUT2,">$opt_d";} else{$opt_d=$opt_f.'.domaindraw';open OUT2,">$opt_d";}

if(!defined $opt_l){print "Provide the length of the sequence:"; $opt_l=<STDIN>; chomp $opt_l;}

$/="\/\/";
print OUT"Seq\tstart\tend\tsize\tSum_GC\t\%GC\tObs\/Exp" if defined $opt_o;
print "Seq\tstart\tend\tsize\tSum_GC\t\%GC\tObs\/Exp" if !defined $opt_o;


if($opt_o){
	while(<FILE>){
		my $ID=();
		my$length=0;
		my@lines=split(/\n/,$_);
			foreach my$new_line(@lines){
				if($new_line=~/ID\s+([\w\d\W\D]+)\s+(\d+)\s+BP/){$ID=$1;$length=$2;}
				elsif($new_line=~/CpG\s+island\s+(\d+)\.+(\d+)/){print OUT"\n$ID\t$1\t$2";
				my$newID=substr($ID,0,10);
				print OUT2"\n$newID\&CGI\&$1:$2\&S\&\#ff0000\&CGI\&\&Y\&$opt_l\&";}
				elsif($new_line=~/size=(\d+)/){print OUT"\t$1";}
				elsif($new_line=~/Sum C\+G=(\d+)/){print OUT"\t$1";}
				elsif($new_line=~/Percent CG=([.\d]+)/){print OUT"\t$1";}
				elsif($new_line=~/ObsExp=([.\d]+)/){print OUT"\t$1";}
				else{next;}
			}	

	}
}

else{
	while(<FILE>){
		my $ID=();
		my$length=0;
		my@lines=split(/\n/,$_);
			foreach my$new_line(@lines){
				if($new_line=~/ID\s+([\w\d\W\D]+)\s+(\d+)\s+BP/){$ID=$1;$length=$2;}
				
				elsif($new_line=~/CpG\s+island\s+(\d+)\.+(\d+)/){
				print "\n$ID\t$1\t$2";
				print OUT2"\n$ID\&CGI\&$1:$2\&S\&\#ff0000\&CGI\&\&Y\&$opt_l\&";
				}

				elsif($new_line=~/size=(\d+)/){print "\t$1";}
				elsif($new_line=~/Sum C\+G=(\d+)/){print "\t$1";}
				elsif($new_line=~/Percent CG=([.\d]+)/){print "\t$1";}
				elsif($new_line=~/ObsExp=([.\d]+)/){print "\t$1";}
				else{next;}
			}	

	}
}


#########################################################################################
# create inputfilr for domaindraw. This file can be used as input to draw domains at following
# website: http://domaindraw.imb.uq.edu.au./
# format required by the website is: X&CGI&2:25&S&#ff0000&CGI&&Y&100&
#########################################################################################


print "\n";