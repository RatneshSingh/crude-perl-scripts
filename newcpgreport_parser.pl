#!/usr/bin/perl -w
use strict;
use Getopt::Std;


our($opt_f,$opt_o);
getopt('fo');

my$usage="\n\nThis script is to parse newcpgreport output to the table\n\nusage:\nperl script -f newcpgreport-file -o output_file\n";

open FILE,"$opt_f" or die "$usage\n";
open OUT,">$opt_o" if defined $opt_o;

$/="\/\/";
print OUT"Seq\tstart\tend\tsize\tSum_GC\t\%GC\tObs\/Exp" if defined $opt_o;
print "Seq\tstart\tend\tsize\tSum_GC\t\%GC\tObs\/Exp" if !defined $opt_o;


if($opt_o){
	my$totalCGI=0;
	my$numberCGI=0;
	while(<FILE>){
		my $ID=();
		my$length=0;
		my@lines=split(/\n/,$_);
			foreach my$new_line(@lines){
				if($new_line=~/ID\s+([\w\d\W\D]+)\s+(\d+)\s+BP/){$ID=$1;$length=$2;}
				elsif($new_line=~/CpG\s+island\s+(\d+)\.+(\d+)/){print OUT"\n$ID\t$1\t$2";}
				elsif($new_line=~/size=(\d+)/){print OUT"\t$1";$totalCGI=$totalCGI+$1;$numberCGI++;}
				elsif($new_line=~/Sum C\+G=(\d+)/){print OUT"\t$1";}
				elsif($new_line=~/Percent CG=([.\d]+)/){print OUT"\t$1";}
				elsif($new_line=~/ObsExp=([.\d]+)/){print OUT"\t$1";}
				else{next;}
			}	

	}
	print "\nTotal number of CGIs:$numberCGI\n";
	print "Total CGI length:$totalCGI\n";

}

else{
	my$totalCGI=0;
	my$numberCGI=0;
	while(<FILE>){
		my $ID=();
		my$length=0;
		my@lines=split(/\n/,$_);
			foreach my$new_line(@lines){
				if($new_line=~/ID\s+([\w\d\W\D]+)\s+(\d+)\s+BP/){$ID=$1;$length=$2;}
				elsif($new_line=~/CpG\s+island\s+(\d+)\.+(\d+)/){print "\n$ID\t$1\t$2";}
				elsif($new_line=~/size=(\d+)/){print "\t$1"; $totalCGI=$totalCGI+$1;$numberCGI++;}
				elsif($new_line=~/Sum C\+G=(\d+)/){print "\t$1";}
				elsif($new_line=~/Percent CG=([.\d]+)/){print "\t$1";}
				elsif($new_line=~/ObsExp=([.\d]+)/){print "\t$1";}
				else{next;}
			}	

	}
	print "\nTotal number of CGIs:$numberCGI\n";
	print "Total CGI length:$totalCGI\n";
	
}

print "\n";