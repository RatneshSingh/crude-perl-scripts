#usr/bin/perl -w


#=====================================================================================
# This Script Reads Contig Assembly information from ContigExpress from VectorNTI
# suite. To get the info Double click contig and select Infopane window from 
# the bottom ribbon--> expand Fagment by clicking '+'--> go to 'Edit'--> 
#  Select 'Camera'-->Select Range 'All' and Copy to 'File' --> SClick 'OK'-->
# save as a Text file. This will be your infopane file. Get info files from all 
# the contigs in the assembly. Provide this file as Input file to the script. You 
# can provide multiple files as input file with names seperated by space.
# 
# Usage: perl Script_name InfopaneFile1 InfopaneFile2 .................
#
# Author: Ratnesh Singh
# Version: 1.0  
# Report any Bug to: ratnesh@hawaii.edu
#=====================================================================================



use strict;
use List::Util 'min';
use List::Util 'max';
my ($infopanefile,%read,@readlist);

if (defined $ARGV[0]){$infopanefile=$ARGV[0];} else {print "\nEnter file name containing info pane data from ContigExpress of Vector NTI\n";$infopanefile=<STDIN>;}




open OUT,">$infopanefile._MatePairSummary.txt";
#-------------Read infopane results files and assign start and orientation of reads on assembly to each read.----------------------------
foreach $infopanefile(@ARGV){
	
	
	chomp $infopanefile;

	if (-e $infopanefile){open INFOPANE,"<$infopanefile";} else {print "\n$infopanefile does not exist. Please check the file name\n"; next; }
	
	#--------Assign contig name from the first line of the file. Should not be empty---------------------------------------
	
	my $contig=<INFOPANE>;
	$contig=~s/\([\s\d\w]+?\)//;
	
	while(<INFOPANE>){

		my ($reads,$cordinate)= split(/:/);
	
		my($R_start,$R_end)=split(/[»«]+/,$cordinate);
	
		$reads=~s/\s*$//;
		$reads=~s/^\s*//;
	
		#-------remove non digit characters from co-ordinates.--------------------
		if($R_start=~/([\d]+)/){$R_start=$1;}
		if($R_end=~/([\d]+)/){$R_end=$1;}
	
		#-------Assign co-ordinates to respective reads in hash---------------------
		chomp ($reads,$R_start,$R_end,$cordinate,$contig);
		#$read{$reads}{'contig'}=$contig;
		#$read{$reads}{'Start'}=$R_start;
		#$read{$reads}{'Send'}=$R_end;
		
		if($cordinate=~/»/){$read{$reads}{'orientation'}='Plus';} 
		elsif($cordinate=~/«/){$read{$reads}{'orientation'}='Minus';}
		
		$read{$reads}{	'contig'=>$contig,
						'Start'=> $R_start,
						'Send'=>$R_end
					};

		
		$reads=~s/\.r1$|\.f1$//;
	
		#-------make list of all the reads in file--------------------------------
		push(@readlist,$reads);
	
	}
}



# collect uniques from reads list to array @reads2 to remove duplicates
my %hsh;
undef @hsh{@readlist};
my @reads2 = keys %hsh;


#find mate pairs

print OUT"Query\tF_on_contig\tR_on_contig\tF_start\tF_End\tR_start\tR_End\tF_orientation\tR_orientation\tDistance\torientationError\n";


foreach my $read2(@reads2){

chomp($read2);
my $F_read2=$read2.'.f1';
my $R_read2=$read2.'.r1';
my $distance=();
if($read{$F_read2}{'contig'} eq $read{$R_read2}{'contig'} ){


 my$distance_max = max($read{$F_read2}{'Start'},$read{$F_read2}{'Send'},$read{$R_read2}{'Start'},$read{$R_read2}{'Send'});
 my$distance_min = min($read{$F_read2}{'Start'},$read{$F_read2}{'Send'},$read{$R_read2}{'Start'},$read{$R_read2}{'Send'});
 
 $distance=$distance_max-$distance_min;
 
 }


elsif (!defined $read{$F_read2}{'contig'}) {$distance= 'F_read missing';}
elsif (!defined $read{$R_read2}{'contig'}) {$distance= 'R_read missing';}
else{$distance= 'Not on same contig';}
my $orien_error=();
if($read{$F_read2}{'orientation'} eq $read{$R_read2}{'orientation'}){$orien_error='Same Orientation';}
else{$orien_error='OK';}


print OUT"$read2\t$read{$F_read2}{'contig'}\t$read{$R_read2}{'contig'}\t$read{$F_read2}{'Start'}\t$read{$F_read2}{'Send'}\t$read{$R_read2}{'Start'}\t$read{$R_read2}{'Send'}\t$read{$F_read2}{'orientation'}\t$read{$R_read2}{'orientation'}\t$distance\t$orien_error\n";
print "$read2\t$read{$F_read2}{'contig'}\t$read{$R_read2}{'contig'}\t$read{$F_read2}{'Start'}\t$read{$F_read2}{'Send'}\t$read{$R_read2}{'Start'}\t$read{$R_read2}{'Send'}\t$read{$F_read2}{'orientation'}\t$read{$R_read2}{'orientation'}\t$distance\t$orien_error\n";


}


