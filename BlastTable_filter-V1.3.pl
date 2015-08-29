#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Std;

#####################################################################################################
# This script is made to filter blast tabular output based on criteria for different blast values   #
#                                                                                                   #
# Author : Ratnesh Singh                                                                            #
# version 1.2                                                                                       #
# contact for bugs: ratnesh@hawaii.edu                                                              #
# last updated: 03/02/2011
#####################################################################################################


getopt('boqs');
our ($opt_b,$opt_o,$opt_q,$opt_s);
my(%qlength,%slength);

# like the shell getopt, "d:" means d takes an argument
print "-blast file/sequence: $opt_b\n" if defined $opt_b;
print "-print output as: $opt_o\n" if defined $opt_o;
print "Unprocessed by Getopt::Std:\n" if $ARGV[0];
foreach (@ARGV) {
  print "-- $_\n";
}

my $help="\n\nThis script will read the blast file in table format and 
and filter results based on provided criteria. 

usage:perl script -b blast file -o outputfile

-q table with query lengths.
-s table with subject lengths.

\n";


die "\nThere is no blast file specified with -b \n $help" if !defined ($opt_b);

# Collect options automatically from commandline for filtering blast results

# still to write............










# Collect options manually from user for filtering blast results

print"Select criteria for filtering. Select in order of use\n

0	Query	
1	Subject	
2	% Identity	
3	Aln length	
4	Mismatch	
5	Gap	
6	Q start	
7	Q end	
8	S start	
9	S end	
10	E value	
11	Bit score
12	Query coverage(in percent)
13	Subject coverage (in percent)
\n";

print "how many criterias you want to use for filtering (please type number):";
 my $numbers_of_Criteria=<STDIN>;
 chomp($numbers_of_Criteria);
 my %criteria;
 my $conditions;

 my@order=qw(Blank First Second Third Fourth Fifth Sixth Seventh Eighth Ninth Tenth Eleventh Twelfth Thirteenth);
 my@blast_options=qw(Query Subject Identity Aln_length Mismatch Gap Q_start Q_end S_start S_end E_value Bit_score Query_coverage_in_percent Subject_Coverage_in_percent);


 	my$Qcoverage='FALSE';
 	my$Scoverage='FALSE';

 for(my$i=1;$i<=$numbers_of_Criteria;$i++){
 	print "\nPlease print the number for criteria you want to use as your $order[$i] ($i) criteria from above list:";
 	my$crit=<STDIN>;
 	chomp ($crit);
 	
 	$crit=~s/\D//g;
 	$criteria{$i}{'criteria'}=$crit;
 	$Qcoverage='TRUE' if $crit==12;
 	$Scoverage='TRUE' if $crit==13;
 	
 	
 	# Collect criteria to filter. convert to digits. 
 	#if($crit > 1){
 	#	$crit=$crit*1;
 	#	$criteria{$i}{'criteria'}=$crit;
	#}
	
	#else{
	#	$criteria{$i}{'criteria'}=$crit;
	#}
	
	
	# collect cutoff value for above value to filter
 	print "\nPlease print the cutoff value of $blast_options[$crit] :";
	my$crit2=<STDIN>;
	chomp($crit2);
	if($crit > 1){
		$crit2=$crit2*1;
		$criteria{$i}{'cutoff'}=$crit2;
	}
	else{
	$criteria{$i}{'cutoff'}=$crit2;

	}

 }


$opt_o = 'filtered_'.$opt_b if !defined ($opt_o);

print "saving output in:$opt_o\n";
open OUT,">$opt_o" or die "cannot open Output file\n";


#---------------------------------------------------------
# collect sequence lengths from table for queries and subjects.

if(defined $opt_q){

	open QLENGTH,$opt_q;
	while(<QLENGTH>){
		(my$qheader,my$qlength)=split(/\s+/,$_);
		chomp($qheader);
		chomp($qlength);
		$qlength{$qheader}=$qlength;
	}
close QLENGTH;
}

if(defined $opt_s){

	open SLENGTH,$opt_s;
	while(<SLENGTH>){
		(my$sheader,my$slength)=split(/\s+/,$_);
		chomp($sheader);
		chomp($slength);

		$slength{$sheader}=$slength;
	}
close SLENGTH;

}




###################################################################
# parse information from blast file                               #
###################################################################

open BLAST,"$opt_b" or die "cannot read blast file \n";
my$count_out_put=0;	
while(<BLAST>){
			my $line=$_;
			chomp($line);
			if ($line=~/^\s+$/){ next ;};
			if ($line=~/query/i or /match/i or /score/ or /gap/ or /mismatch/){ next; } ;
			my @line_info= split(/\t/,$line);
			chomp(@line_info);
			$conditions=0;
			
			
			#calculate query coverage if user selected criteria 12 or 13
			if($Qcoverage eq 'TRUE'){
				$line_info[12]=sprintf "%.2f",(abs($line_info[6]-$line_info[7])+1)*100/$qlength{$line_info[0]}; 
				#print "$line_info[0] $line_info[6] - $line_info[7] * 100 \\ $qlength{$line_info[0]} \t $line_info[12]\n";
				if($line_info[12] >100){print "coverage values are more than 100. \nPlease check if you provided correct query and subject length tables with correct -s and -q flags\n";}
			}
 			if($Scoverage eq 'TRUE'){
 				$line_info[13]=sprintf "%.2f",(abs($line_info[8]-$line_info[9])+1)*100/$slength{$line_info[1]};
 				if($line_info[13] >100){print "coverage values are more than 100. \nPlease check if you provided correct query and subject length tables with correct -s and -q flags\n";}
 				#print "$line_info[1] $line_info[8]- $line_info[9] *100 \\ $slength{$line_info[1]}\n";
 			}

			if($Qcoverage eq 'TRUE' && $Scoverage eq 'TRUE'){$line=$line."\t".$line_info[12]."\t".$line_info[13];}
			elsif($Qcoverage eq 'TRUE'){$line=$line."\t".$line_info[12]."\t".'Not_Calculated';}
			elsif($Scoverage eq 'TRUE'){$line=$line."\t".'Not_Calculated'."\t".$line_info[13];}
			
			
			for(my$i=1;$i<=$numbers_of_Criteria;$i++){
				
				chomp($line_info[$criteria{$i}{'criteria'}]);
				chomp ($criteria{$i}{'cutoff'});
				
				#treat as number if not given Query or subject as criteria
				if($criteria{$i}{'criteria'} > 1){
					$line_info[$criteria{$i}{'criteria'}]=~s/[^\d\.-Ee]//g;
					$criteria{$i}{'cutoff'}=~s/[^\d\.-Ee]//g;
				}
				#Evaluate numbers if Query  or subject is not given as criteria for filtering
				if($criteria{$i}{'criteria'} > 1){

					if($line_info[$criteria{$i}{'criteria'}] >= $criteria{$i}{'cutoff'}){
						$conditions++;
					}
				}
				
				else{
					if($line_info[$criteria{$i}{'criteria'}] =~ m/$criteria{$i}{'cutoff'}/i){
						$conditions++;
					}
				}
			}
			
			if($conditions == $numbers_of_Criteria){
				print OUT"$line\n";
				$count_out_put++;	
			}
			else{ next;}
			
}			

print "$count_out_put lines passed the criterians and are written in output file\n";

close(OUT);
close(BLAST);
exit;


