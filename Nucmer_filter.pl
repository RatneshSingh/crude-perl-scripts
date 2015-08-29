#!/usr/bin/perl -w
use strict;
use Getopt::Std;

our($opt_f,$opt_m,$opt_x,$opt_l,$opt_q,$opt_r,$opt_b,$opt_o,$opt_t,$opt_s,$opt_p,$opt_h);
our(@selected,%uniq_pair);
getopt('fmxlqrbotsph');

my$usage='
Please run show-coords utility with options; -c -l -q   on nucmer output (.delta) before feeding to this script.
usage: perl script options....

-f	Name of show-coords output file.
-m	Minimum percent to start binning[1]
-x	maximum percent to stop binning at[100]
-l	Alignment length cutoff for filtered sequences[1] 
-q	Query coverage cutoff for filtering[1]
-r	Reference coverage cutoff for filtering[1]
-b	bin size[1]
-o	outputfile
-t	yes/no. Print filtered table [no]
-s	yes/no/solo. Include self hits?
	yes: Count hits to self.e.g.  Include alignmens of "A" to "A".[no]
	no : Do  not count all the hits to self.
	solo: Count the hits to self only for those showing non-self hits also.
-p	yes/no. Use only uniq pairs for calculations. Removes reverse pairs from counting.[no]
	e.g. Between alignments of "A" to "B" and "B" to "A" only one alignment
	( with the largest Aln length) is included in the calculations.
-h	string. Print help and exit.
';

if($opt_h){die $usage}


$opt_m=1 if !defined $opt_m;
$opt_x=100 if !defined $opt_x;
$opt_l=1 if !defined $opt_l;
$opt_q=1 if !defined $opt_q;
$opt_r=1 if !defined $opt_r;
$opt_b=1 if !defined $opt_b;
$opt_t='no' if !defined $opt_t;
$opt_s='no' if !defined $opt_s;
$opt_p='no' if !defined $opt_p;


# Acquire input to parse
if($opt_f){open INFILE, $opt_f}
elsif($ARGV[0]){open INFILE,$ARGV[0]}
else{die "Provide a show-coords output file to parse\n\n$usage\n";}


if(!$opt_o){open OUT,">Alignment_binned_.$opt_f" ;}
else {open OUT,">$opt_o" ;}


my(%percent_freq);
for(my$i=$opt_m+$opt_b;$i<=$opt_x;$i=$i+$opt_b){
	
	print OUT "\t",sprintf("%0.2f",$i);
	print "\t",sprintf("%0.2f",$i);

}


print OUT"\n";
print "\n";
# initialize bins and set value to zero;
for(my$i=$opt_m;$i<=$opt_x;$i=$i+$opt_b){$percent_freq{$i}=0;}
$.=0;
#skip initial 5 lines and then parse
do <INFILE> while $.<=5;

while(<INFILE>){
	
			# next if header line or an empty line found.
			if($_=~/^\s+$/){next;}
	
			# split line and catch values.
			# header is : [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS
			my($S1_E1,$S2_E2,$LEN1_LEN2,$PER_IDY,$LENR_LENQ,$COVR_COVQ,$NAMER_NAMEQ)=split(/\|/,$_);
			chomp($S1_E1,$S2_E2,$LEN1_LEN2,$PER_IDY,$LENR_LENQ,$COVR_COVQ,$NAMER_NAMEQ);
			
			#Remove trailing spaces at the start of each variable to assist in proper splitting usng space as delimiter
			foreach($S1_E1,$S2_E2,$LEN1_LEN2,$PER_IDY,$LENR_LENQ,$COVR_COVQ,$NAMER_NAMEQ){s/^\s+//g}
			
			
			my($S1,$E1)=split(/\s+/,$S1_E1);
			my($S2,$E2)=split(/\s+/,$S2_E2);
			my($LEN1,$LEN2)=split(/\s+/,$LEN1_LEN2);
			#my($PER,$IDY)=split(/\s+/,$PER_IDY);
			my($LENR,$LENQ)=split(/\s+/,$LENR_LENQ);
			my($COVR,$COVQ)=split(/\s+/,$COVR_COVQ);
			my($NAMER,$NAMEQ)=split(/\s+/,$NAMER_NAMEQ);
			
			# Remove all the white spaces from name and values.
			foreach($S1,$E1,$S2,$E2,$LEN1,$LEN2,$PER_IDY,$LENR,$LENQ,$COVR,$COVQ,$NAMER,$NAMEQ){s/\s+//g}
			
			#print "$S1\t$E1\t$S2\t$E2\t$LEN1\t$LEN2\t$PER_IDY\t$LENR\t$LENQ\t$COVR\t$COVQ\t$NAMER\t$NAMEQ\n";
			
			# if $opt_S then filter self hits.
			if($opt_s eq 'no' && $opt_p eq 'no'){
				if(($LEN1>=$opt_l || $LEN2>=$opt_l) && $COVR>=$opt_r && $COVQ>=$opt_q && $PER_IDY>=$opt_m && $PER_IDY<=$opt_x && $NAMER ne $NAMEQ){
				
				#print "$S1\t$E1\t$S2\t$E2\t$LEN1\t$LEN2\t$PER_IDY\t$LENR\t$LENQ\t$COVR\t$COVQ\t$NAMER\t$NAMEQ\n" if defined $opt_t;
				push(@selected,"$S1\t$E1\t$S2\t$E2\t$LEN1\t$LEN2\t$PER_IDY\t$LENR\t$LENQ\t$COVR\t$COVQ\t$NAMER\t$NAMEQ") if $opt_t ne 'no';
				for(my$i=$opt_m;$i<=$opt_x;$i=$i+$opt_b){
					
						if ($PER_IDY>$i && $PER_IDY <=($i+$opt_b)){$percent_freq{$i+$opt_b}++;}
					}	
				}
			
			}
			
			elsif($opt_p ne 'no' && $opt_s eq 'no'){
				
				# if $opt_S then filter self hits.
				if(($LEN1>=$opt_l || $LEN2>=$opt_l) && $COVR>=$opt_r && $COVQ>=$opt_q && $PER_IDY>=$opt_m && $PER_IDY<=$opt_x && $NAMER ne $NAMEQ){
				
					#print "$S1\t$E1\t$S2\t$E2\t$LEN1\t$LEN2\t$PER_IDY\t$LENR\t$LENQ\t$COVR\t$COVQ\t$NAMER\t$NAMEQ\n" if defined $opt_t;
							
					push(@selected,"$S1\t$E1\t$S2\t$E2\t$LEN1\t$LEN2\t$PER_IDY\t$LENR\t$LENQ\t$COVR\t$COVQ\t$NAMER\t$NAMEQ") if $opt_t ne 'no';
							
					# Arange the namer and nameq in alphabetical order.
					my$uniq_name=$NAMER lt $NAMEQ ? ($NAMER.$NAMEQ):($NAMEQ.$NAMER);
							
					if (!exists $uniq_pair{'aln_length'}{$uniq_name}){
						$uniq_pair{'aln_length'}{$uniq_name}=($LEN1+$LEN2)/2 ;
						$uniq_pair{'per_idy'}{$uniq_name}=$PER_IDY
					}
					elsif(exists $uniq_pair{'aln_length'}{$uniq_name}){
								
						# Select the Alignment length and Per_IDY fo rthe longest alignment if more than one alignment for same pair are detected.
						$uniq_pair{'aln_length'}{$uniq_name}=($LEN1+$LEN2)/2 if $uniq_pair{'aln_length'}{$uniq_name} < ($LEN1+$LEN2)/2;
						$uniq_pair{'per_idy'}{$uniq_name}=$PER_IDY if $uniq_pair{'aln_length'}{$uniq_name} < ($LEN1+$LEN2)/2 
					}
				}	
			}	
			elsif($opt_p ne 'no' && $opt_s ne 'no'){
				if(($LEN1>=$opt_l || $LEN2>=$opt_l) && $COVR>=$opt_r && $COVQ>=$opt_q && $PER_IDY>=$opt_m && $PER_IDY<=$opt_x ){
		
				#print "$S1\t$E1\t$S2\t$E2\t$LEN1\t$LEN2\t$PER_IDY\t$LENR\t$LENQ\t$COVR\t$COVQ\t$NAMER\t$NAMEQ\n" if defined $opt_t;
				
				push(@selected,"$S1\t$E1\t$S2\t$E2\t$LEN1\t$LEN2\t$PER_IDY\t$LENR\t$LENQ\t$COVR\t$COVQ\t$NAMER\t$NAMEQ") if $opt_t ne 'no';
				
				# Arange the namer and nameq in alphabetical order.
				my$uniq_name=$NAMER lt $NAMEQ ? ($NAMER.$NAMEQ):($NAMEQ.$NAMER);
				
				if (!exists $uniq_pair{'aln_length'}{$uniq_name}){
					$uniq_pair{'aln_length'}{$uniq_name}=($LEN1+$LEN2)/2 ;
					$uniq_pair{'per_idy'}{$uniq_name}=$PER_IDY
				}
				elsif(exists $uniq_pair{'aln_length'}{$uniq_name}){
					
					# Select the Alignment length and Per_IDY fo rthe longest alignment if more than one alignment are detected.
					$uniq_pair{'aln_length'}{$uniq_name}=($LEN1+$LEN2)/2 if $uniq_pair{'aln_length'}{$uniq_name} < ($LEN1+$LEN2)/2;
					$uniq_pair{'per_idy'}{$uniq_name}=$PER_IDY if $uniq_pair{'aln_length'}{$uniq_name} < ($LEN1+$LEN2)/2 
					
				}
					}					
			}
			
					
			else{
				if(($LEN1>=$opt_l || $LEN2>=$opt_l) && $COVR>=$opt_r && $COVQ>=$opt_q && $PER_IDY>=$opt_m && $PER_IDY<=$opt_x){
				
				push(@selected,"$S1\t$E1\t$S2\t$E2\t$LEN1\t$LEN2\t$PER_IDY\t$LENR\t$LENQ\t$COVR\t$COVQ\t$NAMER\t$NAMEQ") if $opt_t ne 'no';
				#print "$S1\t$E1\t$S2\t$E2\t$LEN1\t$LEN2\t$PER_IDY\t$LENR\t$LENQ\t$COVR\t$COVQ\t$NAMER\t$NAMEQ\n" if defined $opt_t;
					for(my$i=$opt_m;$i<=$opt_x;$i=$i+$opt_b){
					
						if ($PER_IDY>$i && $PER_IDY <=($i+$opt_b)){$percent_freq{$i+$opt_b}++;}
					}	
				}
			
			}
			
			
			
			
	}




	if($opt_p ne 'no'){
		
		foreach(keys %{$uniq_pair{'aln_length'}}){
			for(my$i=$opt_m;$i<=$opt_x;$i=$i+$opt_b){
				if ($uniq_pair{'per_idy'}{$_}>$i && $uniq_pair{'per_idy'}{$_} <=($i+$opt_b)){$percent_freq{$i+$opt_b}++;}
			}		
				
		}		
		
	}








	if ($opt_t ne 'no'){foreach(@selected){print "$_\n" }}
	for(my$i=$opt_m;$i<=$opt_x;$i=$i+$opt_b){
		print OUT"\t$percent_freq{$i+$opt_b}" if defined $percent_freq{$i+$opt_b};
		print "\t$percent_freq{$i+$opt_b}" if defined $percent_freq{$i+$opt_b};
	}

print "\n";