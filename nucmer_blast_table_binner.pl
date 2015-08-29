#!/usr/bin/perl -w
use strict;
use Getopt::Std;

our($opt_f,$opt_m,$opt_x,$opt_l,$opt_q,$opt_r,$opt_b,$opt_o,$opt_t,$opt_s,$opt_u,$opt_p,$opt_h);
our(@selected,%uniq_pair);
getopt('fmxlqrbotsuph');

my$usage='
Please run show-coords utility with options; -c -l -q   on nucmer output (.delta) before feeding to this script.
usage: perl script options....

-f	Name of inputfile. This file will be output of blast or show-coordsprogram.
-p	The name of program from which input is obtained. blast,show-coords[show-coords]
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
-u	yes/no. Use only uniq pairs for calculations. Removes reverse pairs from counting.[no]
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
$opt_u='no' if !defined $opt_u;
$opt_p='show-coords' if !defined $opt_p;
#$opt_p='blast' if !defined $opt_p;
if($opt_p eq 'show-coords' || $opt_p eq 'blast' || $opt_p eq 'blastn'){}else{die "Please provide proper value for -p\nAvailable options include: show-coords,blast, blastn\n"}
#$opt_f="test_blastFile.blastn";
# Acquire input to parse
if($opt_f){open INFILE, $opt_f}
elsif($ARGV[0]){open INFILE,$ARGV[0]}
else{die "Provide a input file (table format from blast or show-coords) to parse\n\n$usage\n";}


if(!$opt_o){open OUT,">Alignment_binned_.$opt_f" ;}
else {open OUT,">$opt_o" ;}

print OUT "$0   $opt_f    -m $opt_m   -x $opt_x    -l $opt_l    -q $opt_q    -r $opt_r    -b $opt_b    -s $opt_s    -u $opt_u    -p $opt_p\n\n";

# print the headers
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
if ($opt_p eq 'show-coords'){do <INFILE> while $.<=5}

while(<INFILE>){

			# next if empty line.
			if($_=~/^\s*$/){next}
			my($LENR,$LENQ,$NAMEQ,$NAMER,$PER_IDY,$LEN1,$LEN2,$MISMATCH,$GAP,$S1,$E1,$S2,$E2,$EVALUE,$BITSCORE,$COVR,$COVQ);
			if(lc$opt_p eq 'show-coords'){
				# split line and catch values.
				# header is : [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS
				my($S1_E1,$S2_E2,$LEN1_LEN2,$PER_IDY_temp,$LENR_LENQ,$COVR_COVQ,$NAMER_NAMEQ)=split(/\|/,$_);
				chomp($S1_E1,$S2_E2,$LEN1_LEN2,$PER_IDY_temp,$LENR_LENQ,$COVR_COVQ,$NAMER_NAMEQ);
				$PER_IDY=$PER_IDY_temp;
				#Remove leading spaces at the start of each variable to assist in proper splitting usng space as delimiter
				foreach($S1_E1,$S2_E2,$LEN1_LEN2,$PER_IDY,$LENR_LENQ,$COVR_COVQ,$NAMER_NAMEQ){s/^\s+//g}


				($S1,$E1)=split(/\s+/,$S1_E1);
				($S2,$E2)=split(/\s+/,$S2_E2);
				($LEN1,$LEN2)=split(/\s+/,$LEN1_LEN2);
				#($PER,$IDY)=split(/\s+/,$PER_IDY);
				($LENR,$LENQ)=split(/\s+/,$LENR_LENQ);
				($COVR,$COVQ)=split(/\s+/,$COVR_COVQ);
				($NAMER,$NAMEQ)=split(/\s+/,$NAMER_NAMEQ);

				# Remove all the white spaces from name and values.
				foreach($S1,$E1,$S2,$E2,$LEN1,$LEN2,$PER_IDY,$LENR,$LENQ,$COVR,$COVQ,$NAMER,$NAMEQ){s/\s+//g}

				#print "$S1\t$E1\t$S2\t$E2\t$LEN1\t$LEN2\t$PER_IDY\t$LENR\t$LENQ\t$COVR\t$COVQ\t$NAMER\t$NAMEQ\n";
			}
			elsif(lc$opt_p eq 'blastn' || lc$opt_p eq 'blast'){

				# split line and catch values.
				# header is : $queryId, $subjectId, $percIdentity, $alnLength, $mismatchCount, $gapOpenCount, $queryStart, $queryEnd, $subjectStart, $subjectEnd, $eVal, $bitScore, $Query_coverage
				($NAMEQ,$NAMER,$PER_IDY,$LEN1,$MISMATCH,$GAP,$S1,$E1,$S2,$E2,$EVALUE,$BITSCORE,$COVQ)=split(/\s+/,$_);
				chomp($NAMEQ,$NAMER,$PER_IDY,$LEN1,$MISMATCH,$GAP,$S1,$E1,$S2,$E2,$EVALUE,$BITSCORE,$COVQ);
				$LEN2=$LEN1;
				$COVR=sprintf("%.2f",$COVQ);
				$COVQ=$COVR;
				$LENR=$LENQ='NA';
				#Remove trailing spaces at the start of each variable to assist in proper splitting usng space as delimiter
				foreach($NAMEQ,$NAMER,$PER_IDY,$LEN1,$MISMATCH,$GAP,$S1,$E1,$S2,$E2,$EVALUE,$BITSCORE,$COVQ){s/\s+//g}

			}



			# next if filtering options are not met
			if(!(($LEN1>=$opt_l || $LEN2>=$opt_l) && $COVR>=$opt_r && $COVQ>=$opt_q && $PER_IDY>=$opt_m && $PER_IDY<=$opt_x)){next}

			# Dont proceed if self hit and -s no
			if($opt_s eq 'no' && $NAMER eq $NAMEQ){next}



			# if no uniq pairs and no self hits needed
			if(lc$opt_u eq 'no'){
				#if(($LEN1>=$opt_l || $LEN2>=$opt_l) && $COVR>=$opt_r && $COVQ>=$opt_q && $PER_IDY>=$opt_m && $PER_IDY<=$opt_x && $NAMER ne $NAMEQ){

					#print "$S1\t$E1\t$S2\t$E2\t$LEN1\t$LEN2\t$PER_IDY\t$LENR\t$LENQ\t$COVR\t$COVQ\t$NAMER\t$NAMEQ\n" if defined $opt_t;
					#push(@selected,"$S1\t$E1\t$S2\t$E2\t$LEN1\t$LEN2\t$PER_IDY\t$LENR\t$LENQ\t$COVR\t$COVQ\t$NAMER\t$NAMEQ") if $opt_t ne 'no';
					push(@selected,$_);# if lc$opt_t eq 'yes';
					for(my$i=$opt_m;$i<=$opt_x;$i=$i+$opt_b){

							if ($PER_IDY>$i && $PER_IDY <=($i+$opt_b)){$percent_freq{$i+$opt_b}++;}
					}


			}
			# if uniq pairs needed but no self hits
			elsif(lc$opt_u eq 'yes'){


					#print "$S1\t$E1\t$S2\t$E2\t$LEN1\t$LEN2\t$PER_IDY\t$LENR\t$LENQ\t$COVR\t$COVQ\t$NAMER\t$NAMEQ\n" if defined $opt_t;

					#push(@selected,"$S1\t$E1\t$S2\t$E2\t$LEN1\t$LEN2\t$PER_IDY\t$LENR\t$LENQ\t$COVR\t$COVQ\t$NAMER\t$NAMEQ") if $opt_t ne 'no';
					push(@selected,$_);# if lc$opt_t eq 'yes';
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
			#elsif(lc$opt_u eq 'yes' && lc$opt_s eq 'yes'){
			#
			#
			#	#print "$S1\t$E1\t$S2\t$E2\t$LEN1\t$LEN2\t$PER_IDY\t$LENR\t$LENQ\t$COVR\t$COVQ\t$NAMER\t$NAMEQ\n" if defined $opt_t;
			#
			#	#push(@selected,"$S1\t$E1\t$S2\t$E2\t$LEN1\t$LEN2\t$PER_IDY\t$LENR\t$LENQ\t$COVR\t$COVQ\t$NAMER\t$NAMEQ") if $opt_t ne 'no';
			#	push(@selected,$_) if lc$opt_t eq 'yes';
			#	# Arange the namer and nameq in alphabetical order.
			#	my$uniq_name=$NAMER lt $NAMEQ ? ($NAMER.$NAMEQ):($NAMEQ.$NAMER);
			#
			#	if (!exists $uniq_pair{'aln_length'}{$uniq_name}){
			#		$uniq_pair{'aln_length'}{$uniq_name}=($LEN1+$LEN2)/2 ;
			#		$uniq_pair{'per_idy'}{$uniq_name}=$PER_IDY
			#	}
			#	elsif(exists $uniq_pair{'aln_length'}{$uniq_name}){
			#
			#		# Select the Alignment length and Per_IDY fo rthe longest alignment if more than one alignment are detected.
			#		$uniq_pair{'aln_length'}{$uniq_name}=($LEN1+$LEN2)/2 if $uniq_pair{'aln_length'}{$uniq_name} < ($LEN1+$LEN2)/2;
			#		$uniq_pair{'per_idy'}{$uniq_name}=$PER_IDY if $uniq_pair{'aln_length'}{$uniq_name} < ($LEN1+$LEN2)/2
			#
			#	}
			#
			#}


			else{
				next;
				##push(@selected,"$S1\t$E1\t$S2\t$E2\t$LEN1\t$LEN2\t$PER_IDY\t$LENR\t$LENQ\t$COVR\t$COVQ\t$NAMER\t$NAMEQ") if $opt_t ne 'no';
				#push(@selected,$_) if lc$opt_t eq 'yes';
				##print "$S1\t$E1\t$S2\t$E2\t$LEN1\t$LEN2\t$PER_IDY\t$LENR\t$LENQ\t$COVR\t$COVQ\t$NAMER\t$NAMEQ\n" if defined $opt_t;
				#	for(my$i=$opt_m;$i<=$opt_x;$i=$i+$opt_b){
				#
				#		if ($PER_IDY>$i && $PER_IDY <=($i+$opt_b)){$percent_freq{$i+$opt_b}++;}
				#	}


			}




	}




	if($opt_u eq 'yes'){

		foreach(keys %{$uniq_pair{'aln_length'}}){
			for(my$i=$opt_m;$i<=$opt_x;$i=$i+$opt_b){
				if ($uniq_pair{'per_idy'}{$_}>$i && $uniq_pair{'per_idy'}{$_} <=($i+$opt_b)){$percent_freq{$i+$opt_b}++;}
			}

		}

	}







	print OUT"Frequency";

	for(my$i=$opt_m;$i<=$opt_x;$i=$i+$opt_b){
		print OUT"\t$percent_freq{$i+$opt_b}" if defined $percent_freq{$i+$opt_b};
		print "\t$percent_freq{$i+$opt_b}" if defined $percent_freq{$i+$opt_b};
	}

print
"\n\n\n*********************************************************************************************************************
		Filtered table
*********************************************************************************************************************\n" if lc$opt_t eq 'yes';
if (lc$opt_t eq 'yes'){foreach(@selected){print "$_\n" }}



print "\n\n Total number of hits passing the filter criteria:".scalar@selected."\n\n";