#!/usr/perl
use strict;
use warnings;
die "\n No DE result file provided on commandline." if !$ARGV[0];


open my $RES,"<", $ARGV[0];
my$up_P=my$up_Q=0;
my$down_P=my$down_Q=0;
my$tot_Q=my$tot_P=0;

	while(<$RES>){
			next if /^\s*$/;
	    my@elm=split /\s+/;
		  if($#elm < 9){  ### for older version of results
				if($elm[4] >= 1 || $elm[4] <= -1){$tot_P++ if $elm[7] <=0.05;$tot_Q++ if $elm[8] <= 0.05}
				if($elm[4] >= 1 ){$up_P++ if $elm[7] <=0.05;$up_Q++ if $elm[8] <= 0.05;}
				if($elm[4] <= -1){$down_P++ if $elm[7] <=0.05;$down_Q++ if $elm[8] <= 0.05;}
				}
		  else{  ### newer version(2.24) has two extra columns
		  	if($elm[6] >= 1 || $elm[6] <= -1){$tot_P++ if $elm[9] <=0.05;$tot_Q++ if $elm[10] <= 0.05}
				if($elm[6] >= 1 ){$up_P++ if $elm[9] <=0.05;$up_Q++ if $elm[10] <= 0.05;}
				if($elm[6] <= -1){$down_P++ if $elm[9] <=0.05;$down_Q++ if $elm[10] <= 0.05;}
			}
	}
my$of1=$ARGV[0];
$of1 =~ s/DE_results/DE_summary.txt/g;

my$of=$ARGV[1]?$ARGV[1]:$of1;
$of1=$ARGV[0];
die "\n\nInput:$ARGV[0]  and  Output:$of are same. Cannot overwrite\n\n" if $of eq $ARGV[0];
open my $OUT,">",$of;
##$of1="wAnnot.RSEM_express.matrix_2_trimmed_reads_OnLAP.US56.combined_NormBy.TMM_LeafStem.genes.counts.matrix.StemHighSV_vs_StemLowSV.DESeq2.DE_results"
$of1=~s/\S+\///g;
$of1=~s/wAnnot.*genes.counts.matrix.//g ;
$of1=~s/\.[^\.]+\.DE_results//g ;



print $OUT "File\tComparison\tDE>2x(pvalue<=0.05)\tUP(pvalue<=0.05)\tDown(pvalue<=0.05)\tDE>2x(qvalue<=0.05)\tUP(qvalue<=0.05)\tDown(qvalue<=0.05)";
print $OUT "\n$ARGV[0]\t$of1\t$tot_P\t$up_P\t$down_P\t$tot_Q\t$up_Q\t$down_Q";

close($OUT);
