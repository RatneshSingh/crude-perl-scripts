use strict;
use warnings;
use Getopt::Long;
our($evm,$weight,$out,$score_min,$score_max,$no_single,$outfile,$help,$verbose);


GetOptions ("evm=s" => \$evm,    # numeric
            "weight=s"   => \$weight,      # Weight file
	    "score_min|sm=i"   => \$score_min,      # Remove gene modeles with scores lower than this Minimum score 
	    "score_max|xm=i"   => \$score_max,      # Remove gene modeles with scores lower than this Maximum score
	    "no_single|ns"   => \$no_single,      # Exclde single exon gene models
            "outfile=s"   => \$outfile,      # Save filtered outputs.
	    "verbose"  => \$verbose,
	    "help"=>\$help
);  # 


my$usage="
perl $0 -options....
Options:
-evm	evm
-weight		Weight file
-score_min|sm	Remove gene modeles with scores lower than this Minimum score 
-score_max|xm	Remove gene modeles with scores lower than this Maximum score
-no_single|ns	Exclde single exon gene models
-outfile	Save filtered outputs.
-verbose	verbose
-help		Print usage einfo and exit; 
";

##################################################################
### Obituaries.
die "\n$usage\n" if $help;
die "No evm file found\n$usage\n" if !$evm;
die "No weight file detected\n$usage\n" if !$evm;

###################################################################
## defaults
$outfile=$outfile?$outfile:"$evm.filtered.out";

###################################################################
### Fatal Opening Ceremonies
open WEIGHT,"$weight" or die "Unable to open $weight";
open EVM,"$evm" or die "Unable to open $evm";
open OUT, ">$outfile" or die "unable to open output file";


##################################################################
print "Read filenames :\nWEIGHT:  $weight\nEVM output file $evm\n";
my %evcategory;
my%ev_cat;
my %summary;
while(<WEIGHT>){
s/^\s+//g;
next if m/^#/;  ## ignore commented lines.
next if m/^\s*$/;  ### ignore empty lines
my@cats=split /\s+/;
    $evcategory{$cats[1]}=$cats[0];
    $ev_cat{$cats[0]}=1;
}


print "Processing $evm\n";
$/="\n#";
my$block=0;
while(<EVM>){
    next if m/^!/;
    next if m/^\s*$/;
    #print "\nFound single: $_ " if ($no_single && m/\d+\s+\d+\s+single[\-\+]\s+/);
    next if ($no_single && m/\d+\s+\d+\s+single[\-\+]\s+/);
    print "\nThis line is not supposed to print if when single is found:" if ($no_single && m/\d+\s+\d+\s+single[\-\+]\s+/);
    chomp();
    $block++;   

    print "Finished Processing: $block blocks\n" if $block%1000 ==0;
    #print;# if m/^!/;
    my@block=split /\n/;
    my %evidence;
    my($Mode,$S_ratio,$coord,$orient,$score,$noncoding_equivalent,$raw_noncoding,$offset);              #(4.00)

    foreach my $line(@block){
	$line =~ s/^\s+//g;
	next if $line =~ m/^!|^#/;
	next if $line =~ m/^\s*$/;
	if ($line =~ m/^EVM\s+prediction\:/){
	    $Mode=$1 if $line =~ m/Mode\:\s*([\S]+)/;
	    $S_ratio=$1 if $line =~ m/S\-ratio\:\s*([\d\.]+)/;
	    $coord=$1 if $line =~ m/\s+([\d]+\-[\d]+)\s+orient/;
	    $orient=$1 if $line =~ m/\s+orient\s*\(([\+\-]+)\)\s*/;
	    $score=$1 if $line =~ m/\s+score\s*\(([\d\.]+)\)/;
	    $noncoding_equivalent=$1 if $line =~ m/\s+noncoding_equivalent\s*\(([\d\.]+)\)/;
	    $raw_noncoding=$1 if $line =~ m/\s+raw_noncoding\s*\(([\d\.]+)\)/;
	    $offset=$1 if $line =~ m/\s+offset\(([\d\.]+)\)/;
	    #print "\n*****Found values in: $line" if ($Mode);
	    #print join("\t","****\n",$Mode,$S_ratio,$coord,$orient,$score,$noncoding_equivalent,$raw_noncoding,$offset);
	}
	
	next if ($score_min && $score < $score_min);
	next if ($score_max && $score > $score_max);
	next if ($line =~ m/^EVM\s+prediction\:\s+Mode/);
	#print join("\t","\n****",$Mode,$S_ratio,$coord,$orient,$score,$noncoding_equivalent,$raw_noncoding,$offset);


	my@col=split /\s+/,$line;
	#print "\n\n\n*****This column have 2: $line ***\n" if $col[2];
	my @evms=split /,/,$col[5] if $col[2] !~ m/INTRON/;
    	@evms=split /,/,$col[3] if $col[2] =~ m/INTRON/;
	#print join "*","\nPrinting ind evd:",@evms,"\n";
	foreach my$ev(@evms){
	    next if $ev =~ m/^\s*$/;
	    $ev =~ s/\{|\}//g;
	    #print "\n*** New ev $ev";
	    my@enames=split /;/,$ev;
    	    print "\n**** Could not find $ev in line: $line \n IN BLOCK \n $_" if !$evcategory{$enames[-1]};
	    $evidence{$evcategory{$enames[-1]}}++;
	    $evidence{$enames[-1]}++;
	    #print "*** Processes $evname for category $evcategory{$evname}\n";
	}
	


    }
    #print "\nThis block has evidences from ";
    foreach my $catgs(keys %ev_cat, keys %evcategory){
	$summary{$catgs}{"NO"}++ if (!$evidence{$catgs} || $evidence{$catgs} ==0);
	$summary{$catgs}{"YES"}++ if ($evidence{$catgs} && $evidence{$catgs}> 0);
	#print "$catgs:",$evidence{$catgs}?$evidence{$catgs}:0;
    	#print "\t";
    }

    
    #foreach my $catgs(keys %evcategory){
    #    $summary{$catgs}{"NO"}++ if (!$evidence{$catgs} || $evidence{$catgs} ==0);
    #    $summary{$catgs}{"YES"}++ if ($evidence{$catgs} && $evidence{$catgs}> 0); 
    #    #print "$catgs:",$evidence{$catgs}?$evidence{$catgs}:0;
    #    #        #print "\t";
    #}     
print OUT "\n#$_" if ($evidence{TRANSCRIPT} && $evidence{TRANSCRIPT} > 0);


}


foreach my$catgs(keys %evcategory){
print "\nEntries with No $catgs evidences:",$summary{$catgs}{"NO"}?$summary{$catgs}{"NO"}:0;
print "\nEntries with $catgs evidences:", $summary{$catgs}{"YES"}?$summary{$catgs}{"YES"}:0;

}



print "\n\n\nSummary of Major Categories:";
foreach my$catgs(keys %ev_cat){
print "\nEntries with No $catgs evidences:",$summary{$catgs}{"NO"}?$summary{$catgs}{"NO"}:0;
print "\nEntries with $catgs evidences:", $summary{$catgs}{"YES"}?$summary{$catgs}{"YES"}:0;

}









