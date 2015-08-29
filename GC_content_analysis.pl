#!/usr/bin/perl -w
use strict;
use Getopt::Long;
#use PDL::Stats;
#use ScatterPlot;
use GD::Graph::points;
use Statistics::LineFit;
use Data::Dumper;
$| = 1;

our (
    $opt_m,      $opt_n, $opt_b,       $opt_s,        $opt_o,        $opt_d,    $opt_f,      $opt_r,
    $opt_l,      $opt_w, $opt_slide,   $trim_at_stop, $start_at_atg, $frame,    $count_this, $which_base,
    $num_codons, $group, @group_lower, @group_upper,  $min_length,,  $group_by, $opt_plot,   ,$opt_stat,
	$help
);

$opt_m        = 100;
$opt_b        = 1;
$opt_n        = 1;
$opt_d        = 'relative';
$opt_f        = 'gcbin';
$opt_r        = 50000;
$opt_l        = 1000;
$opt_slide    = 100;
$opt_w        = 100;
$start_at_atg = 'yes';
$trim_at_stop = 'yes';
$count_this   = 'GC';
$which_base   = 3;
$frame        = 1;
$num_codons   = 300;
$min_length   = 600;
$group_by     = "GC3";

#getopt('smnbofdrl');

my $result = GetOptions(
    "s=s"                      => \$opt_s,
    "m=i"                      => \$opt_m,
    "n=i"                      => \$opt_n,
    "b=i"                      => \$opt_b,
    "o=s"                      => \$opt_o,
    "f=s"                      => \$opt_f,
    "d=s"                      => \$opt_d,
    "r=i"                      => \$opt_r,
    "l=i"                      => \$opt_l,
    "w=i"                      => \$opt_w,
    "slide=i"                  => $opt_slide,
    "start_at_atg|sat=s"       => \$start_at_atg,
    "trim_at_stop|tas=s"       => \$trim_at_stop,
    "frame=i"                  => \$frame,
    "count_this|ct=s"          => \$count_this,
    "which_base|wb=i"          => \$which_base,
    "num_codon|num_nt|nc|nn=i" => \$num_codons,
    "group=i"                  => \$group,
    "group_lower|gl=f"         => \@group_lower,
    "group_upper|gu=f"         => \@group_upper,
    "group_by|gb=s"            => \$group_by,
    "min_length|ml=i"          => \$min_length,
	"plot"					   => \$opt_plot,
	"stat"					   => \$opt_stat,
    "help"                     => \$help

);

if(uc$opt_f eq "GC3_VS_GC"){$opt_f ="GC_VS_GC3"}

my $usage = '
This script can calculate GC content of all the sequences or individual sequences in one or multiple files.
It can also be used to calculate GC content of randomly selected sequences of specified length from one or multiple files.
GC bin result can be absolute frequency or relative frequency of GC contents normalized to the number of sequences in a file.


usage: perl script options file_names
-s	Sequences file. If using multiple files put the names at
	the end of command line after all the options .

-f	gcbin|randomgcbin|gcperfile|gcperseq|gc_sliding_window|GC3|GC3_positional|GC3_positional_grouped|codon_usage|GC_vs_length|GC3_vs_length|nt|GC3_grouped|GC_vs_GC3
	gcbin: Report Frequency (relative or absolute) of GC content of sequences in each bin
		based on Max(-m),Min(-n) and Binsize(-b).
	randomgcbin: Randomly pick the number(-r) sequence of specified (-l) size from the
		file and report GC content in each bin
	gcperfile: Report GC content of concatenated sequences in each file.
	gcperseq: Report the GC content for each sequence in a file.
    gc_sliding_window   Report GC content for each window of size(-w) which slides every (-slide) bases.
	GC3:	Calculate GC3 for coding sequences.
	GC3_positional:	Calculate Average GC3 for each codon starting from ATG in all sequence.
	GC3_positional_grouped:	Calculate Average GC3 for each codon starting from ATG in sequences in a group of low and high GC3 frequency.
	codon_usage: print codon usage table for all the species.
	GC_vs_length: calculate fraction of sequences in GC% bin per sequence length.
	NT:	Calculate frequency of any nucleotide at any position in a codons.
	GC3_grouped: GC3 distribution of sequences grouped by GC3 or GC content.
	GC_vs_GC3: Calculate GC content and GC3 content for each sequence.


++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
options for gcbin|randomgcbin|gcperfile|gcperseq|gc_sliding_window
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
-d	relative|absolute|. Report frequency normalized to number of sequences(relative) of absolute values [relative]
-m	Max bin [100].
-n	Min bin [1]
-b	Bin size [1]
-w  Window size[100]
-slide  step size for the sliding window[100]
-o	Output file name. If not provided program will generate the name automatically.
-r	Number of random sequences to be used for randomgcbin.
-l	Length of randomly selected sequence to be used for randomgcbin.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
options for GC3.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
-start_at_atg|sat	yes|no. Find ORF starting at ATG [yes]
-trim_at_stop|tas	yes|no. Trim ORF at stop codon if internal stop codons were found [yes]
-frame				1|2|3. Which frame to look into for GC3 calculations[1]
-count_this|ct			GC|AT. count frequency for either GC or AT.
-which_base|wb		1|2|3. Which base position to look at for counting.[3]
-min_length|ml		Min length of a cds to be used for GC3 analysis[600]


+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
options for NT. Frequency of A,T,G, or C in synonymous codons at any position.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
-start_at_atg|sat	yes|no. Find ORF starting at ATG [yes]
-trim_at_stop|tas	yes|no. Trim ORF at stop codon if internal stop codons were found [yes]
-frame				1|2|3. Which frame to look into for GC3 calculations[1]
-count_this			A|T|G|C. count frequency for A,T,G, or C or multiple combination [ATGC].
-which_base|wb		1|2|3. Which base position to look at for counting.[3]
-min_length|ml		Min length of a cds to be used for GC3 analysis[600]



+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
options for GC3_positional.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
-start_at_atg|sat	yes|no. Find ORF starting at ATG [yes]
-trim_at_stop|tas	yes|no. Trim ORF at stop codon if internal stop codons were found [yes]
-frame				1|2|3. Which frame to look into for GC3 calculations[1]
-count_this			GC|AT. count frequency for either GC or AT.
-which_base|wb		1|2|3. Which base position to look at for counting.[3]
-num_codons|nc		Till How many codons from ATG to print in GC3_positional analysis[300]
-min_length|ml		Min length of a cds to be used for GC3 analysis[600]


+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
options for GC3_positional_grouped.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
-start_at_atg|sat	yes|no. Find ORF starting at ATG [yes]
-trim_at_stop|tas	yes|no. Trim ORF at stop codon if internal stop codons were found [yes]
-frame				1|2|3. Which frame to look into for GC3 calculations[1]
-count_this			GC|AT. count frequency for either GC or AT.
-which_base|wb		1|2|3|123. Which base position to look at for counting.[3]
					123: means all the position.
-num_codons|nc		Till How many codons from ATG to print in GC3_positional analysis[300]
-group				How many group to create if range need to be entered manually.
					Program will ask for the lower and upper range GC/GC3 values for each group.
					use -group_lower and -group_upper for autmatic creation of groups instead..
-group_lower|gl		lower percent range of GC/GC3 groups. Accepts multiple values seperated with space e.g. 20 30 40.
-group_upper|gu		Upper percent range of GC/GC3 groups. Accepts multiple values seperated with space e.g 30 40 50.
-group_by			GC|GC3. Group sequences by GC or GC3.[GC3]
-min_length|ml		Min length of a cds to be used for GC3 analysis[600]



+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
options for GC3_grouped.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
-start_at_atg|sat	yes|no. Find ORF starting at ATG [yes]
-trim_at_stop|tas	yes|no. Trim ORF at stop codon if internal stop codons were found [yes]
-frame				1|2|3. Which frame to look into for GC3 calculations[1]
-count_this			GC|GC3|AT. count frequency for either GC,GC3 or AT.
-which_base|wb		1|2|3. Which base position to look at for counting.[3]
-group				How many group to create if range need to be entered manually.
					Program will ask for the lower and upper range GC/GC3 values for each group.
					use -group_lower and -group_upper for autmatic creation of groups instead..
-group_lower|gl		lower percent range of GC/GC3 groups. Accepts multiple values seperated with space e.g. 20 30 40.
-group_upper|gu		Upper percent range of GC/GC3 groups. Accepts multiple values seperated with space e.g 30 40 50.
-group_by			GC|GC3. Group sequences by GC or GC3.[GC3]
-min_length|ml		Min length of a cds to be used for GC3 analysis[600]

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
options for codon_usage
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
-d	relative|absolute|. Report absolute number of codons or relative fractions normalized to number of synonymous codons [relative]



+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
options for GC_vs_length|GC3_vs_length
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-num_nt|nn	Calculate the data for sequence till this base number [300].



+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
options for GC3_vs_length
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
-start_at_atg|sat	yes|no. Find ORF starting at ATG [yes]
-trim_at_stop|tas	yes|no. Trim ORF at stop codon if internal stop codons were found [yes]
-frame				1|2|3. Which frame to look into for GC3 calculations[1]
-count_this			GC|AT. count frequency for either GC or AT.
-which_base|wb		1|2|3. Which base position to look at for counting.[3]
-num_nt|nn	Maximum sequence length to calculate GC3 for [300].


+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
options for GC_vs_GC3
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
-plot	Plot Scatter plot for each species[FALSE].
-stat	Write statistics fo correlation in a file [FALSE]



';

die "\n$usage\n" if $help;

my ( @GC_range, %GC, @files );
#############################################################################################################################
# Check for files or multiple files and create the array of the file name.
if ( !$opt_s && !$ARGV[0] ) { die "\n\n\nError:Please provide sequence file in fasta format\n\n$usage\n"; }

# add filenames into an array @files.
else {
    if ( $ARGV[0] ) {
        foreach (@ARGV) { $_ =~ s/\s+//g; push( @files, $_ ) }
    }
    unshift( @files, $opt_s ) if $opt_s;
}

#########################################################################################################################
# Create output filenames if not provided by user.
my $other_files_num = scalar @files - 1;
my $filename;
if   ($opt_s) { $filename = getfilename( remove_extension($opt_s) ) }
else          { $filename = getfilename( remove_extension( $ARGV[0] ) ) }
# temp string to add group info in output file name
my$temp_string="_";
for(my$i=0; $i < min( scalar @group_lower, scalar @group_upper ); $i++){$temp_string.="$group_lower[$i]-$group_upper[$i]"; $temp_string.="_";}


if ( !$opt_o ) {
    if ( scalar @files > 1 ) {
        $opt_o = $opt_f . "_" . $filename . "_and_" . $other_files_num . "_other_sequences.$opt_d.table";
        $opt_o = $opt_f . "_by_$group_by" .$temp_string. "_" . $filename . "_and_" . $other_files_num . "_other_sequences.$opt_d.table" if uc $opt_f eq "GC3_POSITIONAL_GROUPED"||uc $opt_f eq "GC3_GROUPED";
    }
    else {
        $opt_o = $opt_f . "_" . $filename . $opt_d . ".table";
        $opt_o = $opt_f . "_by_$group_by" .$temp_string."_" . $filename . $opt_d . ".table" if uc $opt_f eq "GC3_POSITIONAL_GROUPED"||uc $opt_f eq "GC3_GROUPED";
    }
}
print "\nCalculating GC contents for sequences from file:" . join( "\n", @files ) . "\n";

############################################################################################################################
#open output file and print header and create GC range if asked for gcbin
open OUT, ">$opt_o";
if (   uc $opt_f eq "GCBIN"
    || uc $opt_f eq "RANDOMGCBIN"
    || uc $opt_f eq "GC_SLIDING_WINDOW"
    || uc $opt_f eq "GC3"
    || uc $opt_f eq "GC3_POSITIONAL"
    || uc $opt_f eq "GC3_POSITIONAL_GROUPED"
    || uc $opt_f eq "NT"
    || uc $opt_f eq "GC3_GROUPED")
{

    if ( uc $opt_f eq "GC3_POSITIONAL" || uc $opt_f eq "GC3_POSITIONAL_GROUPED" ) { $opt_m = $num_codons; $opt_b = 1; $opt_n = 1; }
    for ( my $i = $opt_n; $i <= $opt_m; $i = $i + $opt_b ) { push( @GC_range, $i ); print OUT "\t$i"; print "\t$i"; }
}
elsif ( uc $opt_f eq "GCPERFILE" || uc $opt_f eq "GCPERSEQ" ) { print OUT "\tGC"; print "\tGC"; }

if(uc$opt_f eq "GC_VS_GC3"){open OUT2,">Stat_$opt_o"; print OUT2 "\tGC\tintercept\tslope\trSquared\tmeanSquaredError\tdurbinWatson\tsigma\ttStatIntercept\ttStatSlope\tvarianceIntercept\tvarianceSlope"}
############################################################################################################################
# create group labels if GC3_positional_grouped was selected as option.
if ( uc $opt_f eq "GC3_POSITIONAL_GROUPED" ) {

    if ( @group_lower && @group_upper ) {
        if ( scalar @group_lower ne scalar @group_upper ) {
            print "\n****Warning: Number of values in lower and upper are not equal. Only following " . min( scalar @group_lower, scalar @group_upper ) . " Groups will be created";
            for ( my $i = 0; $i < min( scalar @group_lower, scalar @group_upper ); $i++ ) {
                print "\nGroup", $i + 1, ": $group_lower[$i] - $group_upper[$i]";
            }
            print "\n";
        }

    }
    elsif ($group) {
        for ( my $i = 1; $i <= $group; $i++ ) {

            # collect lower range for group.
            my $entry = 0;
            print "\nEnter Lower range for group $i:";
            $entry = <STDIN>;
            push( @group_lower, $entry );

            # collect upper range for group.
            $entry = 0;
            print "\nEnter upper range for group $i:";
            $entry = <STDIN>;
            push( @group_upper, $entry );

        }
    }
    else { $group_lower[0] = 0; $group_upper[0] = 100; }
}

############################################################################################################################
# calculate GC for each file.
############################################################################################################################

# keep the codon usage hash outside of the loop so it can accumulate al the values in all the files. will help to print later.
my %codon_usage;
my %GC_vs_length;
my $gcvslen_max_len = 1;
my @nt = split( //, $count_this );
foreach (@nt) { $_ =~ s/[^\w]//g }

#print "\nintercept\tslope\trSquared\tmeanSquaredError\tdurbinWatson\tsigma\ttStatIntercept\ttStatSlope\tvarianceIntercept\tvarianceSlope" if uc$opt_f eq "GC_VS_GC3";

foreach my $file (@files) {

    # keep following hashes inside of the foreach loop so that they empty out everytime new file is read.
    #keeping them out of the loop will add on the values from previous fles to the next one.
    my %GC3_positional;
    my %GC3_positional_grouped;
    my @GC_contents;
    my $print = 'no';
    #######################################################################################
    # reset the counting of bins for each file to 0
    for ( my $i = $opt_n; $i <= $opt_m; $i = $i + $opt_b ) {
        $GC{$file}{$i} = 0 if ( uc $opt_f ne "GC3_POSITIONAL_GROUPED" && uc $opt_f ne "NT" );
        if ( uc $opt_f eq "NT" ) {
            foreach my $base (@nt) { $GC{$file}{$base}{$i} = 0 }
        }
    }

    # get sequence hash for gcbin and randomgcbin
    my $hSequence = ReadFastaTohash( $file, 1, $opt_f );
    my $current_filename = remove_extension( getfilename($file) );
    ################################################################################################
    # if gcbin or randomgcbin was asked. calculate GC% and store the count in %GC{$filename}{$gcbin}
    if ( uc $opt_f eq "GCBIN" || uc $opt_f eq "RANDOMGCBIN" ) {

        # sort array in numerical order
        #@GC_range = sort { $a <=> $b } @GC_range;
        #my $count = 0;

        # Calculate GC content for each sequence
        #my @GC_contents;
        foreach ( keys %{$hSequence} ) {
            my $GC_content = GC_content_percent( $$hSequence{$_} );
            push( @GC_contents, $GC_content ) if $GC_content > 0;
            for ( my $i = $opt_n; $i <= $opt_m; $i = $i + $opt_b ) {
                if ( $GC_content <= $i && $GC_content > $i - $opt_n ) { $GC{$file}{$i}++; }
            }
        }
        $print = 'yes';
    }

    elsif ( uc $opt_f eq "GC_VS_LENGTH" || uc $opt_f eq "GC3_VS_LENGTH" ) {

        # Calculate GC content for each sequence
        #my @GC_contents;
		my$num_seq=0;
        foreach ( keys %{$hSequence} ) {
            my $GC_content = GC_content_percent( $$hSequence{$_} ) if uc $opt_f eq "GC_VS_LENGTH";
            $GC_content = count_GCn( $$hSequence{$_}, $frame, $start_at_atg, $trim_at_stop, $count_this, $which_base ) if uc $opt_f eq "GC3_VS_LENGTH";
            my $seq_length = length( $$hSequence{$_} );
            push( @GC_contents, $GC_content ) if $GC_content > 0;
            push( @{ $GC_vs_length{$current_filename}{$seq_length}{GC} }, $GC_content );
			$num_seq++;
        }

        foreach my $seq_length ( keys %{ $GC_vs_length{$current_filename} } ) {

            $GC_vs_length{$current_filename}{$seq_length}{average} = average( \@{ $GC_vs_length{$current_filename}{$seq_length}{GC} } );
            $GC_vs_length{$current_filename}{$seq_length}{stdeve}  = stdev( \@{ $GC_vs_length{$current_filename}{$seq_length}{GC} } );
			$GC_vs_length{$current_filename}{$seq_length}{proportion}  = scalar (@{ $GC_vs_length{$current_filename}{$seq_length}{GC} })*100/$num_seq;
            $gcvslen_max_len = $seq_length if $seq_length > $gcvslen_max_len;
        }

    }

	#===========================================================================================================================
	################################################################################################
    # if gcbin or randomgcbin was asked. calculate GC% and store the count in %GC{$filename}{$gcbin}
    elsif ( uc $opt_f eq "GC_VS_LENGTH_GROUPED" || uc $opt_f eq "GC3_VS_LENGTH_GROUPED" ) {

        # Calculate GC content for each sequence
        #my @GC_contents;
        foreach ( keys %{$hSequence} ) {
            my $GC_content = GC_content_percent( $$hSequence{$_} ) if uc $opt_f eq "GC_VS_LENGTH_GROUPED";
            $GC_content = count_GCn( $$hSequence{$_}, $frame, $start_at_atg, $trim_at_stop, $count_this, $which_base ) if uc $opt_f eq "GC3_VS_LENGTH_GROUPED";
            my $seq_length = length( $$hSequence{$_} );
            push( @GC_contents, $GC_content ) if $GC_content > 0;
            push( @{ $GC_vs_length{$current_filename}{$seq_length}{GC} }, $GC_content );
        }

        foreach my $seq_length ( keys %{ $GC_vs_length{$current_filename} } ) {

            $GC_vs_length{$current_filename}{$seq_length}{average} = average( \@{ $GC_vs_length{$current_filename}{$seq_length}{GC} } );
            $GC_vs_length{$current_filename}{$seq_length}{stdeve}  = stdev( \@{ $GC_vs_length{$current_filename}{$seq_length}{GC} } );
			$GC_vs_length{$current_filename}{$seq_length}{proportion}  = scalar @{ $GC_vs_length{$current_filename}{$seq_length}{GC} };


            $gcvslen_max_len = $seq_length if $seq_length > $gcvslen_max_len;
        }

    }

	#===========================================================================================================================






    #################################################################################################
    # if gc_sliding_window was asked, calculate GC% in each window and store info in %GC{$filename}{$gcbin}
    elsif ( uc $opt_f eq "GC_SLIDING_WINDOW" ) {

        # sort array in numerical order
        @GC_range = sort { $a <=> $b } @GC_range;

        my $count = 0;

        # Calculate GC content for each substr of sequence
        #my @GC_contents;
        foreach ( keys %{$hSequence} ) {
            my $seq_len = length( $$hSequence{$_} );

            for ( my $j = 0; $j <= $seq_len - $opt_w - 1; $j = $j + $opt_slide ) {

                my $sebseq = substr( $$hSequence{$_}, $j, $opt_w );
                my $GC_content = GC_content_percent($sebseq);
                push( @GC_contents, $GC_content ) if $GC_content > 0;
                for ( my $i = $opt_n; $i <= $opt_m; $i = $i + $opt_b ) {
                    if ( $GC_content <= $i && $GC_content > $i - $opt_n ) { $GC{$file}{$i}++; }
                }
            }
        }
        $print = 'yes';
    }
    ####################################################
    # for calculating GC3 in coding sequences
    elsif ( uc $opt_f eq "GC3" ) {

        # sort array in numerical order
        #@GC_range = sort { $a <=> $b } @GC_range;
        #my $count = 0;

        # Calculate GC content for each sequence
        foreach ( keys %{$hSequence} ) {
            next if $_ =~ /^\s*$/;
            next if length( $$hSequence{$_} ) < $min_length;
            my $GC_content = count_GCn( $$hSequence{$_}, $frame, $start_at_atg, $trim_at_stop, $count_this, $which_base );
            push( @GC_contents, $GC_content ) if $GC_content > 0;
            for ( my $i = $opt_n; $i <= $opt_m; $i = $i + $opt_b ) {
                if ( $GC_content <= $i && $GC_content > $i - $opt_n ) { $GC{$file}{$i}++; }
            }
        }
        $print = 'yes';
    }

    ####################################################
    # for calculating NT in coding sequences
    elsif ( uc $opt_f eq "NT" ) {

        my %num_seq;

        # Calculate NT content for each sequence
        foreach my $bases (@nt) {
            foreach ( keys %{$hSequence} ) {
                next if $_ =~ /^\s*$/;
                next if length( $$hSequence{$_} ) < $min_length;
                my $nt_content = count_nt( $$hSequence{$_}, $frame, $start_at_atg, $trim_at_stop, $bases, $which_base );
                push( @{ $num_seq{$bases} }, $nt_content ) if $nt_content > 0;
                for ( my $i = $opt_n; $i <= $opt_m; $i = $i + $opt_b ) {
                    if ( $nt_content <= $i && $nt_content > $i - $opt_n ) { $GC{$file}{$bases}{$i}++; }
                }
            }
            $print = 'no';
        }

        # print the result for this file.
        my $filename = remove_extension( getfilename($file) );

        #print "Printing GC frequencies to output file\n";

        foreach my $base ( keys %{ $GC{$file} } ) {
            my $num_GC_contents = scalar @{ $num_seq{$base} };
            print OUT "\n$filename/$base";
            print "\n$filename/$base";
            for ( my $i = $opt_n; $i <= $opt_m; $i = $i + $opt_b ) {

                printf OUT "\t%0.3f", $GC{$file}{$base}{$i} / $num_GC_contents if lc $opt_d eq 'relative';
                printf "\t%0.3f", $GC{$file}{$base}{$i} / $num_GC_contents if lc $opt_d eq 'relative';
                printf OUT "\t%0.3f", $GC{$file}{$base}{$i} if lc $opt_d eq 'absolute';
                printf "\t%0.3f", $GC{$file}{$base}{$i} if lc $opt_d eq 'absolute';
            }
        }
        print "\nFinished Calculating $opt_f for file $file.\n"

    }

    ####################################################
    # for calculating GC3 in coding sequences
    elsif ( uc $opt_f eq "GC3_POSITIONAL" ) {

        # sort array in numerical order
        #@GC_range = sort { $a <=> $b } @GC_range;
        # Collect codons from each position of each sequence in an array and put it as a value of Hash.
        # HASH WILL BE PROCESSED BELOW.
        foreach ( keys %{$hSequence} ) {
            next if $_ =~ /^\s*$/;

            my $orf = get_ORF( $$hSequence{$_}, 1, 'yes', 'yes' );
            next if length($orf) < $min_length;
            collect_positional_codons( $orf, $frame, \%GC3_positional ) if $which_base <=3;
			collect_positional_NT( $orf, \%GC3_positional ) if $which_base >3;
        }


		foreach my $base_positions ( sort { $a <=> $b } keys %GC3_positional ) {

            # create virtusl sequence from collected codons and count GCn in that.
            # my$virt_seq=join(@{$GC3_positional{$base_positions}},"");
            my $virt_seq = join( "", @{ $GC3_positional{$base_positions} } );

            # print "\nbase_positions:$base_positions\nVirtual Sequence:$virt_seq\n";#print (@{$hash{$_}})
            my $GC_content = count_GCn( $virt_seq, $frame, 'no', 'no', $count_this, $which_base )  if $which_base <=3;
		    $GC_content = GC_content_percent($virt_seq)  if $which_base <=3;
            push( @GC_contents, $base_positions );
            $GC{$file}{$base_positions} = $GC_content;

            $print = 'no';
        }

	}

    ####################################################
    # for calculating GC3 in coding sequences
    elsif ( uc $opt_f eq "GC3_POSITIONAL_GROUPED" ) {

        #get the group with smallest number of values.
        my $group_length = min( scalar @group_lower, scalar @group_upper );

        #print "\nGroup_length minimum:$group_length\n";
        # sort array in numerical order
        #@GC_range = sort { $a <=> $b } @GC_range;
        my $count = 0;

        # Collect codons from each position of each sequence in an array and put it as a value of Hash.
        # HASH WILL BE PROCESSED BELOW.
        foreach ( keys %{$hSequence} ) {
            next if $_ =~ /^\s*$/;
            my $orf = get_ORF( $$hSequence{$_}, 1, 'yes', 'yes' );
            next if length($orf) < $min_length;

            # use GC percent as criteria to group if $group_by is GC else use GC3.
            my $current_GCn = count_GCn( $orf, $frame, 'no', 'no', $count_this, 3 );
            $current_GCn = GC_content_percent($orf) if uc $group_by eq "GC";

            for ( my $i = 0; $i < $group_length; $i++ ) {

                #print "i:$i\t";
                if ( $current_GCn > $group_lower[$i] && $current_GCn <= $group_upper[$i] ) {
                    collect_positional_codons_grouped( $orf, $frame, \%GC3_positional_grouped, $group_lower[$i], $group_upper[$i] ) if $which_base <=3 ;
					collect_positional_NT_grouped( $orf, \%GC3_positional_grouped, $group_lower[$i], $group_upper[$i] ) if $which_base >3;

                }
            }
        }

		# process and print
		foreach my $group_number ( keys %GC3_positional_grouped ) {
            foreach my $base_positions ( sort { $a <=> $b } keys %{ $GC3_positional_grouped{$group_number} } ) {

                # print "$group_number:$base_positions\t";
                # create virtusl sequence from collected codons and count GCn in that.
                # my$virt_seq=join(@{$GC3_positional{$base_positions}},"");
                my $virt_seq = join( "", @{ $GC3_positional_grouped{$group_number}{$base_positions} } );

                # print "\nbase_positions:$base_positions\nVirtual Sequence:$virt_seq\n";#print (@{$hash{$_}})
                my $GC_content = count_GCn( $virt_seq, $frame, 'no', 'no', $count_this, $which_base )if $which_base <=3 ;
				my $GC_content = GC_content_percent($virt_seq)if $which_base >3 ;
                #push( @GC_contents, $base_positions );
                #print "Group_number:$group_number\tBase_positions:$base_positions\tGC content:$GC_content\n";
                $GC{$file}{$group_number}{$base_positions} = $GC_content;

                #print "GC{file}{group_number}{base_positions}:$GC{$file}{$group_number}{$base_positions}\n";

                $print = 'no';
            }
        }



    }

    ####################################################
    # for calculating GC3 in coding sequences grouped by GC3 or GC
    elsif ( uc $opt_f eq "GC3_GROUPED" ) {
        my $group_length = min( scalar @group_lower, scalar @group_upper );
        my $count = 0;

        # Collect codons from each position of each sequence in an array and put it as a value of Hash.
        # HASH WILL BE PROCESSED BELOW.
        foreach ( keys %{$hSequence} ) {
            next if $_ =~ /^\s*$/;
            my $orf = get_ORF( $$hSequence{$_}, 1, 'yes', 'yes' );
            next if length($orf) < $min_length;

            # use GC percent as criteria to group if $group_by is GC else use GC3.
			my $current_GCn = GC_content_percent($$hSequence{$_});
			   $current_GCn = count_GCn( $orf, $frame, $start_at_atg, $trim_at_stop, $count_this, $frame )if uc $group_by eq "GC3";


            for ( my $i = 0; $i < $group_length; $i++ ) {

                #print "i:$i\t";
                if ( $current_GCn > $group_lower[$i] && $current_GCn <= $group_upper[$i] ) {
                    count_GCn_grouped( $orf, $frame, \%GC3_positional_grouped, $group_lower[$i], $group_upper[$i], $count_this );

                }
            }
        }


		# process and print.
		 my $num_GC_contents = 0;
        foreach my $group_number ( keys %GC3_positional_grouped ) {

			# initialize hash elements to zero.
            for ( my $i = $opt_n; $i <= $opt_m; $i = $i + $opt_b ) {
                    $GC{ $file . $group_number }{$i} = 0;
            }

			foreach ( @{ $GC3_positional_grouped{$group_number} } ) {
				my $GC_content = $_;
				$num_GC_contents++ if $GC_content > 0;
				for ( my $i = $opt_n; $i <= $opt_m; $i = $i + $opt_b ) {
					if ( $GC_content <= $i && $GC_content > $i - $opt_n ) { $GC{ $file . $group_number }{$i}++; }
				}
			}
		}


		foreach my $group_number ( keys %GC3_positional_grouped ) {

			#print "Printing GC frequencies to output file\n";
			my $tempfilename = remove_extension( getfilename($file) );
			print OUT "\n$tempfilename.$group_number";
			print "\n$tempfilename.$group_number";
			for ( my $i = $opt_n; $i <= $opt_m; $i = $i + $opt_b ) {

				printf OUT "\t%0.3f", $GC{ $file . $group_number }{$i} / $num_GC_contents if lc $opt_d eq 'relative';
				printf "\t%0.3f", $GC{ $file . $group_number }{$i} / $num_GC_contents if lc $opt_d eq 'relative';
				printf OUT "\t%0.3f", $GC{ $file . $group_number }{$i} if lc $opt_d eq 'absolute';
				printf "\t%0.3f", $GC{ $file . $group_number }{$i} if lc $opt_d eq 'absolute';
			}
		}
		$print = 'no';

    }

    #####################################################
    # For calculating codon usage
    elsif ( uc $opt_f eq "CODON_USAGE" ) {

        foreach ( keys %{$hSequence} ) {
            next if $_ =~ /^\s*$/;
            my $orf = get_ORF( $$hSequence{$_}, 1, 'yes', 'yes' );
            next if length($orf) < $min_length;
            count_codon( $current_filename, $$hSequence{$_}, 1, \%codon_usage, 'no' );

        }

    }
    # Calculate GC and GC3 for each sequence
	elsif(uc$opt_f eq "GC_VS_GC3"){

		my@GC_cont=();
		my@GC3_cont=();

		foreach ( keys %{$hSequence} ) {
            next if $_ =~ /^\s*$/;
			my$orf=get_ORF($$hSequence{$_},$frame,$start_at_atg,$trim_at_stop);

            next if length( $orf ) < $min_length;
            my $GC3_content = count_GCn($orf, $frame, $start_at_atg, $trim_at_stop, $count_this, $which_base );
			my $GC_content  = GC_content_percent($orf);
			push( @GC_cont, $GC_content) if $GC_content > 0 && $GC3_content >0;
			push( @GC3_cont, $GC3_content) if $GC_content > 0 && $GC3_content >0;

			push( @GC_contents, [$GC_content,$GC3_content]) if $GC_content > 0 && $GC3_content >0;
			#push( @GC_contents, $GC3_content ) if $GC3_content > 0;

         }



		print OUT "\n$file GC";
		foreach(@GC_cont){print OUT "\t$_"}
		print OUT "\n$file GC3";
		foreach(@GC3_cont){print OUT "\t$_"}

	if(uc$opt_plot){
		# plot using chart::plot
		#use Chart::Plot_RSMod;
		use Plot_RSMod2;
		my$image_width=400;
		my$image_height=300;
		my$xmin=0;
		my$ymin=0;
		my$xmax=100;
		my$ymax=100;

		my $img = Plot_RSMod2->new();
		my $anotherImg = Plot_RSMod2->new ($image_width, $image_height);

		$img->setData (\@GC_cont,\@GC3_cont,'Points Noline Black') or die( $img->error() );
		#$img->setData (\@xdataset, \@ydataset);
		#$img->setData (\@anotherdataset, 'red_dashedline_points');
		#$img->setData (\@xanotherdataset, \@yanotherdataset,
		#               'Blue SolidLine NoPoints');

		($xmin, $ymin, $xmax, $ymax) = $img->getBounds();
		my($name,@rest)=split(/_/,$current_filename);
		my$gn=substr($name,0,1);
		my$sp=substr($name,1);
		$img->setGraphOptions ('horGraphOffset' => 25,
								'vertGraphOffset' => 25,
								'title' => "$gn. $sp GC vs GC3",
								'horAxisLabel' => '%GC',
								'vertAxisLabel' => '%GC3' );

		open(IMG, ">$current_filename.GC_vs_GC3.png") or die $!;
		binmode IMG;
		print IMG $img->draw();
		close IMG;
	}
	# calculate R square of correlation coefficient.
	if($opt_stat){
		my$lineFit = Statistics::LineFit->new();
		$lineFit->setData (\@GC_cont,\@GC3_cont) or die "Invalid data";
		my($intercept, $slope) = $lineFit->coefficients();
		defined $intercept or die "Can't fit line if x values are all equal";
		my$rSquared = $lineFit->rSquared();
		my$meanSquaredError = $lineFit->meanSqError();
		my$durbinWatson = $lineFit->durbinWatson();
		my$sigma = $lineFit->sigma();
		my($tStatIntercept, $tStatSlope) = $lineFit->tStatistics();
		my@predictedYs = $lineFit->predictedYs();
		my@residuals = $lineFit->residuals();
		my($varianceIntercept, $varianceSlope) = $lineFit->varianceOfEstimates();
		# GC per file
		my$fileGC_content=roundoff(GC_per_file($file),2);
		print OUT2 "\n$current_filename\t$fileGC_content\t$intercept\t$slope\t$rSquared\t$meanSquaredError\t$durbinWatson\t$sigma\t$tStatIntercept\t$tStatSlope\t$varianceIntercept\t$varianceSlope";

	}


	}



	####################################################
    # else calculate the GC content for each sequence seperately in a file.
    else {

        #die "-f $opt_f function is not correct function to perform.";
        foreach ( keys %{$hSequence} ) {

            my $GC_content = GC_content_percent( $$hSequence{$_} );
            print OUT "\n$_";
            print "\n$_";
            printf OUT "\t%0.3f", $GC_content;
            printf "\t%0.3f", $GC_content;
        }

    }

###################################################################################################
# Print the rsults depending on what options were selected
###################################################################################################
# print the GC content in bin in output file.
if ( $print eq "yes" ) {

    my $num_GC_contents = scalar @GC_contents;
    my $filename        = remove_extension( getfilename($file) );

    #print "Printing GC frequencies to output file\n";
    print OUT "\n$filename";
    print "\n$filename";
    for ( my $i = $opt_n; $i <= $opt_m; $i = $i + $opt_b ) {

        printf OUT "\t%0.3f", $GC{$file}{$i} / $num_GC_contents if lc $opt_d eq 'relative';
        printf "\t%0.3f", $GC{$file}{$i} / $num_GC_contents if lc $opt_d eq 'relative';
        printf OUT "\t%0.3f", $GC{$file}{$i} if lc $opt_d eq 'absolute';
        printf "\t%0.3f", $GC{$file}{$i} if lc $opt_d eq 'absolute';
    }
}

##########################################################
# print the result for GC3_window calculations.
if ( uc $opt_f eq "GC3_POSITIONAL" ) {

    #my $num_base_positions = scalar @GC_contents;
    my $filename1 = remove_extension( getfilename($file) );

    #print "Printing GC frequencies to output file\n";
    print OUT "\n$filename1";
    print "\n$filename1";
    @GC_contents = sort { $a <=> $b } @GC_contents;
    my $count = 0;
    foreach (@GC_contents) {

        #printf OUT "\t%0.3f", $GC{$file}{$_} / $num_base_positions if lc $opt_d eq 'relative';
        #printf "\t%0.3f", $GC{$file}{$_} / $num_base_positions if lc $opt_d eq 'relative';
        printf OUT "\t%0.3f", $GC{$file}{$_};    # if lc $opt_d eq 'absolute';
        printf "\t%0.3f", $GC{$file}{$_};        # if lc $opt_d eq 'absolute';
        last if $count >= $opt_m;
        $count++;
    }
}

##########################################################
# print the result for GC3_window calculations.
if ( uc $opt_f eq "GC3_POSITIONAL_GROUPED" ) {

    #my $num_base_positions = scalar @GC_contents;
    my $filename1 = remove_extension( getfilename($file) );

    #print "Printing GC frequencies to output file\n";
    foreach my $group_number ( keys %{ $GC{$file} } ) {

        #print "Group_number:$group_number\n";
        print OUT "\n$filename1.$group_number";
        print "\n\n\n$filename1.$group_number";
        my $count = 0;
        foreach my $base_positions ( sort { $a <=> $b } keys %{ $GC{$file}{$group_number} } ) {
            next if $count >= $opt_m;

            #print "\t$base_positions:";
            printf OUT "\t%0.3f", $GC{$file}{$group_number}{$base_positions};    # if lc $opt_d eq 'absolute';
            printf "\t%0.3f", $GC{$file}{$group_number}{$base_positions};        # if lc $opt_d eq 'absolute';

            $count++;

        }
    }
}
}

################################################################################################
# print result for codon usage
if ( uc $opt_f eq "CODON_USAGE" ) {

  my $Ref_codon_hash = codon_hash();

  # create hash with the total number of synonymous codons per amino acid for each file.
  my %count_synonymous;
  my %codon_array;
  my %file_array;
  foreach my $codon ( sort keys %codon_usage ) {
      $codon =~ s/\s+//g;
      $codon_array{$codon} = 1;
      foreach my $file_codon ( sort keys %{ $codon_usage{$codon} } ) {
          $file_codon =~ s/\s+//g;
          $file_array{$file_codon} = 1;

          #print "Keys: $_\n";
          $count_synonymous{$file_codon}{ $$Ref_codon_hash{$codon}{amino} } += $codon_usage{$codon}{$file_codon}{1}{count};

      }
  }

  my @codon_array;
  my @file_array;

  foreach ( sort keys %codon_array ) { push( @codon_array, $_ ) }
  foreach ( sort keys %file_array )  { push( @file_array,  $_ ) }

  #print the header line
  print OUT "Codon\tAmino_Acid\tExpected_frequency";
  foreach my $file_codon (@file_array) { print OUT "\t$file_codon"; }

  #print the codon fraction/frequency
  foreach my $codon (@codon_array) {
      $codon =~ s/\s+//g;
      print OUT "\n$codon\t$$Ref_codon_hash{$codon}{amino}\t$$Ref_codon_hash{$codon}{exp_freq}";
      foreach my $file_codon (@file_array) {
          $file_codon =~ s/\s+//g;

          #print "Num_Codon and Synon_Num_codon \t$codon_usage{$codon}{$file_codon}{1}{count}\t$count_synonymous{$file_codon}{ $$Ref_codon_hash{$codon}{amino} }\n";
          print OUT "\t" . roundoff( $codon_usage{$codon}{$file_codon}{1}{count} / $count_synonymous{$file_codon}{ $$Ref_codon_hash{$codon}{amino} }, 2 )
            if uc $opt_d eq "RELATIVE";
          print OUT "\t" . $codon_usage{$codon}{$file_codon}{1}{count} if uc $opt_d eq "ABSOLUTE";

      }
  }
}

#_________________________________________________________________________________________________________________________________________

################################################################################################
# print result for GC_VS_LENGTH
if ( uc $opt_f eq "GC_VS_LENGTH" || uc $opt_f eq "GC3_VS_LENGTH" ) {
  print "\t\t";

  foreach ( 1 .. $num_codons ) {
      print OUT "\t$_";
  }

  foreach my $temp_filename ( sort keys %GC_vs_length ) {
      print OUT "\n$temp_filename\tAverage_GC";
      foreach my $temp_seqlen ( 1 .. $num_codons ) {
          if ( $GC_vs_length{$temp_filename}{$temp_seqlen}{average} ) { print OUT "\t$GC_vs_length{$temp_filename}{$temp_seqlen}{average}" }
          else { print OUT "\t0" }
      }
      print OUT "\n$temp_filename\tStdev_GC";
      foreach my $temp_seqlen ( 1 .. $num_codons ) {
          if ( $GC_vs_length{$temp_filename}{$temp_seqlen}{average} ) { print OUT "\t$GC_vs_length{$temp_filename}{$temp_seqlen}{stdeve}" }
          else { print OUT "\t0" }
      }

	  print OUT "\n$temp_filename\tProportion_GC";
      foreach my $temp_seqlen ( 1 .. $num_codons ) {
          if ( $GC_vs_length{$temp_filename}{$temp_seqlen}{proportion} ) { print OUT "\t$GC_vs_length{$temp_filename}{$temp_seqlen}{proportion}" }
          else { print OUT "\t0" }
      }


  }
}

print "\nFinished calculation.\nResult has been saved to fie:$opt_o\n\n";

#__________________________________________________________________________________________________________________________________________________

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###################################################################################
# subroutines
###################################################################################
sub ReadFastaTohash {    # returns reference to hash header as key and sequence as value

    my $seqfile     = shift;
    my $amb         = shift;
    my $func        = shift;
    my $concatenate = 'no';
    if ( lc $func eq "randomgcbin" || lc $func eq "gcperfile" ) { $concatenate = 'yes'; }
    my ( $header, @sequence );
    chomp $seqfile;
    open FASTA, "$seqfile";

    #print "\nreading Sequences from input file.....Plz wait...\n";
    my %seq_hash = ();
    $/ = "\n>";    # Change record seperator to read Fasta
    if ( $concatenate eq 'no' ) {
        my $last_N = 1;
        while (<FASTA>) {
            chomp;
            ( $header, @sequence ) = split( "\n", $_ );

            $header =~ s/>//;       # Remove Leading > from Header
            $header =~ s/\s*$//;    # Remove trailing spaces from header
            $header =~ s/^\s*//;    # Remove Leading spaces from Header

            my $sequence = join( "", @sequence );
            $sequence =~ s/\s//g;
            $sequence =~ s/\n//g;
            if ( $amb == 1 ) { $sequence =~ s/[^ATGC]//gi; }
            if ( $header =~ /^\s*$/ ) { next; }
            ############################################################################################################
            # Make sure no duplicate header exists else they will be overwritten and lost during hashing due to same key.
            if ( !exists $seq_hash{$header} ) {
                $seq_hash{$header} = $sequence;    #feed headers and sequences in hash.
                                                   #$seq_hash{'RS_Concatenated'}=$seq_hash{'RS_Concatenated'}.$Ns.$sequence;
            }
            else {

                # find a uniq header name by adding a number at the end. If header still exists, increase the number by one
                while ( exists $seq_hash{$header} ) { $header = $header . $last_N; $last_N++; }

                $seq_hash{$header} = $sequence;
            }
            ############################################################################################################
        }

        my @seq_count = keys(%seq_hash);
        my $seq_count = @seq_count;

        #print "Done....\nNumber of sequences read form input file = $seq_count\n\n";
        @seq_count = ();
        $/         = "\n";    # Record seperator set back to default newline.
    }

    else {
        my $sequence = ConcatenateFasta( $seqfile, $amb );    # concatenate sequence
        my $rand_range = length($sequence) - $opt_l;

        if ( $func eq "gcperfile" ) { $seq_hash{$seqfile} = $sequence; }

        elsif ( $func eq "randomgcbin" ) {

            # create $opt_r number of sequences of $opt_l length
            for ( my $i = 1; $i <= $opt_r; $i++ ) {
                my $random_number = int( rand($rand_range) );
                my $sub_seq       = substr( $sequence, $random_number, $opt_l );
                my $header        = remove_extension($seqfile);
                $header .= "Coordinates: ";
                $header .= $random_number + 1;
                $header .= "..";
                $header .= $random_number + $opt_l + 1;
                $seq_hash{$header} = $sub_seq;
            }
        }
    }

    return ( \%seq_hash );
}

sub ConcatenateFasta {    # to read fasta format files into hash. joins all the sequences in one file into one sequence and returns joined sequence as string.

    my $seqfile = shift;
    my $remove_amb     = shift;

    my ( $header, @sequence );
    my $sequence2 = ();
    my $sequence  = ();
    chomp $seqfile;
    open FASTA, "$seqfile";

    #	print "reading Sequences from input file.....Plz wait...\n";
    my %seq_hash = ();
    $/ = "\n>";    # Change record seperator to read Fasta

    while (<FASTA>) {
        chomp;
        ( $header, @sequence ) = split( "\n", $_ );

        $header =~ s/>//;       # Remove Leading > from Header
        $header =~ s/\s*$//;    # Remove trailing spaces from header
        $header =~ s/^\s*//;    # Remove Leading spaces from Header

        $sequence = join( "", @sequence );
        $sequence =~ s/\s//g;
        $sequence =~ s/\n//g;
        if ( $remove_amb == 1 ) { $sequence =~ s/[^ATGC]//gi; }    # remove ambigous bases if asked for
        if ( $header =~ /^\s*$/ ) { next; }
        $sequence2 .= $sequence;
    }

    $/ = "\n";                                              # Record seperator set back to default newline.
    $sequence2 =~ s/\s+//g;
    return ($sequence2);

}

# Calculate GC content and return
sub GC_content {
    my $sequence   = shift;
    my $GCnumber   = $sequence =~ tr/GCgc/GCgc/;
    my $length_seq = $sequence =~ tr/ATGCatgc/ATGCatgc/;

    my $GCcontent = $GCnumber / $length_seq;
    return ($GCcontent);
}

sub GC_content_percent {
    my $sequence   = shift;
    my $GCnumber   = $sequence =~ tr/GCgc/GCgc/;
    my $length_seq = $sequence =~ tr/ATGCatgc/ATGCatgc/;

    my $GCcontent = $GCnumber / $length_seq;
    return ( $GCcontent * 100 );
}

sub GC_per_file{
	my$filename=shift;
	my$sequence=ConcatenateFasta($filename,1);
	return(GC_content_percent($sequence));


}

sub min {
    @_ = sort { $a <=> $b } @_;
    return $_[0];

}

sub max {

    @_ = sort { $a <=> $b } @_;
    return $_[-1];

}

sub roundoff {
    my $float          = shift;
    my $place          = shift;
    my $before_decimal = index( $float, "." );
    my $temp_float     = substr( $float, 0, $place + $before_decimal + 2 );
    my $add            = "0." . ( "0" x $place ) . "5";
    $add *= 1;
    $temp_float += $add;
    my $new_float = substr( $temp_float, 0, $place + $before_decimal + 1 );
    $new_float =~ s/\.// if $place == 0;
    return ($new_float);
}

sub log_base {
    my $n    = shift;
    my $base = shift;
    return log($n) / log($base);
}

sub revcomp {
    my $seq = shift;
    $seq = reverse $seq;
    $seq =~ tr/ATGCNatgcn/TACGNtacgn/;
    return $seq;
}

sub getfilename {
    my $rawfilename = shift;

    #remove folder names from rawfilename names.
    my @names = split( /[\/\\]/, $rawfilename );
    return $names[-1];

}

sub remove_extension {
    my $filename = shift;
    my @names = split( /\./, $filename );
    $names[-1] = "";
    return join( "", @names );

}
######################################
sub codon_hash {
    my $codon_hash;
    my %g;
    $g{TCA}{amino} = 'S';
    $g{TCC}{amino} = 'S';
    $g{TCG}{amino} = 'S';
    $g{TCT}{amino} = 'S';
    $g{TTC}{amino} = 'F';
    $g{TTT}{amino} = 'F';
    $g{TTA}{amino} = 'L';
    $g{TTG}{amino} = 'L';
    $g{TAC}{amino} = 'Y';
    $g{TAT}{amino} = 'Y';
    $g{TAA}{amino} = '*';
    $g{TAG}{amino} = '*';
    $g{TGC}{amino} = 'C';
    $g{TGT}{amino} = 'C';
    $g{TGA}{amino} = '*';
    $g{TGG}{amino} = 'W';
    $g{CTA}{amino} = 'L';
    $g{CTC}{amino} = 'L';
    $g{CTG}{amino} = 'L';
    $g{CTT}{amino} = 'L';
    $g{CCA}{amino} = 'P';
    $g{CCC}{amino} = 'P';
    $g{CCG}{amino} = 'P';
    $g{CCT}{amino} = 'P';
    $g{CAC}{amino} = 'H';
    $g{CAT}{amino} = 'H';
    $g{CAA}{amino} = 'Q';
    $g{CAG}{amino} = 'Q';
    $g{CGA}{amino} = 'R';
    $g{CGC}{amino} = 'R';
    $g{CGG}{amino} = 'R';
    $g{CGT}{amino} = 'R';
    $g{ATA}{amino} = 'I';
    $g{ATC}{amino} = 'I';
    $g{ATT}{amino} = 'I';
    $g{ATG}{amino} = 'M';
    $g{ACA}{amino} = 'T';
    $g{ACC}{amino} = 'T';
    $g{ACG}{amino} = 'T';
    $g{ACT}{amino} = 'T';
    $g{AAC}{amino} = 'N';
    $g{AAT}{amino} = 'N';
    $g{AAA}{amino} = 'K';
    $g{AAG}{amino} = 'K';
    $g{AGC}{amino} = 'S';
    $g{AGT}{amino} = 'S';
    $g{AGA}{amino} = 'R';
    $g{AGG}{amino} = 'R';
    $g{GTA}{amino} = 'V';
    $g{GTC}{amino} = 'V';
    $g{GTG}{amino} = 'V';
    $g{GTT}{amino} = 'V';
    $g{GCA}{amino} = 'A';
    $g{GCC}{amino} = 'A';
    $g{GCG}{amino} = 'A';
    $g{GCT}{amino} = 'A';
    $g{GAC}{amino} = 'D';
    $g{GAT}{amino} = 'D';
    $g{GAA}{amino} = 'E';
    $g{GAG}{amino} = 'E';
    $g{GGA}{amino} = 'G';
    $g{GGC}{amino} = 'G';
    $g{GGG}{amino} = 'G';
    $g{GGT}{amino} = 'G';

    # expested frequence
    $g{GCT}{exp_freq} = 0.25;
    $g{GCA}{exp_freq} = 0.25;
    $g{GCC}{exp_freq} = 0.25;
    $g{GCG}{exp_freq} = 0.25;
    $g{AGA}{exp_freq} = 0.17;
    $g{AGG}{exp_freq} = 0.17;
    $g{CGA}{exp_freq} = 0.17;
    $g{CGT}{exp_freq} = 0.17;
    $g{CGG}{exp_freq} = 0.17;
    $g{CGC}{exp_freq} = 0.17;
    $g{AAT}{exp_freq} = 0.5;
    $g{AAC}{exp_freq} = 0.5;
    $g{GAT}{exp_freq} = 0.5;
    $g{GAC}{exp_freq} = 0.5;
    $g{TGT}{exp_freq} = 0.5;
    $g{TGC}{exp_freq} = 0.5;
    $g{TGA}{exp_freq} = 0.33;
    $g{TAA}{exp_freq} = 0.33;
    $g{TAG}{exp_freq} = 0.33;
    $g{CAA}{exp_freq} = 0.5;
    $g{CAG}{exp_freq} = 0.5;
    $g{GAA}{exp_freq} = 0.5;
    $g{GAG}{exp_freq} = 0.5;
    $g{GGA}{exp_freq} = 0.25;
    $g{GGT}{exp_freq} = 0.25;
    $g{GGG}{exp_freq} = 0.25;
    $g{GGC}{exp_freq} = 0.25;
    $g{CAT}{exp_freq} = 0.5;
    $g{CAC}{exp_freq} = 0.5;
    $g{ATT}{exp_freq} = 0.33;
    $g{ATC}{exp_freq} = 0.33;
    $g{ATA}{exp_freq} = 0.33;
    $g{TTG}{exp_freq} = 0.17;
    $g{CTT}{exp_freq} = 0.17;
    $g{CTG}{exp_freq} = 0.17;
    $g{CTC}{exp_freq} = 0.17;
    $g{TTA}{exp_freq} = 0.17;
    $g{CTA}{exp_freq} = 0.17;
    $g{AAG}{exp_freq} = 0.5;
    $g{AAA}{exp_freq} = 0.5;
    $g{ATG}{exp_freq} = 1;
    $g{TTT}{exp_freq} = 0.5;
    $g{TTC}{exp_freq} = 0.5;
    $g{CCA}{exp_freq} = 0.25;
    $g{CCT}{exp_freq} = 0.25;
    $g{CCC}{exp_freq} = 0.25;
    $g{CCG}{exp_freq} = 0.25;
    $g{TCT}{exp_freq} = 0.17;
    $g{TCA}{exp_freq} = 0.17;
    $g{AGT}{exp_freq} = 0.17;
    $g{TCC}{exp_freq} = 0.17;
    $g{AGC}{exp_freq} = 0.17;
    $g{TCG}{exp_freq} = 0.17;
    $g{ACA}{exp_freq} = 0.25;
    $g{ACT}{exp_freq} = 0.25;
    $g{ACC}{exp_freq} = 0.25;
    $g{ACG}{exp_freq} = 0.25;
    $g{TGG}{exp_freq} = 1;
    $g{TAT}{exp_freq} = 0.5;
    $g{TAC}{exp_freq} = 0.5;
    $g{GTT}{exp_freq} = 0.25;
    $g{GTG}{exp_freq} = 0.25;
    $g{GTC}{exp_freq} = 0.25;
    $g{GTA}{exp_freq} = 0.25;

    return ( \%g );

}

# Calculate GC content for each codon in a coding sequence from start to stop codon in given frame.
# need sequence, base# in codon and reading frame.
sub count_GCn {
    my $sequence     = shift;
    my $frame        = shift;
    my $start_at_ATG = shift;
    my $trim_at_stop = shift;
    my $count_this   = shift;
    my $which_base   = shift;

    my $synonymos_codon = 0;
    my $synonymos_GC3   = 0;
    my $orf;
    my $base1 = "G";
    my $base2 = "C";

    if ( uc $count_this eq "AT" ) { $base1 = "A"; $base2 = "T" }

    # Get ORF from the sequence provided. If not asked, use full sequence.
    if   ( lc $start_at_ATG eq 'no' && lc $trim_at_stop eq 'no' ) { $orf = $sequence }
    else                                                          { $orf = get_ORF( $sequence, 1, $start_at_ATG, $trim_at_stop ); }

    $orf = uc $orf;

    # count the number of synonymos codons.
    for ( my $i = $frame - 1; $i <= length($orf) - 3; $i = $i + 3 ) {
        my $codon    = substr( $orf,   $i, 3 );
        my $lastbase = substr( $codon, 2,  1 );
        $synonymos_codon++ if ( $codon ne "ATG" || $codon ne "TGG" );
        $synonymos_GC3++ if ( ( $codon ne "ATG" && $codon ne "TGG" ) && ( $lastbase eq $base1 || $lastbase eq $base2 ) );
    }

    # return fraction of synonymos codons contining G or C at last base.
    return ( $synonymos_GC3 * 100 / $synonymos_codon );
}

########################################################################################################
# Calculate GC content for each codon in a coding sequence from start to stop codon in given frame.
# need sequence, base# in codon and reading frame.
sub count_nt {
    my $sequence     = shift;
    my $frame        = shift;
    my $start_at_ATG = shift;
    my $trim_at_stop = shift;
    my $base         = shift;
    my $which_base   = shift;

    my $synonymos_codon = 0;
    my $synonymos_nt    = 0;
    my $orf;

    $base =~ s/[^\w]//g;

    # Get ORF from the sequence provided. If not asked, use full sequence.
    if   ( lc $start_at_ATG eq 'no' && lc $trim_at_stop eq 'no' ) { $orf = $sequence }
    else                                                          { $orf = get_ORF( $sequence, 1, $start_at_ATG, $trim_at_stop ); }

    $orf  = uc $orf;
    $base = uc $base;

    # count the number of synonymos codons.
    for ( my $i = $frame - 1; $i <= length($orf) - 3; $i = $i + 3 ) {
        my $codon    = substr( $orf,   $i, 3 );
        my $lastbase = substr( $codon, 2,  1 );
        $synonymos_codon++ if ( $codon ne "ATG" || $codon ne "TGG" );
        $synonymos_nt++ if ( ( $codon ne "ATG" && $codon ne "TGG" ) && ( $lastbase eq $base ) );
    }

    # return fraction of synonymos codons contining G or C at last base.
    return ( $synonymos_nt * 100 / $synonymos_codon );
}

# check if a cds is full length starting with ATG and ending with stop codon and there is no stop codon in between.
# usage; get_ORF(sequence,$frame,trim_at_stop(yes|no),start_at_ATG(yes|no))
# returns ORF as specified
sub get_ORF {
    my $sequence     = shift;
    my $frame        = shift;
    my $start_at_ATG = shift;
    my $trim_at_stop = shift;

    $sequence = uc $sequence;
    my $start_codon = 'no';
    my $stop_codons = 0;
    my $ORF;

    # find the first ATG
    my $base = -1;
    if ( $sequence =~ /ATG/i ) {
        do { $base++ } while ( substr( $sequence, $base, 3 ) ne 'ATG' );
    }
    else {
        $base = 0;
    }

    #start ORF from ATG if asked to do.
    $frame = $base + 1 if $start_at_ATG eq 'yes';

    for ( my $i = $frame - 1; $i <= length($sequence); $i = $i + 3 ) {
        my $codon = substr( $sequence, $i, 3 );
        $codon =~ s/\s+//g;

        #find ATG and start joining codons to create ORF
        if ( $codon eq "ATG" && $start_at_ATG eq 'yes' ) { $start_codon = 'yes' }
        elsif ( $start_at_ATG eq 'no' || $base == 0 ) { $start_codon = 'yes' }

        # create ORF
        if ( $start_codon eq 'yes' ) { $ORF .= $codon }

        # stop when meet stop codon.
        last if ( ( $codon eq "TAA" || $codon eq "TAG" || $codon eq "TGA" ) && $trim_at_stop eq 'yes' );
    }

    $ORF =~ s/\s+//g;
    return ($ORF);
}
#####################################################################################
sub collect_positional_codons {

    my $sequence = shift;
    my $frame    = shift;
    my $Ref_hash = shift;

    for ( my $i = $frame - 1; $i * 3 <= length($sequence) - 3; $i++ ) {

        my $codon = substr( $sequence, $i * 3, 3 );
        push( @{ $$Ref_hash{$i} }, $codon );
    }
}

sub collect_positional_NT {

    my $sequence = shift;
    my $Ref_hash = shift;

    for ( my $i = 0; $i <= length($sequence) - 1; $i++ ) {

        my $base = substr( $sequence, $i , 1 );
        push( @{ $$Ref_hash{$i} }, $base );
    }
}


#####################################################################################
sub collect_positional_codons_grouped {

    my $sequence  = shift;
    my $frame     = shift;
    my $Ref_hash  = shift;
    my $low_range = shift;
    my $up_range  = shift;

    my $group_range = $low_range . "--" . $up_range;
    $group_range =~ s/\s+//g;
    for ( my $i = $frame - 1; $i * 3 <= length($sequence) - 3; $i++ ) {

        my $codon = substr( $sequence, $i * 3, 3 );
        push( @{ $$Ref_hash{$group_range}{$i} }, $codon );
    }
}



sub collect_positional_NT_grouped{
	my $sequence  = shift;
    my $Ref_hash  = shift;
    my $low_range = shift;
    my $up_range  = shift;

    my $group_range = $low_range . "--" . $up_range;
    $group_range =~ s/\s+//g;
    for ( my $i = 0; $i <= length($sequence) - 1; $i++ ) {

        my $base = substr( $sequence, $i , 1 );
        push( @{ $$Ref_hash{$group_range}{$i} }, $base );
    }


}

#####################################################################################
sub count_GCn_grouped {

    my $sequence        = shift;
    my $frame           = shift;
    my $Ref_hash        = shift;
    my $low_range       = shift;
    my $up_range        = shift;
    my $temp_count_this = shift;
    my $GCn             = 0;

    my $group_range = $low_range . "--" . $up_range;
    $group_range =~ s/\s+//g;
    $GCn = count_GCn( $sequence, $frame, 'yes', 'yes', 'GC', 3 ) if ( uc $temp_count_this eq "GC3" );
    $GCn = GC_content_percent($sequence) if ( uc $temp_count_this eq "GC" );
    $GCn = count_GCn( $sequence, $frame, 'yes', 'yes', 'AT', 3 ) if ( uc $temp_count_this eq "AT3" );

    push( @{ $$Ref_hash{$group_range} }, $GCn );

}

######################################################################################
# count codons and return hash $$Ref_hash{$codon}{$filename}{$frame}{'count'}
sub count_codon {
    my $filename     = shift;
    my $sequence     = shift;
    my $frame        = shift;
    my $Ref_hash     = shift;
    my $trim_at_stop = shift;
    $sequence = uc $sequence;
    my $stop_codons = 0;

    #print "\n\nSequence:$sequence\nLength:".length($sequence)."\n" if($stop_codons >1);
    for ( my $i = $frame - 1; $i <= length($sequence); $i = $i + 3 ) {

        my $codon = substr( $sequence, $i, 3 );
        $codon =~ s/\s+//g;
        next if $codon =~ m/[^ATGC]/i;
        next if length($codon) < 3;

        if ( $stop_codons == 0 && lc $trim_at_stop eq "yes" ) { $$Ref_hash{$codon}{$filename}{$frame}{'count'}++ }
        elsif ( lc $trim_at_stop eq "no" ) { $$Ref_hash{$codon}{$filename}{$frame}{'count'}++ }

        $stop_codons++ if ( $codon eq "TAA" || $codon eq "TAG" || $codon eq "TGA" );

        #if($stop_codons >1){if ($codon eq "TAA"||$codon eq "TAG"||$codon eq "TGA"){print "*$codon"}else{print "$codon"}}

    }

    if ( lc $trim_at_stop eq "yes" && $stop_codons > 1 ) { print "\nWarning: Sequence has $stop_codons stop codon\nTruncating sequence to first stop codon\n\n" }

}

###########################################################################################
sub average {
    my ($data) = @_;
    if ( not @$data ) {
        die("Empty array\n");
    }
    my $total = 0;
    foreach (@$data) {
        $total += $_;
    }
    my $average = $total / @$data;
    return $average;
}

sub stdev {
    my ($data) = @_;
    if ( @$data == 1 ) {
        return 0;
    }
    my $average = &average($data);
    my $sqtotal = 0;
    foreach (@$data) {
        $sqtotal += ( $average - $_ )**2;
    }
    my $std = ( $sqtotal / ( scalar @$data ) )**0.5;
    return $std;
}
##########################################################################################
