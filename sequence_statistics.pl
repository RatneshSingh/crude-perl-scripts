#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use GD::Graph::lines;
use Data::Dumper;
our ( $opt_s,$opt_sample,@files );
#getopt('s');
$opt_sample=1; # print each data value in OUT3
my$results=GetOptions( "sequence|s=s"=>\$opt_s,
					  "sample=i"=>\$opt_sample
);

my$usage="

perl script -options sequence_files........

options:
-sample	N.	Sample every Nth value while calculating cumulative distribution.
			Output files end having very large number of cumulative values
			if sequence have very large number of contigs. Chossing larger samples
			number will reduce the number of values in output file and still give similar
			graph. By default, 1st,3rd, 5th, 10th, 20th, 50th,75th and 100th value will be written
			no matter what sample value you choose.[1]



";


##############################################################################################################################################################
# get the names of files from command line
if ($opt_s) { push( @files, $opt_s ); }
if ( $ARGV[0] ) {
    foreach my $file (@ARGV) { $file =~ s/\s+//g; push( @files, $file ) if $file !~ /^\s*$/; }
}

die "No Sequence file found" if scalar @files < 1;

# clean the file name  and create output file name
my $filename_noext = remove_extension( $files[0] );
my $outfile;
if   ( scalar @files > 1 ) { $outfile = "SeqStat_" . $filename_noext . ".Multiple.SeqStat.table" }
else                       { $outfile = "SeqStat_" . $filename_noext . ".table" }
open OUT,  ">$outfile";
open OUT2, ">$outfile.Ndist.table";
open OUT3, ">$outfile.CumulativeDistribution.table";
################################################################################################################################################################
# Calculate and print N75, N50, and N25
print OUT2 "\t";
for ( my $i = 1; $i <= 100; $i++ ) { print OUT2 "$i\t" }

print OUT
  "\tTotal_number_of_sequences\tTotal_number_of_bases\tAverage_Length\tN50\tN25\tN75\tlargest_seq\tsmallest_seq\tseqs_gt_1000\tseq_gt_5000\tseqs_gt_10000\n";

my @final_Nlens_arrays;    # to collect N length values for graph generation
my %cumulative_values;     # to collect the cumulative lengths from largest to smallest sequence
my $cumulative_max = 0;    # To collect the maximum number of sequence among all the file.
my @list_filenames;        # to collect file names to keep them in order.
foreach my $file (@files) {
    my $filename = remove_extension($file);
    open FASTA, "$file";
    $/ = "\n>";
    my $count = 0;
    my @seq_lengths;
    my @rev_seq_lengths;
    my $total_length   = 0;
    my $count_gt_1000  = 0;
    my $count_gt_5000  = 0;
    my $count_gt_10000 = 0;

    # read fasta file and collect the sequence lengths
    while (<FASTA>) {
        $count++;
        ( my $header, my @sequence ) = split( /\n/, $_ );
        my $sequence = join( "", @sequence );
        $sequence =~ s/\s+//g;
        my $length = length $sequence;
        push( @seq_lengths, $length );
        $count_gt_1000++  if $length >= 1000;
        $count_gt_5000++  if $length >= 5000;
        $count_gt_10000++ if $length >= 10000;
        $total_length += $length;
    }

    my $average = int( $total_length / $count );

    # sort the array and create reverse sorted array
    @seq_lengths     = sort { $a <=> $b } @seq_lengths;
    @rev_seq_lengths = sort { $b <=> $a } @seq_lengths;

    # calculate distribution of N
    my $Ref_array_Nlens = N_distribution( \@rev_seq_lengths );
    push( @final_Nlens_arrays, \@{$Ref_array_Nlens} );

    my $N50 = $$Ref_array_Nlens[50];
    my $N25 = $$Ref_array_Nlens[25];
    my $N75 = $$Ref_array_Nlens[75];

#print "Total number of sequences:$count \nTotal number of bases:$total_length \nAverage Length:$average\nN50=$seq_lengths[int($count/2)]\nN25=$seq_lengths[int($count/4)]\nN75=$seq_lengths[int($count*3/4)] \n\n";
    print OUT getfilename($file),
      "\t$count\t$total_length\t$average\t$N50\t$N25\t$N75\t$seq_lengths[-1]\t$seq_lengths[0]\t$count_gt_1000\t$count_gt_5000\t$count_gt_10000\n";

    # Draw Cumulative graph
    my $cumulative = 0;
    $cumulative_values{$filename}{'count'} = scalar @rev_seq_lengths;
    $cumulative_max = $count if $cumulative_max < $count;
    for ( my $i = 0; $i < scalar @rev_seq_lengths; $i++ ) {

        $cumulative += $rev_seq_lengths[$i];
        push( @{ $cumulative_values{$filename}{'values'} }, $cumulative );
    }
    push( @list_filenames, $filename );

}

####################################################################################################################################
# Print cumulative data in a file to draw graph by using excel.

foreach my $file (@list_filenames) { print OUT3 "\t$file"; }
for ( my $i = 0; $i < $cumulative_max; $i++ ) {

	if($i%$opt_sample==0||$i==0||$i==2||$i==4||$i==9||$i==19||$i==49||$i==74||$i==99){
		print OUT3 "\n",$i+1;
		foreach my $file (@list_filenames) {
			print OUT3 "\t";
			if ( defined ${ $cumulative_values{$file}{'values'} }[$i] ) {
				print OUT3 ${ $cumulative_values{$file}{'values'} }[$i];
			}
		}
	}
}
######################################################################################################################################
#draw_line_chart(\@final_Nlens_arrays );

my @data_list;

# undef the values till the largest count so that the number of data points are equal in all the arrays. If not, Graph module will throw an error..
foreach my $filenam (%cumulative_values) {
    for ( my $i = $cumulative_values{$filenam}{'count'}; $i < $cumulative_max; $i++ ) {
        ${ $cumulative_values{$filenam}{'values'} }[$i] = undef;
    }
    push( @data_list, \@{ $cumulative_values{$filenam}{'values'} } );

}

draw_line_chart( \@data_list );

#    print Dumper @final_Nlens_arrays;
#	my $graph = GD::Graph::lines->new(400, 300);
#	my $gd = $graph->plot(\@final_Nlens_arrays) or die $graph->error;
#	open(IMG, '>file.png') or die $!;
#	binmode IMG;
#	print IMG $gd->png;

###############################################################
sub remove_extension {

    my $full_path = shift;
    my $filename  = getfilename($full_path);
    my @names     = split( /\./, $filename );
    $names[-1] = "";
    return join( "", @names );
}

sub getfilename {
    my $rawfilename = shift;

    #print "get file name:$rawfilename\n";

    #remove folder names from rawfilename names.
    my @names = split( /[\/\\]/, $rawfilename );
    return $names[-1];

}

#####################################################################################
sub N_distribution {

    my $Ref_array_seq_length = shift;

    # Arrange the length in reverse order (largest first)
    my @rev_seq_lengths = sort { $b <=> $a } @{$Ref_array_seq_length};

    my $total_length = sum_array( \@rev_seq_lengths );
    my $count        = scalar @rev_seq_lengths;
    my $average      = int( $total_length / $count );

    # calculate the Sequence length for each N values. Total will be 100 such values
    my %dist_limits;

    #my %cumulative_lengths;
    for ( my $i = 0; $i <= 100; $i++ ) {
        $dist_limits{$i}{'cum_limit'} = int( $total_length * $i / 100 );

        #print "\ndist_limits:$i:$dist_limits{$i}{'cum_limit'}\n";
    }

    my $cumulative_len = 0;
    foreach my $seq_len (@rev_seq_lengths) {

        $cumulative_len += $seq_len;

        #print "\n cumulative_len:$cumulative_len\n";
        for ( my $i = 100; $i >= 1; $i-- ) {
            $dist_limits{$i}{'N_len'} = $seq_len
              if ( $cumulative_len <= $dist_limits{$i}{'cum_limit'}
                && $cumulative_len > $dist_limits{ $i - 1 }{'cum_limit'} );
        }

    }

    print OUT2 "\nN\t";
    for ( my $i = 1; $i <= 100; $i++ ) { print OUT2 "$dist_limits{$i}{'N_len'}\t" if defined $dist_limits{$i}{'N_len'} }

    my $Nlens = [];
    for ( my $i = 0; $i <= 100; $i++ ) {
        if   ( defined $dist_limits{$i}{'N_len'} ) { $$Nlens[$i] = $dist_limits{$i}{'N_len'} }
        else                                       { $$Nlens[$i] = undef }

    }
    return $Nlens;
}

##############################################################################################

sub sum_array {
    my $Ref_array_nums = shift;
    my $sum            = 0;
    foreach ( @{$Ref_array_nums} ) { $sum += $_; }
    return $sum;
}

#############################################################################################
sub draw_line_chart {
    my $Ref_array_data = shift;

    # set default values
    my $line_width    = 2;
    my $x_label       = 'X-axis';
    my $y_label       = 'Y-axis';
    my $title         = 'Line Graph';
    my $y_max_value   = 100;
    my $y_tick_number = 10;
    my $y_label_skip  = 2;

    # change the values if provided by the user.
    #$line_width     = shift;
    #$x_label        = shift;
    #$y_label        = shift;
    #$title          = shift;
    #$y_max_value    = shift;
    #$y_tick_number  = shift;
    #$y_label_skip   = shift;

    my $graph = GD::Graph::lines->new( 400, 300 );
    $graph->set(
        line_width => $line_width,
        x_label    => $x_label,
        y_label    => $y_label,

        # title         => $title,
        # y_max_value   => $y_max_value,
        #y_tick_number => $y_tick_number,
        #y_label_skip  => $y_label_skip
    );
    my $gd = $graph->plot($Ref_array_data) or die $graph->error;
    open( IMG, '>file.png' ) or die $!;
    binmode IMG;
    print IMG $gd->png;
}
