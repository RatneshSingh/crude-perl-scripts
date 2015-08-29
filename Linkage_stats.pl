#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Math::CDF;
use Math::NumberCruncher;
use Math::Random qw(random_poisson);

our ( $file, %lgmarkers, $window, $start, $lg, $slide, %mark_per_window, @all_mark_perwin, $bin_start, $bin_size, $bin_end, $out, $help, %lg_info, $max_clust, $meanperwin,
    $per_lg, $fh, $pabs, $pvalue, $print_raw, $print_minclust );

#### get options
my $result = GetOptions(
    "file|f=s"           => \$file,
    "window|wl=i"        => \$window,
    "initiate|wi=i"      => \$start,
    "slide|ws=i"         => \$slide,
    "bin_start|bi=i"     => \$bin_start,
    "bin_end|be=i"       => \$bin_end,
    "bin_size|bs=i"      => \$bin_size,
    "maxclust|mc=i"      => \$max_clust,
    "meanperwin|mw"      => \$meanperwin,
    "pabs|pa"            => \$pabs,
    "perlg|pl"           => \$per_lg,
    "pvalue_cutoff|pv=f" => \$pvalue,
    "print_raw|pr"       => \$print_raw,
    "print_minclust|pm"  => \$print_minclust,
    "out|o=s"            => \$out,
    "help|h"             => \$help
);

my $usage = "

$0 -options

options:
    file|f        File containing marker names(Col1) and distance in cM (Col2). Linkage groups are inserted in rows between
    window|wl     Length of sliding window [10]
    initiate|wi   Start window at [0]
    slide|ws      Slide window by [1]
    bin_start|bi  start bin for histogram at [0].
    bin_end|be    Stop binning histogram at[100]
    bin_size|bs   Bin size of histogram[1].
    maxclust|mc   exclude if markers/per window is largerthan this.
    meanperwin|mw Calculate mean marker/window
    pabs|pa       Print number of windows per in each group
    print_raw|pr  Print raw per window data in seperate file[linkage_file_name.raw]
    perlg|pl      Calculate stat per linkage group.
    pvalue_cutoff|pv count Num cluster when probability is less than this p-value[0.01],
    out|o         Output file to save results
    help|h        Print help and exit


";

## check for input file or help tag. die if required.
die "$usage" if $help;
open FILE, "$file" or die "\n\nERROR: Input file not found.\n\tPlease provide file of linkage markers\n\n$usage\n";

## prepare to print results
if ($out) { open $fh, ">$out"; }
else      { $fh = \*STDOUT }

### set default values
$bin_start = $bin_start ? $bin_start : 0;
$bin_end   = $bin_end   ? $bin_end   : 100;
$bin_size  = $bin_size  ? $bin_size  : 1;
$window    = $window    ? $window    : 10;
$start     = $start     ? $start     : 0;
$slide     = $slide     ? $slide     : 1;
$max_clust = $max_clust ? $max_clust : 1e10;
$pvalue    = $pvalue    ? $pvalue    : 0.01;



#open RAW, ">${file}.raw.ws$window.wl$slide.table" if $print_raw;
open RAW, ">${file}.raw.table" if $print_raw;

print $fh "\nParameters:: Window_size:$window\tStart_at:$start\tSlide_by:$slide";

### Read the marker information from file having two columns: 1)Marker_name 2) Marker distance in CentiMorghan
while (<FILE>) {
    $lg = $1 if m/^\s*Group\s+LG(\d+)/;
    $lg_info{$lg}{'len'} = $lg_info{$lg}{'len'} ? $lg_info{$lg}{'len'} : 0;
    $lg_info{$lg}{'num'} = $lg_info{$lg}{'num'} ? $lg_info{$lg}{'num'} : 0;
    next if m/Group/i;
    s/^\s+//g;
    my ( $mname, $mdist ) = split(/\s+/);
    $mdist=0.001 if $mdist ==0;
    $mdist =~ s/\s+//g;

    push( @{ $lgmarkers{$lg} }, $mdist ) if $mdist;
    $lg_info{$lg}{'len'} = $mdist if $mdist > $lg_info{$lg}{'len'};
    $lg_info{$lg}{'num'}++;

}

### count the numbers of markers in sliding window
my $all_len = 0;
my $all_num = 0;
my %raw_data;
print $fh "\n#################################################################\n" if $per_lg;
foreach my $lg ( sort { $a <=> $b } keys %lgmarkers ) {

    #print $fh "\nLength of linkage group $lg is: $lg_info{$lg}{'len'}"            if $per_lg;
    #print $fh "\nNumber of markers on linkage group $lg is: $lg_info{$lg}{'num'}" if $per_lg;
    $all_len += $lg_info{$lg}{'len'};    ### calculate total length of linkage groups.
    $all_num += $lg_info{$lg}{'num'};    ### calculate total number markers in linkage groups.

    my @mark_per_window = ();
    for ( $start = 0; $start <= int( max( @{ $lgmarkers{$lg} } ) ); $start += $slide ) {
        my $count = 0;

        #print "\n following centimorghans are in window $start - ",$start+$window,":";
        foreach my $cm ( @{ $lgmarkers{$lg} } ) {
            $count++ if ( $cm >= $start && $cm < $start + $window );

            #print "\t$cm" if ($cm > $start && $cm < $start+$window);
        }

        push( @mark_per_window, $count ) if $count < $max_clust;

        $raw_data{$start}{$lg} = $count;

    }

    push( @all_mark_perwin, @mark_per_window );
    $mark_per_window{$lg} = \@mark_per_window;

    #print "\nNum Markers in Linkage group $lg per window of $window:", join(",",@mark_per_window);
}

### print markers per cM to a file.
if ($print_raw) {
    print RAW "CM";
    foreach my $lnkg ( sort { $a <=> $b } keys %lg_info ) {
        print RAW "\tLG_$lnkg";
    }
    foreach my $start ( sort { $a <=> $b } keys %raw_data ) {
        print RAW "\n$start";
        foreach my $lg ( sort { $a <=> $b } keys %lg_info ) {
            print RAW "\t", $raw_data{$start}{$lg} ? $raw_data{$start}{$lg} : 0;
        }
    }
}

########==============================================================================================
#### stat with per linkage group.
my%min_clust;
foreach my $linkage ( sort { $a <=> $b } keys %mark_per_window ) {

    #foreach(@all_mark_perwin){print "\n$_\t".Math::CDF::ppois($_,6.646789);}
    my $n            = @{ $mark_per_window{$linkage} };
    my $chisq_sum    = 0;
    my $num_clust    = 0;
    my $mark_clust   = 0;
    $min_clust{$linkage}=1e10;
    my $ref_bin_hash = bin( $mark_per_window{$linkage}, $bin_start, $bin_end, $bin_size );
    my $mean_lg      = $lg_info{$linkage}{'num'} * $window / $lg_info{$linkage}{'len'};
    $mean_lg = Math::NumberCruncher::Mean( $mark_per_window{$linkage} ) if $meanperwin;
    print $fh "\n#################################################################\n" if $per_lg;
    print $fh
      "\nTotal Linkage length:$lg_info{$linkage}{'len'}\nTotal number of markers:$lg_info{$linkage}{'num'}\nAverage number of markers in Linkgae Group $linkage per $window cM",
      $meanperwin ? "\(Mean per window\)" : "\(NumMarker/LnkgLen\)", ":$mean_lg"
      if $per_lg;

    print $fh "\n\nLinkageGroup\tMarkers/$window.cM\tExpected_Freq\tObserved_freq" if $per_lg;
    print $fh "\tNumber_of_Windows" if ( $pabs && $per_lg );

    foreach my $clust ( sort { $a <=> $b } keys %$ref_bin_hash ) {

        #my $ppois = $clust < $mean_lg ? Math::CDF::ppois( $clust, $mean_lg ) : 1 - Math::CDF::ppois( $clust, $mean_lg );
        my $exp = exp_ppois( $mean_lg, $clust );

        print $fh "\n$linkage\t$clust\t$exp\t$$ref_bin_hash{$clust}{'rel'}" if $per_lg;
        print $fh "\t$$ref_bin_hash{$clust}{'abs'}" if ( $pabs && $per_lg );

        $mark_clust += ($$ref_bin_hash{$clust}{'abs'}* $clust) if ( $exp <= $pvalue && $clust > $mean_lg );
        $min_clust{$linkage}=$clust if ($min_clust{$linkage} > $clust && $exp <= $pvalue && $clust > $mean_lg);
        #$mark_clust+=($raw_data{$clust}{$linkage}?$raw_data{$clust}{$linkage}:0) if ( $exp <= $pvalue && $clust > $mean_lg );
        $num_clust++ if ( $exp <= $pvalue && $clust > $mean_lg && $$ref_bin_hash{$clust}{'abs'} > 0 );

        #print "\n$clust\t$$ref_rand_bin{$clust}{'rel'}\t$$ref_bin_hash{$clust}{'rel'} ";

        $chisq_sum += ( $$ref_bin_hash{$clust}{'rel'} - $exp )**2 / $exp if $exp > 0;

    }
    print $fh "\n Total number of windows with clusters in LG:$linkage with (p-value <= $pvalue):$num_clust\n"             if $per_lg;
    print $fh "\n Total number of Markers in windows with clusters in LG:$linkage with (p-value <= $pvalue):$mark_clust\n" if $per_lg;

}

if($print_minclust){
open MINCLUST,">$file.MinClustSizePerLG.table";
print MINCLUST "\tMinMark";
foreach my$lnk(sort{$a<=>$b} keys %min_clust){
  print MINCLUST "\n$lnk\t$min_clust{$lnk}"

}
}


#########=============================================================================================
my $ref_bin_hash = bin( \@all_mark_perwin, $bin_start, $bin_end, $bin_size );
my $mean_all = $all_num * $window / $all_len;
$mean_all = Math::NumberCruncher::Mean( \@all_mark_perwin ) if $meanperwin;
my $chisq_sum = 0;
print $fh "\n#################################################################\n";

print $fh "\nTotal Linkage length:$all_len\nTotal number of markers:$all_num\nAverage number of markers per $window cM ",
  $meanperwin ? "\(Mean per window\)" : "\(LgLen/NumMarker\)", ":$mean_all";
print $fh "\n\nLinkageGroup\tMarkers/$window.cM\tExpected_Freq\tObserved_freq";
print $fh "\tNumber_of_Windows" if $pabs;
my $num_clust  = 0;
my $mark_clust = 0;
foreach my $clust ( sort { $a <=> $b } keys %$ref_bin_hash ) {
    my $exp = exp_ppois( $mean_all, $clust );
    print $fh "\nWhole\t$clust\t$exp\t$$ref_bin_hash{$clust}{'rel'}";
    print $fh "\t$$ref_bin_hash{$clust}{'abs'}" if $pabs;

    ### figure out to get number of markers in clusters when window slide is overlapping
    $mark_clust += ($$ref_bin_hash{$clust}{'abs'} * $clust) if ( $exp <= $pvalue && $clust > $mean_all );
    $num_clust++ if ( $exp <= $pvalue && $clust > $mean_all && $$ref_bin_hash{$clust}{'abs'} > 0 );

    $chisq_sum += ( $$ref_bin_hash{$clust}{'rel'} - $exp )**2 / $exp if $exp > 0;
}
print $fh "\n Total number of windows with clusters in all linkage groups with (p-value <= $pvalue):$num_clust\n";
print $fh "\n Total number of markers in windows with clusters in all linkage groups with (p-value <= $pvalue):$mark_clust\nPS: Overlapping windows will show inflated marker numbers. For accurate number of markers, use non-overlapping window slide.";

print "\n\n\n";

#open OUT, ">$out";
#print OUT join( "\n", @all_mark_perwin );

######
### subroutine
sub min {
    @_ = sort { $a <=> $b } @_;
    return $_[0];

}

sub max {

    @_ = sort { $a <=> $b } @_;
    return $_[-1];
}

sub bin {
    my $ref_data_array = shift;
    my $bin_start      = shift;
    my $bin_end        = shift;
    my $bin_size       = shift;

    my %bin_hash;
    my $num_data = scalar(@$ref_data_array);
    $bin_start = $bin_start ? $bin_start : 0;
    $bin_end   = $bin_end   ? $bin_end   : 100;
    $bin_size  = $bin_size  ? $bin_size  : 10;

    #print "\nNum_Data:$num_data\tBin_start:$bin_start\tBin_end:$bin_end\tBin_size:$bin_size";

    my %bin;
    ## Initialize all the bins windows to 0
    for ( my $i = $bin_start; $i <= $bin_end; $i += $bin_size ) { $bin{$i} = 0 }
    foreach my $value (@$ref_data_array) {

        #print "\ntesting value:$value";
        for ( my $i = $bin_start; $i <= $bin_end; $i += $bin_size ) {

            if ( $value >= $i && $value < $i + $bin_size ) {
                $bin{$i}++;

                #print "\n The value $value is in bin:$i"
            }
            elsif ( $value > $bin_end ) {
                $bin{$bin_end}++;
            }

        }
    }

    foreach my $freq ( keys %bin ) {

        $bin_hash{$freq}{'abs'} = $bin{$freq};
        $bin_hash{$freq}{'rel'} = $bin{$freq} / @$ref_data_array;

        #print "$freq";

    }

    return ( \%bin_hash );

}

sub mean {
    return sum(@_) / @_;
}

sub exp_ppois {
    my $mean  = shift;
    my $value = shift;

    my $minus_mean = "-$mean";

    my $P = Math::CDF::ppois( $value, $mean ) - ( $value - $bin_size > 0 ? ( Math::CDF::ppois( $value - $bin_size >= 0 ? $value - $bin_size : 0, $mean ) ) : 0 ) if $bin_size > 1;

    #$P=Math::CDF::ppois( $value, $mean );
    $P = ( ( 2.71828**$minus_mean ) * ( $mean**$value ) ) / factorial($value) if $bin_size == 1;

    return ($P);
}

sub factorial {
    my $n    = shift;
    my $prod = 1;
    $prod *= $n-- while $n > 0;
    return $prod;
}
