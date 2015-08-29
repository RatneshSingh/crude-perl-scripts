use warnings;
use strict;
use Getopt::Long;
use Parallel::ForkManager;

our ( $aln, $out, $miss, $rlt, $rlc, $help, $MAX_PROCESSES );

my $options = GetOptions(
    "alignment|a=s"    => \$aln,
    "out|o=s"          => \$out,
    "missing|m=s"      => \$miss,
    "rel_2_total|rlt"  => \$rlt,
    "rel_2_common|rlc" => \$rlc,
    "cpu|c=i"          => \$MAX_PROCESSES,
    "help|h"           => \$help

);

my $usage = "

perl $0 -options...

-alignment|a  Alignment file name.
-out|o        Outputfilename to save results.
-missing|m    Missing values are represented as. [N]
-rlt          Calculate percents against total number of sites[default]
-rlc          Calculate percents against total common sites used (Total - missing)
-cpu|c        Max number of processes to run in parallel.
-help|h       Print help and exit.
";

die "\n\n$usage\n\n" if $help;


###### Assign default values for some parameters
if(!$out){$out = $rlc?"$aln.alnstat.rel2usable.table":"$aln.alnstat.rel2total.table";}
$rlt           = $rlc           ? undef          : 1;
$miss          = $miss          ? $miss          : "N";
$miss =~ s/\s+//g;
$MAX_PROCESSES = $MAX_PROCESSES ? $MAX_PROCESSES : 10;
my $tot = $rlc ? "Usable_sites" : "Total_sites";

##### Open files for reading alignment and writing results
open FASTA, "$aln"  or die "Cannot open alignment file $aln\n";
open OUT,   ">$out" or die "cannot open output file for writing.\n";

####### Read sequences from file and add to hash.
my %seq;
$/ = "\n>";
while (<FASTA>) {
    chomp;

    my ( $header, @sequence ) = split( /\n/, $_ );
    $header =~ s/>//;
    $header = break_at_space($header);

    my $sequence = join( "", @sequence );
    $sequence =~ s/\s+|>//g;

    my @bases = split( //, $sequence );

    $seq{$header} = \@bases;
}
close(FASTA);
$/ = "\n";


###### Create pairs from hash and compute statistics for individual sequences.
print OUT"\nSeq\tTotal_sites\tHomozygous_Sites/$tot\tHeterozygous_Sites/$tot\tMissing_data/Total_Sites";
foreach my $header ( keys %seq ) {
    print OUT"\n$header\t", scalar @{ $seq{$header} };
    print OUT"\t",          individual_stat( \@{ $seq{$header} } );

}

##### Pairwise stat calculations
my $pm = Parallel::ForkManager->new($MAX_PROCESSES);  ### to create parallel runs of subroutines.

print OUT
  "\n\nSequence_1\tSequence_2\tHomozygous_Sites/$tot\tHomozygous_sites_with_ATGC/$tot\tIdentical_but_Ambiguous/$tot\tHeterozygous_Sites/$tot\theterozygous_with_ATGC/$tot\theterozygous_both_Ambiguous/$tot\theterozygous_one_Ambiguous/$tot\tSites_with_Missing_data/Total_sites\tSites_with_usable_data/Total_sites";
foreach my $header ( keys %seq ) {
    foreach my $head2 ( keys %seq ) {

        ###### Run pairwise stat steps in parallel using Parallel::ForkManager
        my $pid = $pm->start and next; ### Forks and returns the pid for the child:
        print OUT"\n$header\t$head2\t", pairwise_stat( $seq{$header}, $seq{$head2} );
        $pm->finish;    ### Terminates the child process

    }

    delete $seq{$header};

}
$pm->wait_all_children;
print "\n\n";

##############################################

sub individual_stat {

    my $array_ref = shift;
    my $hmz_num   = 0;
    my $het_num   = 0;
    my $miss_num  = 0;
    my $tot_num   = scalar @{$array_ref};
    foreach ( @{$array_ref} ) {
        $hmz_num++  if $_ =~ /[ATGC]/i;
        $het_num++  if $_ =~ /[KMRYSWBVHD]/i;
        $miss_num++ if $_ =~ /$miss/i;
    }
    my$usable_sites=$tot_num;
    $usable_sites=$tot_num-$miss_num if $rlc;
#return "hmz_num:$hmz_num(", ( sprintf "%.2f", ( $hmz_num * 100 / $tot_num ) ), "%)\thet_num:$het_num(", ( sprintf "%.2f", ( $het_num * 100 / $tot_num ) ),"%)\tmiss_num:$miss_num(", ( sprintf "%.2f", ( $miss_num * 100 / $tot_num ) ), "%)";
    return "$hmz_num/$usable_sites(", ( sprintf "%.2f", ( $hmz_num * 100 / $usable_sites ) ), "%)\t$het_num/$usable_sites(", ( sprintf "%.2f", ( $het_num * 100 / $usable_sites ) ), "%)\t$miss_num/$tot_num(",
      ( sprintf "%.2f", ( $miss_num * 100 / $tot_num ) ), "%)";
}

sub pairwise_stat {
    my $seq1     = shift;
    my $seq2     = shift;
    my $seq1_ind = $#{$seq1};
    my $seq2_ind = $#{$seq2};

    return 0 if $seq1_ind != $seq2_ind;

    my $homo            = 0;
    my $hetero          = 0;
    my $missing         = 0;
    my $homo_ATGC       = 0;
    my $homo_Ambi       = 0;
    my $hetero_ATGC     = 0;
    my $hetero_AmbiBoth = 0;
    my $hetero_AmbiOne  = 0;
    for ( 0 .. $seq1_ind ) {
        #### Count the number of homozygous sites (Both sequence share same bases at site.)
        if ( $seq1->[$_] eq $seq2->[$_] && $seq1->[$_] !~ /$miss/i && $seq2->[$_] !~ /$miss/i ) {
            if    ( $seq1->[$_] =~ /[ATGC]/i       && $seq2->[$_] =~ /[ATGC]/i )       { $homo_ATGC++ }
            elsif ( $seq1->[$_] =~ /[KMRYSWBVHD]/i && $seq2->[$_] =~ /[KMRYSWBVHD]/i ) { $homo_Ambi++ }
            $homo++

        }
        #### Count the number of heterozygous sites (Both sequence share same bases at site.)
        if ( $seq1->[$_] ne $seq2->[$_] && $seq1->[$_] !~ /$miss/i && $seq2->[$_] !~ /$miss/i ) {
            if    ( $seq1->[$_] =~ /[ATGC]/i       && $seq2->[$_] =~ /[ATGC]/i )       { $hetero_ATGC++ }
            elsif ( $seq1->[$_] =~ /[KMRYSWBVHD]/i && $seq2->[$_] =~ /[KMRYSWBVHD]/i ) { $hetero_AmbiBoth++ }
            elsif ( $seq1->[$_] =~ /[KMRYSWBVHD]/i || $seq2->[$_] =~ /[KMRYSWBVHD]/i ) { $hetero_AmbiOne++ }
            $hetero++;
        }
        $missing++ if ( $seq1->[$_] eq "N" || $seq2->[$_] eq "N" );
    }

    my $divider = $rlc ? $seq1_ind - $missing + 1 : $seq1_ind + 1;

# return "Homozygous sites:$homo(", ( sprintf "%.2f", ( $homo * 100 / ($seq1_ind+1) ) ), "%)\thet_num:$hetero(", ( sprintf "%.2f", ( $hetero * 100 / ($seq1_ind+1) ) ),"%)\tmiss_num:$missing(", ( sprintf "%.2f", ( $missing * 100 / ($seq1_ind+1) ) ), "%)";
    return "$homo/$divider(", ( sprintf "%.2f", ( $homo * 100 / ($divider) ) ),
      "%)\t$homo_ATGC/$divider(",       ( sprintf "%.2f", ( $homo_ATGC * 100 /       ($divider) ) ),
      "%)\t$homo_Ambi/$divider(",       ( sprintf "%.2f", ( $homo_Ambi * 100 /       ($divider) ) ),
      "%)\t$hetero/$divider(",          ( sprintf "%.2f", ( $hetero * 100 /          ($divider) ) ),
      "%)\t$hetero_ATGC/$divider(",     ( sprintf "%.2f", ( $hetero_ATGC * 100 /     ($divider) ) ),
      "%)\t$hetero_AmbiBoth/$divider(", ( sprintf "%.2f", ( $hetero_AmbiBoth * 100 / ($divider) ) ),
      "%)\t$hetero_AmbiOne/$divider(",  ( sprintf "%.2f", ( $hetero_AmbiOne * 100 /  ($divider) ) ),
      "%)\t$missing/", $seq1_ind + 1, "(", ( sprintf "%.2f", ( $missing * 100 / ( $seq1_ind + 1 ) ) ),
      "%)\t", $seq1_ind - $missing + 1, "/", $seq1_ind + 1, "(", ( sprintf "%.2f", ( ( $seq1_ind - $missing + 1 ) * 100 / ( $seq1_ind + 1 ) ) ), "%)";

}

sub break_at_space {
    my $string = shift;
    $string =~ s/\s+.*//;

    return ($string);
}
