#!/usr/bin/perl
# This script parses the regular blastoutput file and converts into table format.
# It will also calculate the query coverage and subject coverage based on the query and
# subject length. coverage clculations are not possible from table format blast result.
use strict;
use warnings;
use Bio::SearchIO;
use Getopt::Long;

# Usage to be printed out for user.
my $usage = "

This script reads standard blast result and convert it into blast table format.
The blast results can be filtered using several parameter.

Usage: Perl Script -i blast_result_file [options]

-l	Alignment length to filter. Results above this threshold will be reported [All]
-p	Percent value to filter. Results above this threshold will be reported [All].
-e	Evalue value to filter. Results below this threshold will be reported [All].
-c	Query Coverage value to filter. Results above this threshold will be reported.
	Coverage is calculated for each hsp individually [All].
-d	Description in hit. Hits matching description will be reported [All].
-o	Output file to save result [STDOUT]
-s	Remove selfhits [FALSE];
-b  subject|query.Pick hits at end of subject or query.
-m Max distance allowed for the hit from the end.
";

our ( $format, $infile, $outfile, $evalue, $coverage, $percent, $description, $alnlength, $OUT, $opt_selfhit, $bac_end, $max_dist ) = ('blast');
GetOptions(
    'f|format:s'      => \$format,
    'i|infile:s'      => \$infile,
    'o|outfile:s'     => \$outfile,
    'e|evalue:f'      => \$evalue,
    'c|coverage:f'    => \$coverage,
    'p|percent:f'     => \$percent,
    'd|description:s' => \$description,
    'l|alnlength:i'   => \$alnlength,
    's|selfhits'      => \$opt_selfhit,
    'b|be|bac_end=s'  => \$bac_end,
    'max_dist|md|m=i' => \$max_dist
);

# Pick the file name from commandline if option flags are are not used.
if (@ARGV) { $infile = shift; }

# print Usage instructions when blast file name is not provided.
print "$usage" if !defined $infile;

my $searchIO = Bio::SearchIO->new( -format => $format, -file => $infile );

# Open file to save output or direct to STDOUT when not defined.
if ($outfile) { open( $OUT, '>', $outfile ) || die "Unable to open $outfile for writing\n $usage"; }
else          { $OUT = \*STDOUT; }

#Set parameters to default values for filtering when not defined.
if ( !$coverage )    { $coverage    = 0; }
if ( !$percent )     { $percent     = 0; }
if ( !$evalue )      { $evalue      = 10; }
if ( !$description ) { $description = ""; }
if ( !$alnlength )   { $alnlength   = 1; }

# Process blast results.
while ( my $result = $searchIO->next_result ) {
    while ( my $hit = $result->next_hit ) {
        while ( my $hsp = $hit->next_hsp ) {

            # calculate mismatches. Not provided by SearchIO methods.
            my $mismatch = $hsp->length('total') - $hsp->num_conserved - $hsp->gaps;

            # Calculate query coverage using Alignment length. it includes gaps.
            #my$current_coverage=($hsp->length('total')/$result->query_length)*100;

            # Calculate query coverage using Q start and Q end. Does not includes gaps.
            my $current_coverage = ( $hsp->length('query') / $result->query_length ) * 100;
            my $subject_coverage = ( $hsp->length('hit') / $hit->length ) * 100;

            #print $hsp->query->strand < 0 ? ( $hsp->query->start - $hsp->query->end ):( $hsp->query->end - $hsp->query->start )."\n";

            next if $result->query_name() eq $hit->name() && $opt_selfhit;

				next
              if ( $bac_end && $bac_end eq 'subject' && ( $max_dist < ($hit->length - max( $hsp->hit->end, $hsp->hit->start )) && $max_dist < min( $hsp->hit->end, $hsp->hit->start )) )
              || ( $bac_end && $bac_end eq 'query' &&   ($max_dist  < ($result->query_length - max( $hsp->query->end, $hsp->query->start      )) && $max_dist < min( $hsp->query->end, $hsp->query->start)));

            #print join ("\t",$hsp->length('total'),$result->query_length, $current_coverage, "\n" );

            if (     $hsp->hsp_length >= $alnlength
                  && $hsp->percent_identity >= $percent
                  && $hsp->evalue <= $evalue
                  && $current_coverage >= $coverage
                  && $hit->description =~ m/$description/i )
              {

                  #print $hit->description."\n";
                  print $OUT join(
                      "\t",
                      (
                          $result->query_name,
                          $hit->name,
                          sprintf( "%.2f", $hsp->percent_identity ),
                          $hsp->length('total'),
                          $mismatch,
                          $hsp->gaps('total'),

                          # flip start/end on rev strand
                          $hsp->query->strand < 0 ? ( $hsp->query->end, $hsp->query->start ) : ( $hsp->query->start, $hsp->query->end ),
                          $hsp->hit->strand < 0   ? ( $hsp->hit->end,   $hsp->hit->start )   : ( $hsp->hit->start,   $hsp->hit->end ),
                          $hsp->evalue,
                          $hsp->bits,
                          sprintf( "%.2f", $current_coverage ),
                          sprintf( "%.2f", $subject_coverage ),
                          $result->query_length,
                          $hit->length
                      )
                    ),
                    "\n";

            }
        }
    }
}




sub min{
  @_= sort { $a <=> $b }@_;
   return $_[0];

}

sub max{

 @_ = sort { $a <=> $b }@_;
   return $_[-1];


}
