#/usr/bin/perl -w
use strict;
use Bio::Tools::GFF;
use Getopt::Long;

our($range_start,$range_end,$gff_file,$ver,$help);

my$result=GetOptions(
  "gff=s"=>\$gff_file,
  "ver=i"=>\$ver,
  "start=i"=>\$range_start,
  "end=i"=>\$range_end,
  "help"=>\$help
  );


my$help="

perl $0 -options...
-gff  gff file path.
-ver  gff version[3]
-start  start site on sequence.
-end  end site on the sequence.
-help print help.
";

print $help if (!$gff_file || !$range_start || !$range_end);
$ver=$ver?$ver:3;
my @features_in_range = ( );


my $gffio = Bio::Tools::GFF->new(-file => $gff_file, -gff_version => $ver);

while (my $feature = $gffio->next_feature()) {

    ## What about features that are not contained within the coordinate range but
    ## do overlap it?  Such features won't be caught by this check.
    if (
        ($feature->start() >= $range_start)
        &&
        ($feature->end()   <= $range_end)
       ) {

        push @features_in_range, $feature;

    }

}

my $gff_out = Bio::Tools::GFF->new(-file => ">$gff_file.selected.gff3", -gff_version => $ver);
foreach (@features_in_range){
  $gff_out->write_feature($_);
}





$gffio->close();
