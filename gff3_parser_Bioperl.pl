#/usr/bin/perl -w
use strict;
use Bio::Tools::GFF;
use Bio::SeqIO;
use Getopt::Long;

our($range_start,$range_end,$gff_file,$ver,$help,$fasta);

my$result=GetOptions(
  "gff=s"=>\$gff_file,
  "fasta=s"=>\$fasta,
  "ver=i"=>\$ver,
  "start=i"=>\$range_start,
  "end=i"=>\$range_end,
  "help"=>\$help
  );


$help="

perl $0 -options...
-gff  gff file path.
-fasta fasta file of sequences
-ver  gff version[3]
-start  start site on sequence.
-end  end site on the sequence.
-help print help.
";

print $help if (!$gff_file || !$range_start || !$range_end);
$ver=$ver?$ver:3;
my @features_in_range = ( );


my $gffio = Bio::Tools::GFF->new(-file => $gff_file, -gff_version => $ver);
my$seq= Bio::SeqIO->new(-file=> $fasta, -format=>'fasta');
my$count=0;


while (my $feat = $gffio->next_feature()) {
  if ($feat->primary_tag eq 'cds') {
    print join " ",$feat->start,$feat->end;
  }



        $count++;
}






$gffio->close();
