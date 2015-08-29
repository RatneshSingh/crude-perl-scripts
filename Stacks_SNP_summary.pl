#!/usr/bin/perl -w
use strict;
use Getopt::Long;

our ( $AG, $CT, $AC, $AT, $CG, $GT, $transition, $transversion, $total_loci, $unknown, $count, $percent, $abs, $out, $snp_file, $help, $outfh );

my $result = GetOptions(
    "snp|s=s"   => \$snp_file,
    "abs|a"     => \$abs,
    "percent|p" => \$percent,
    "out|o=s"   => \$out,
    "help|h"    => \$help
);


my $usage="

$0 -options
Options:
-s  snp file from stacks
-a  snp count in absolute number[both]
-p  snp count in percent of total snps [both]
-o  output file to save result[STDOUT]
-h  print this help

";

die "$usage" if $help;




if(!$abs && !$percent){$abs=1;$percent=1}
$snp_file=$snp_file?$snp_file:$ARGV[0];

$AG = $CT = $AC = $AT = $CG = $GT = $transition = $transversion = $total_loci = $unknown = $count = 0;
open SNP, "$snp_file";

while (<SNP>) {
    my @line = split( /\s+/, $_ );
    if    ( $line[5] =~ /A|G/ && $line[6] =~ /A|G/ ) { $transition++;   $AG++; $total_loci++ }
    elsif ( $line[5] =~ /C|T/ && $line[6] =~ /C|T/ ) { $transition++;   $CT++; $total_loci++ }
    elsif ( $line[5] =~ /A|C/ && $line[6] =~ /A|C/ ) { $transversion++; $AC++; $total_loci++ }
    elsif ( $line[5] =~ /A|T/ && $line[6] =~ /A|T/ ) { $transversion++; $AT++; $total_loci++ }
    elsif ( $line[5] =~ /C|G/ && $line[6] =~ /C|G/ ) { $transversion++; $CG++; $total_loci++ }
    elsif ( $line[5] =~ /G|T/ && $line[6] =~ /G|T/ ) { $transversion++; $GT++; $total_loci++ }
    else                                             { $unknown++ }
    $count++;
}


if($out){open $outfh, ">$out"}else{$outfh=\*STDOUT}

print $outfh "file_name\tTotal_num_loci\tunprocessed\tTransitions\tTransversions\tAG\tCT\tAC\tAT\tCG\tGT\n";
print $outfh "$snp_file\t$total_loci\t$unknown\t$transition\t$transversion\t$AG\t$CT\t$AC\t$AT\t$CG\t$GT\n" if $abs;
print $outfh round( $total_loci * 100 / $total_loci ) . "%\t"
  . round( $unknown * 100 / $total_loci ) . "%\t"
  . round( $transition * 100 / $total_loci ) . "%\t"
  . round( $transversion * 100 / $total_loci ) . "%\t"
  . round( $AG * 100 / $total_loci ) . "%\t"
  . round( $CT * 100 / $total_loci ) . "%\t"
  . round( $AC * 100 / $total_loci ) . "%\t"
  . round( $AT * 100 / $total_loci ) . "%\t"
  . round( $CG * 100 / $total_loci ) . "%\t"
  . round( $GT * 100 / $total_loci ) . "%\n"
  if $percent;






sub round {
    my $decimals = 2;
    my ($number) = @_;
    substr( $number + ( '0.' . '0' x $decimals . '5' ), 0, $decimals + length( int($number) ) + 1 );
}
