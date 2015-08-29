#!/usr/bin/perl -w
####
use strict;
use Data::Dumper;
use Getopt::Long;

our($seq_file,$re_site,$help,$total_count,$perseq,$minfrag,$maxfrag,$bin,$histo);

GetOptions(
           'seq|s=s' => \$seq_file,
           're_site|r=s' => \$re_site,
           'perseq|p'=>\$perseq,
		   'maxfrag|x=i'=>\$maxfrag,
		   'minfrag|m=i'=>\$minfrag,
		   'bin|b=i'=>\$bin,
		   "histo|hg"=>\$histo,
           'help|h'=>\$help
           );

my$usage="

Usage:\nperl $0 -s seqfilename  -r RES_SITE

:options:
-s  Sequence file with sequences in fasta format.
-r  Restriction/or site to look for.
-p  print Number of restriction sites found per seq in file.
-hg Print distribution info for restriction fragments created by res site.
-m	Min limit for histogram data[0]
-x	Max limit for histogram data[5000]
-b	Bin size for histogram data[10].
-h  print help.


";

$minfrag=$minfrag?$minfrag:0;
$maxfrag=$maxfrag?$maxfrag:5000;
$bin=$bin?$bin:100;


unless ( $seq_file && $re_site ) {
    print "\n\n$usage\n\n";
    exit();
}

die "\n\n$usage\n" if $help;

#Get the reverse complement of the RE site - for RE sites that are not palindromes


my $re_site_r = reverse($re_site);
$re_site_r =~ tr/ACGTacgt/TGCAtgca/;


#clean sequence file for white spaces.
open FASTA, "$seq_file";
$total_count=0;  ## to store total counts of restriction sites.
my @frags; ## store restriction fragments.

$/="\n>";
while (<FASTA>) {
    chomp;
    my ( $header, @sequence ) = split( /\n/, $_ );
    $header =~ s/>//;

    #if ( defined $opt_g ) { my @names = split( /\s/,       $header ); $header = $names[0]; }
    #if ( defined $opt_d ) { my @names = split( /\Q$opt_d/, $header ); $header = $names[ $opt_c - 1 ]; }
    $header=~s/^\s+//g;
    $header=~s/\s+$//g;

    my $sequence = join( "", @sequence );
    $sequence =~ s/\s+//g;
    $sequence =~ s/\n+//g;

    #$seq{$header}{'sequence'} = $sequence;
    #$seq{$header}{'len'}      = length($sequence);

    #$count = () = $sequence =~ m/$re_site|$re_site_r/gi;
    my (@starts,@ends);
    my$count=0;

    while ( $sequence =~ /$re_site|$re_site_r/gi ) {

       $starts[$count]=$-[0];
       $ends[$count]=$+[0];

      my$fragment_size=$ends[$count] - $starts[$count-1] if $count>0;

      if (defined $fragment_size && $fragment_size > 1){
        push(@frags,$fragment_size) if (defined $fragment_size && $fragment_size > 1);
        #$maxfrag=$fragment_size if (!$maxfrag && $maxfrag <= $fragment_size);
        #$minfrag=$fragment_size if (!$minfrag && $minfrag >= $fragment_size);
      }
       $count++;
       print"\n$header:$count $re_site" if $perseq;

   }

	$total_count+=$count;


}
#$bin=roundup(($maxfrag-$minfrag)/100) if !$bin;


histogram(\@frags,$minfrag,$maxfrag,$bin) if $histo;


print "\n\nTotal number of $re_site : $total_count\n\n";

exit();

sub histogram{
  my$refDataArray=shift;
  my$lowerlim=shift;
  my$upperlim=shift;
  my$bin=shift;
  my %percent_freq;

  open OUT, ">HistogramData.table";


	foreach(@$refDataArray){
		for(my$i=$lowerlim;$i<=$upperlim;$i=$i+$bin){ if ($_ > $i && $_ <= ($i+$bin)){$percent_freq{$i+$bin}++;}}
		if ($_ > $upperlim){$percent_freq{$upperlim+$bin+}++;}

	}




for(my$i=$lowerlim;$i<=$upperlim;$i=$i+$bin){
		print OUT "\t$i-",$i+$bin;
		#print "\t$percent_freq{$i+$bin}" if defined $percent_freq{$i+$bin};
	}
print OUT "\nFrequency";

for(my$i=$lowerlim;$i<=$upperlim;$i=$i+$bin){
		print OUT  "\t$percent_freq{$i+$bin}" if defined $percent_freq{$i+$bin};
		#print "\t$percent_freq{$i+$bin}" if defined $percent_freq{$i+$bin};
	}
}


sub roundup {
    my $n = shift;
    return(($n == int($n)) ? $n : int($n + 1))
}
