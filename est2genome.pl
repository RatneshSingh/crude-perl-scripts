#!/usr/bin/perl -w

use strict;
use Getopt::Std;

my %switches;
my $cdsinmrna = 0;
my $putframe = 0;
my $reverse = 0;
my %estinfo = ();

getopts("mhfre:",\%switches);

if (exists($switches{h})) {
  print "est2genome2gff.pl [-mhfr -e <filename>] <cdsstart> <cdsend> <est2genomeoutputfile> [<refgenecoord>]\n";
  print "\n  this program parses output files from est2genome into gff format. options are:\n";
  print "  -m cds coordinates are mRNA-based (overrides coordinates from refgene file)\n";
  print "  -h this help\n";
  print "  -f include frame\n";
  print "  -r reverse cds coordinates if strand is negative in refgene file\n";
  print "  -e <filename> additional information about the EST in <filename> (GenBank format)\n";
  exit(1);
}

if (exists($switches{m})) {
  $cdsinmrna = 1;
}

if (exists($switches{f})) {
  $putframe = 1;
}

if (exists($switches{r})) {
  $reverse = 1;
}

if (exists($switches{e})) {
  %estinfo = read_estinfo($switches{e});
}

if (scalar(@ARGV) < 3) {
  print "est2genome2gff.pl [-mhfr -e <filename>] <cdsstart> <cdsend> <est2genomeoutpufile> [<refgenecoord>]\n";
  exit(1);
}

if (!open(SEQ,"< $ARGV[2]")) {
  print "est2genome2gff.pl: impossible to open $ARGV[2]\n";
  exit(1);
}

my $cdsSt = $ARGV[0];
my $cdsEn = $ARGV[1];
my $incds = 0;
my $rmndr = 0;

if (scalar(@ARGV) > 3) {
  if (!$cdsinmrna && open(REF,"< $ARGV[3]")) {
    my $buf = <REF>;

    if ($buf =~ m/[\w\d_]+\s+[\w\d]+\s+([\+\-])\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/) {
      my $st = $1;
      my $gS = $2;
      my $gE = $3;
      my $cS = $4;
      my $cE = $5;

      if (!$reverse || $st eq '+') {
        $cdsSt = $cS - $gS;
        $cdsEn = $cE - $gS;
      } else {
        $cdsSt = $gE - $cE;
        $cdsEn = $cdsSt + ($cE - $cS);
      }
    } else {
       print "est2genome.pl: wrong refgene format in $ARGV[3]\n";
       exit(1);
    }
  }
}

if ($cdsSt >= $cdsEn) {
  print "est2genome2gff.pl: wrong cds coordinates <cdsstart> >= <cdsend>\n";
  close(SEQ);
  exit(1);
}

while (<SEQ>) {

  if (!m/([\w\+\-\d_]+)\s+\d+\s+(\d+\.\d+)\s+(\d+)\s+(\d+)\s+[\w\d\._]+\s+(\d+)\s+(\d+)\s+([\w\._]+)\s+([^\n]+)/) {
    #print STDERR "unmatched: $_";
    next;
  }

  my $srce;
  my $atype = $1;
  my $score = $2;
  my $genSt = $3;
  my $genEn = $4;
  my $estSt = $5;
  my $estEn = $6;
  my $estDe1 = $7;
  my $estDe2 = $8;
  my $frgSt = $cdsinmrna ? $estSt : $genSt;
  my $frgEn = $cdsinmrna ? $estEn : $genEn;

  if ($atype eq "Exon") {

    chomp($estDe2);

    my $estTy = "";
    #my $estId = substr($estDe2,0,30);
    my $exSt;
    my $exEn;
    my $fram = ".";
    my $exTy = "";
    my $estId;

    #$estId =~ s/\s//g;

    if ($estDe2 =~ m/3\'/ || $estDe2 =~ m/3 \'/) {
      $estTy = "3p";
    } elsif ($estDe2 =~ m/5\'/ || $estDe2 =~ m/5 \'/) {
      $estTy = "5p";
    }

    if ($estDe1 =~ m/gb_(\w+)/) {
      $estId = $1;
      $srce = $estId;
      if (exists($estinfo{$estId})) {
        $estDe2 = "\"$estinfo{$estId}{organ}\___$estinfo{$estId}{tissu}\___$estId\___$estinfo{$estId}{devst}\"";
      } else {
        $estDe2 = $estId;
      }
    } else {
      $estId = substr($estDe1,0,30);
      $srce = substr($estId,0,10);
    }

    if ($cdsSt >= $frgSt && $cdsEn <= $frgEn) {  # cds within an exon

      if ($frgSt < $cdsSt) { # non-coding
        $exTy = "utr" . $estTy;
        $fram = ".";
        if ($cdsinmrna) {
          $exSt = $genSt;
          $exEn = $cdsSt - $estSt + $genSt - 1;
        } else {
          $exSt = $genSt;
          $exEn = $cdsSt - 1;
        }
        print "$estId\t$srce\t$exTy\t$exSt\t$exEn\t$score\t+\t$fram\t$estDe2\n";
      }

      # coding
      $exTy = "single" . $estTy; # this is in fact a single exon thingy
      $fram = $putframe ? 0 : ".";
      if ($cdsinmrna) {
        $exSt = $cdsSt - $estSt + $genSt;
        $exEn = $cdsEn - $estEn + $genEn;
      } else {
        $exSt = $cdsSt;
        $exEn = $cdsEn;
      }
      print "$estId\t$srce\t$exTy\t$exSt\t$exEn\t$score\t+\t$fram\t$estDe2\n";

      # non-coding
      if ($frgEn > $cdsEn) {
        $exTy = "utr" . $estTy;
        $fram = ".";
        if ($cdsinmrna) {
          $exSt = $cdsEn - $estEn + $genEn + 1;
          $exEn = $genEn;
        } else {
          $exSt = $cdsEn + 1;
          $exEn = $genEn;
        }
        print "$estId\t$srce\t$exTy\t$exSt\t$exEn\t$score\t+\t$fram\t$estDe2\n";
      }

    } elsif ($cdsSt >= $frgSt && $cdsSt <= $frgEn) { # cdsstart within an exon

      # non-coding
      $exTy = "utr" . $estTy;
      $fram = ".";
      if ($cdsinmrna) {
        $exSt = $genSt;
        $exEn = $cdsSt - $estSt + $genSt - 1;
      } else {
        $exSt = $genSt;
        $exEn = $cdsSt - 1;
      }
      print "$estId\t$srce\t$exTy\t$exSt\t$exEn\t$score\t+\t$fram\t$estDe2\n";

      # coding
      $exTy = "initial" . $estTy;
      $fram = $putframe ? 0 : ".";
      if ($cdsinmrna) {
        $exSt = $cdsSt - $estSt + $genSt;
        $exEn = $genEn;
      } else {
        $exSt = $cdsSt;
        $exEn = $genEn;
      }
      print "$estId\t$srce\t$exTy\t$exSt\t$exEn\t$score\t+\t$fram\t$estDe2\n";

      $incds = 1;

      if ($putframe) {
        $rmndr = ( 3 - ( $exEn - ( $exSt + $fram) + 1) % 3 ) % 3;       # rtfm (gff2ps)
      }

    } elsif ($cdsEn >= $frgSt && $cdsEn <= $frgEn) { # cdsend within an exon

      # coding
      $exTy = "terminal" . $estTy;
      $fram = $incds && $putframe ? $rmndr : "." ;
      if ($cdsinmrna) {
        $exSt = $genSt;
        $exEn = $cdsEn - $estSt + $genSt;
      } else {
        $exSt = $genSt;
        $exEn = $cdsEn;
      }
      print "$estId\t$srce\t$exTy\t$exSt\t$exEn\t$score\t+\t$fram\t$estDe2\n";

      # non-coding
      $exTy = "utr" . $estTy;
      $fram = ".";
      if ($cdsinmrna) {
        $exSt = $cdsEn - $estSt + $genSt + 1;
        $exEn = $genEn;
      } else {
        $exSt = $cdsEn + 1;
        $exEn = $genEn;
      }
      print "$estId\t$srce\t$exTy\t$exSt\t$exEn\t$score\t+\t$fram\t$estDe2\n";

      $incds = 0;
      $rmndr = 0;

    } elsif ($cdsSt < $frgSt && $cdsEn > $frgEn) { # exon within cds

      $exTy = "internal" . $estTy;
      $fram = $incds && $putframe ? $rmndr : ".";
      $exSt = $genSt;
      $exEn = $genEn;
      print "$estId\t$srce\t$exTy\t$exSt\t$exEn\t$score\t+\t$fram\t$estDe2\n";

      if ($incds && $putframe) {
        $rmndr = ( 3 - ( $exEn - ( $exSt + $fram) + 1) % 3 ) % 3;       # rtfm (gff2ps)
      }

    } else {                                       # exon out of cds

      $exTy = "utr" . $estTy;
      $fram = ".";
      $exSt = $genSt;
      $exEn = $genEn;
      print "$estId\t$srce\t$exTy\t$exSt\t$exEn\t$score\t+\t$fram\t$estDe2\n";

    }
  }

}


#
# functions
#

sub read_estinfo {
  my %einfo;
  my $fname = $_[0];

  if (!open(FH,"< $fname")) {
    return ();
  }

  my $gbid="x";
  my $organ="org.unk";
  my $devst="devst.unk";
  my $tissu="tiss.unk";

  while (<FH>) {

    if (/\AGenBank Acc:\s+(\w+)/) {
      $gbid=$1;
    } elsif (/\AOrgan:\s+([^\n]+)/) {
      $organ="org.$1";
      $organ=~s/ /_/g;
    } elsif (/\ADevelop\. stage:\s+([^\n]+)/) {
      $devst="devst.$1";
      $devst=~s/ /_/g;
    } elsif (/\ATag Tissue:\s+([^\n]+)/) {
      $tissu="tiss.$1";
      $tissu=~s/ /_/g;
    } elsif (/\A\|\|/) {
      if (exists($einfo{$gbid})) {
        print "$gbid already exists!!!!\n";
      } else {
        $einfo{$gbid}{organ}=$organ;
        $einfo{$gbid}{devst}=$devst;
        $einfo{$gbid}{tissu}=$tissu;
        $gbid = "x";
        $organ = "org.unk";
        $devst = "devst.unk";
        $tissu = "tiss.unk";
      }
    }

  }

  return %einfo;
}
