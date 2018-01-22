#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;


our($cat,$loc,$mstmap,$help,%markers);

GetOptions(
  "catalog|cat|c=s"=>\$cat,
  "loci|l=s"=>\$loc,
  "mstmap|m=s"=>\$mstmap,
  "help|h"=>\$help
);

#my$cat=$ARGV[0];
#my$loc=$ARGV[1];

my$usage="\n$0 -c catalog.tags.tsv  -l joinmap_file.loc


  catalog|cat|c   catalog tags tsv file which contains location information,

  ## input file to transfer location on. use only one.
  loci|l          joinmap loci file
  mstmap|m        mstmap file.


  help|h

\n";

open(CAT,"$cat") or die "\nUnable to open file $cat\n$usage\n";
my$outfile;
if ($loc) {
  open(LOC,"$loc") or die "\nUnable to open file $loc\n$usage\n";
  $outfile="$loc.widLoc.loc";
  open(OUT,">$loc.widLoc.loc") or die "\nUnable to open out file for writing\n";
}
elsif($mstmap){
  open(MSTMAP,"$mstmap") or die "\nUnable to open file $mstmap\n$usage\n";
  $outfile="$mstmap.widLoc.output";
  open(OUT,">$mstmap.widLoc.output") or die "\nUnable to open out file for writing\n";
}


while (<CAT>) {
  s/^\s+//g;
  my@catline=split /\s+/;
  $markers{$catline[2]}{chr}=$catline[3];
  $markers{$catline[2]}{loc}=$catline[4];
}
close CAT;

if ($loc) {
  our($prln);
  while (<LOC>) {
   print OUT $_ if ($. < 5 || m/^\s*$/);
    next if ($. < 5 || m/^\s*$/);
   #print OUT $_ if m/^\s*$/;
    #next if m/^\s*$/;
   s/^\s+//g;
   $prln=2 if m/^individual names:/;

    if (!$prln) {
      my@locline=split /\s+/;
     $locline[0]=~s/_[^\s]+//g;
     my$tmp=join "_",$locline[0],$markers{$locline[0]}{chr},$markers{$locline[0]}{loc};
     $locline[0]=$tmp;
     print OUT join " ",@locline,"\n";
   }
   print OUT $_ if $prln;
  }

}elsif($mstmap){
  my$prln=0;

  while (<MSTMAP>) {
    s/^\s+//g;

    $prln=0 if m/^;ENDOFGROUP/;

    if ($prln==2) {
      my@locline=split /\s+/;
      $locline[0]=~s/_[^\s]+//g;
      my$tmp=join "_",$locline[0],$markers{$locline[0]}{chr},$markers{$locline[0]}{loc};
      $locline[0]=$tmp;
      print OUT join "\t",@locline,"\n";
    }
    print OUT $_ if $prln == 0;
    $prln=2 if m/^;BEGINOFGROUP/;
  }

}


print "\nresults are saved as $outfile\n";
close OUT;
close MSTMAP if $mstmap;
close LOC if $loc;