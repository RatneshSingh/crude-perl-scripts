use strict;
use warnings;

open HAPLO,"$ARGV[0]" ;
our(@headers,%alpersamp);

print "Catalog_Id\tCnt\t$ARGV[0]\n";
while (<HAPLO>){

    if (/^Cat/){
      @headers=split(/\t/);
      s/\s+//g foreach @headers;
      next;
    }
    #next if /-/;
    next if tr/-/-/ > $ARGV[1];
    my$haplo=$_;
    my @line=split(/\s+/,$haplo);
    s/\s+//g foreach @line;
    my%alleles;
    for(my$i=2; $i < scalar @line; $i++){
      next if $line[$i]=~/-/;
      #next if $line[$i]=~/consensus/;
      my$numall=$line[$i]=~ tr/\//\//;
      $alpersamp{$headers[$i]}{$numall+1}+=1;
      $alleles{$_}=1 foreach split /\//,$line[$i];
    }

    my$big_allele=join("/",keys %alleles);
    #next if $big_allele=~/^\s*$/g;
    if ($big_allele=~/^\s*$/g){
    print "\n$line[0]\t$line[1]\t-\t0"
    }else{
    print "\n$line[0]\t$line[1]\t$big_allele\t",scalar keys %alleles;
    }


}
print "\n";
exit;
print "\nSample_ID\tNum_Allele\tNum_Loci\tAverage_allele_per_loci\n";#$ARGV[0]\n";
foreach my$samp(keys %alpersamp){
  my$totalal=0;
  my$totalloc=0;
  print "$samp\t";
  foreach my$alnum(sort keys %{$alpersamp{$samp}}){
   # print "$samp\t",$alnum+1,"\t$alpersamp{$samp}{$alnum}\n";
    $totalal +=  $alpersamp{$samp}{$alnum}*$alnum;
    $totalloc+=$alpersamp{$samp}{$alnum};
  }
  print "\t$totalal\t$totalloc\t",$totalal/$totalloc,"\n";
}