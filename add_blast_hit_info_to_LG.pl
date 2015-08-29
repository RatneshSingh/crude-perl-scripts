#!/usr/bin/perl
use strict;
use Getopt::Long;

our ( $gstab, $gmtab, $lgtab, $help, $out, $blast, $gfffile,$forcegmod,$replace,$with,$col,$printall );
$col=0;
my $result = GetOptions(
    "gstab|s=s" => \$gstab,
    "gmtab|g=s" => \$gmtab,
    "gff|f=s"   => \$gfffile,
    "lgtab|l=s" => \$lgtab,
    "forcegmod|m"=>\$forcegmod,
    "out|o=s"   => \$out,
    "blast|b=s"   => \$blast,
    "replace|r=s"=>\$replace,
    "with|w=s"=>\$with,
    "help|h"    => \$help,
    "col|c=i"   =>\$col,
    "printall|p"=>\$printall
);

my $usage = "
perl $0 -s GenoSeqTable -g GeneModTable -f gff3File -l LinkageTable -b blasttableWidSeqLengths -o outfile

where:
GenoSeqTable contains:
MarkerName\tChromosomeNumer\tlocationOnChromosome

GeneModTable contains:
MarkerName\tGeneModelName

LinkageTable contains:
MarkerName\tCentiMorghan\tLinkageGroupNumber

blasttableWidSeqLengths
BlastTable+QCov\tSCov\tQlen\tSubLen

Options:
-m force genemodel location over Chromosome hit location if marker hits both.
-h print help and exit.
-o output file to save data.
-b reference genome blast output in table format to get Chromosome sizes.
-r replace this with \"with|w\" in gene model name in GFF file. To avoid discrepancies due to dot or underscore in names.
-w Replace with this.
-c Use column num as name if names are seperated by \"|\" in genemod blastfile
-p Print all markers and their mapping info. (Default is to print only the ones on Linkage group)
";

die "$usage" if $help;
die "$usage" if ( !$gstab || !$gmtab || !$lgtab );
$out = $out ? $out : "$lgtab.markinfo.table";

#my$gstab=$ARGV[0]; # genomic sequence
#my$gmtab=$ARGV[1];  # gene models
#my$lgtab=$ARGV[2];  # linkage info
my %markinfo;
my%gffinfo;
my %chr_conv;  ## for conversion of coordinates for plotting
my %chr_size;
my %chrlist;
open GMTAB, "$gmtab" or die "\nCannot open GeneModel table (gmtab):$gmtab\n$usage\n\n";
open GSTAB, "$gstab" or die "\nCannot open GenomicSeq table(gstab):$gstab\n$usage\n\n";
open GFF, "$gfffile" or die "\nCannot open gff file(gff):$gfffile\n$usage\n\n";
open LGTAB, "$lgtab" or die "\nCannot open Linkage Group table(lgtab):$lgtab\n$usage\n\n";
open OUT,   ">$out" or die "\nCannot open Output file:$out\n$0\n$usage\n\n";
open BLAST, "$blast" or die "\nCannot open GenomicBlast blast table(blast):$blast\n$usage\n\n";
open OUT2, ">$out.ModCoord.table";
##Read whole mtab to hash

while (<GSTAB>) {
    next if /^\s+$/;
    my @line = split(/\s+/);
    chomp @line;
    for (@line) { s/\s+//; }

    $markinfo{ $line[0] }{'chr'}    = $line[1];
    $markinfo{ $line[0] }{'chrloc'} = $line[2];
    $chrlist{$line[1]}=1;
}

while (<GFF>) {
    # chr[0]  method[1]  type[2]  start[3] end[4] dot[5] strand[6]  phase[7]  desc[8]
    next if /^\s*$|^#/;
    my @gff = split(/\s+/);
    chomp @gff;
    next if $gff[2] !~ m/mRNA/i;
    for (@gff) { s/\s+//; }
    my$name=$1 if $gff[8]=~/Name\=([^;]+)/;
    my$id=$1 if $gff[8]=~/ID\=([^;]+)/;
    $name=~s/$replace/$with/gi if ($replace && $with);
    $id=~s/$replace/$with/gi if ($replace && $with);
    $gffinfo{$name}{start}=$gff[3];# if $gff[2] =~ m/mRNA/i;
    $gffinfo{$name}{end}=$gff[4];# if $gff[2] =~ m/mRNA/i;
    $gffinfo{$name}{chr}=$gff[0];#  if $gff[2] =~ m/mRNA/i;

    $gffinfo{$id}{start}=$gff[3];#  if $gff[2] =~ m/mRNA/i;
    $gffinfo{$id}{end}=$gff[4];#  if $gff[2] =~ m/mRNA/i;
    $gffinfo{$id}{chr}=$gff[0];#  if $gff[2] =~ m/mRNA/i;
}




while (<GMTAB>) {
    next if /^\s+$/;
    my @line = split(/\s+/);
    chomp @line;
    for (@line) { s/\s+//; }
    my$gm=$line[1];

    my@temp=split(/\|/,$line[1]);


    $gm=$temp[$col-1] if $line[1]=~/\|/;
    $markinfo{ $line[0] }{'gmod'} = $gm;
}







### LGTAB: MarkerName [0]  \t  CentiMorghan [1]  \t  LinkageGroupNumber[2]
my%maxlg;
my%lg_conv;
while (<LGTAB>) {
    next if /^\s+$/;
    my @line = split(/\s+/);
    chomp @line;
    for (@line) { s/\s+//; }

    $markinfo{ $line[0] }{'cm'} = $line[1];
    $markinfo{ $line[0] }{'lg'} = $line[2];

    my$lgnum=$line[2];
    $lgnum=~s/\D//g;
    $lg_conv{$lgnum}=$line[2];
    $maxlg{$line[2]}=0 if !$maxlg{$line[2]};
    $maxlg{$line[2]}=$line[1] if $maxlg{$line[2]}<$line[1];
    #print "\nZero Cm:$line[1] with:$markinfo{ $line[0] }{'cm'} " if $line[1]==0;
}


### find the appropriate order of chromosome to plot.
my%chrorder;
my%homolog;
my@Chr_lg_ord;
foreach my $marker ( keys %markinfo ) {$chrorder{$markinfo{ $marker }{'lg'}}{$markinfo{ $marker }{'chr'}}++;}
foreach my $lg(keys %chrorder){
    next if $lg=~/^\s*$/;
    #print "\nFor LG:$lg\n";
    my$max=0;
    my$chrlg;
    my$numlg=$lg;
    $numlg=~s/[\D]+//g;
    foreach my$chr(keys %{$chrorder{$lg}}){
        next if $chr=~/^\s*$/;
        #print "\n\tFor Chr:$chr freq is:$chrorder{$lg}{$chr}\n";
        if($chrorder{$lg}{$chr} >= $max){
            $max=$chrorder{$lg}{$chr};
            $chrlg=$chr;
        }
        else{next}
    }
    print "\nLinkage group $lg Homologous to $chrlg";
    $homolog{$numlg}=$chrlg;#$lg;

}
##add chr names in order of homology to LGs
my%tempchrlist=%chrlist;
foreach my$numlg(sort {$a<=>$b} keys %homolog){push(@Chr_lg_ord,$homolog{$numlg}) if exists $tempchrlist{$homolog{$numlg}};delete $tempchrlist{$homolog{$numlg}}; }
## add rest of the chromosome for which no homology was found.
foreach my$restchr(sort{(my$x=$a)=~s/\D+//g;(my$y=$b)=~s/\D+//g;($x?$x:1e10) <=> ($y?$y:1e10)} keys  %tempchrlist){next if $restchr=~/^\s*$/; push(@Chr_lg_ord,$restchr)};




### create modified coordinates need to be added to original coord for plotting.
my$prev_lgsize=0;
my%modlg_size;
print OUT2 "Num\tLinkageGroup\tOriginalCM\tAddedToCM\t \tX-axis\tY-axis(LGBoundries)\t \t \t \tLabel(X-axis)\tY-axis";
foreach my$nums(sort{$a<=>$b} keys %lg_conv){
    next if $nums < 1;
    $modlg_size{$lg_conv{$nums}}=$prev_lgsize;
    $prev_lgsize+=$maxlg{$lg_conv{$nums}};
    print OUT2 "\n$nums\t$lg_conv{$nums}\t$maxlg{$lg_conv{$nums}}\t$modlg_size{$lg_conv{$nums}}\t \t0\t$prev_lgsize\t \t \t$lg_conv{$nums}\t",($maxlg{$lg_conv{$nums}}/2)+$modlg_size{$lg_conv{$nums}},"\t0";
}




while (<BLAST>) {
    next if /^\s*$/;
    s/^\s+|\s+$//g;
    my @line = split(/\s+/);
    chomp @line;
    for (@line) { s/\s+//; }

## for conversion of coordinates to plotable data
      my$chr_num=$line[1];
      $chr_num=~s/\D+//g;
      $chr_conv{$chr_num}=$line[1];
      $chr_size{$line[1]}=$line[15];
}

### create modified coordinates need to be added to original coord for plotting.
my$prev_size=0;
my%mod_size;
print OUT2 "\n\nNum\tChromosome\tOriginalSize\t \t \tX-axis\tY-axis(ChrBoundries)\t \t \t \tX-axis\tLabel(Y-axis)";
#foreach my$nums(sort{$a<=>$b} keys %chr_conv){
#    next if $nums < 1;
#    $mod_size{$chr_conv{$nums}}=$prev_size;
#    $prev_size+=$chr_size{$chr_conv{$nums}};
#
#    print OUT2 "\n$nums\t$chr_conv{$nums}\t$chr_size{$chr_conv{$nums}}\t$prev_size\n";
#}
my$chrnum=0;
foreach my$chrs(sort{$a<=>$b} @Chr_lg_ord){
    #next if $chrs < 1;
    #$chrs=~s/\D+//g;
    $mod_size{$chrs}=$prev_size;
    $prev_size+=$chr_size{$chrs};
    $chrnum++;
    print OUT2 "\n$chrnum\t$chrs\t$chr_size{$chrs}\t \t$chrs\t0\t$prev_size\t \t \t$chrs\t0\t",$mod_size{$chrs}+($chr_size{$chrs}/2);
    #print "\n$chrs\t$chrs\t$chr_size{$chrs}\t$prev_size";

}







print OUT "Marker\tLinkageGroup\tCentiMorghan\tChr\tChr_loc\tGeneModel\tGenModOnChrom\tGenModCoord\t \tModLGCoordinates\tModChrCoordinates";
for my $marker ( keys %markinfo ) {


    next if !exists $markinfo{$marker}{'cm'} && !$printall;
    #print "\nMarker:$marker has Zero:$markinfo{$marker}{'cm'} value" if ( $markinfo{$marker}{'cm'}== 0);
    my$modchrsize;
    my$modlgsize;

    $modchrsize=$markinfo{$marker}{'chrloc'}?$markinfo{$marker}{'chrloc'} + $mod_size{$markinfo{$marker}{'chr'}}:" " if ($markinfo{$marker}{'chrloc'});
    $modlgsize=$markinfo{$marker}{'cm'}>=0?$markinfo{$marker}{'cm'} + $modlg_size{$markinfo{$marker}{'lg'}}:" "  if ($markinfo{$marker}{'cm'});
    if($forcegmod){$modchrsize=$gffinfo{$markinfo{$marker}{'gmod'}}{start}?$gffinfo{$markinfo{$marker}{'gmod'}}{start} + $mod_size{$gffinfo{$markinfo{$marker}{'gmod'}}{chr}}:$modchrsize;}

    ## print markers ones with LG and chr info only.
    next if !$printall && ($markinfo{$marker}{'lg'}=~/^\s*$/ || $markinfo{$marker}{'chr'}=~/^\s*$/);

    my$result=join("\t",
                   $marker,
                   $markinfo{$marker}{'lg'}?$markinfo{$marker}{'lg'}:" ",
                   $markinfo{$marker}{'cm'}>=0?$markinfo{$marker}{'cm'}:" ",
                   $markinfo{$marker}{'chr'}?$markinfo{$marker}{'chr'}:" ",
                   $markinfo{$marker}{'chrloc'}?$markinfo{$marker}{'chrloc'}:" ",
                   $markinfo{$marker}{'gmod'}?$markinfo{$marker}{'gmod'}:" ",
                   $gffinfo{$markinfo{$marker}{'gmod'}}{chr}?$gffinfo{$markinfo{$marker}{'gmod'}}{chr}:" ",
                   $gffinfo{$markinfo{$marker}{'gmod'}}{start}?$gffinfo{$markinfo{$marker}{'gmod'}}{start}:" ",
                   " ",
                   $modlgsize,
                   $modchrsize
                   );


    print OUT "\n$result";
    print "\n$result" if ( $markinfo{$marker}{'cm'}== 0)

}
