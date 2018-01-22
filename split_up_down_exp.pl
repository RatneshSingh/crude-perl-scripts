use strict;
use warnings;
use Getopt::Long;
our($file,$type,$val,$help,$up,$down,$namcol,$col,$out,$log,$exclude,$maxnum,$nofc);
$val=1;
$namcol=1;
$col=2;

GetOptions(
  "file|f=s"=>\$file,
  "out|o=s"=>\$out,
  "col|c=i"=>\$col,          ### column number to evaluate expression from
  "ratio|r=f"=>\$val,        ### fold change cutoof value to filter for
  "namecol|nc=i"=>\$namcol,  ### column containing gene name to print.
  "help|h"=>\$help,
  "up|u"=>\$up,              ### filter for upregulated genes
  "down|d"=>\$down,           ### filter for down-regulated genes
  "log|l=i"=>\$log        ,    ### fold change are as log values
  "exclude|e=s"=>\$exclude,    ### exclude genes with this pattern
  "maxnum|m=i"=>\$maxnum,      ### max num of genes to print. top hits will be printed
  "nofc"=>\$nofc
);

my$usage="
usage: perl $0 options....
  file|f      file with expression fold change data,
  namecol|nc  column number containing gene names to print[1].
  col|c       column number to evaluate expression from[2]
  out|o       output file to save[column name+expr_vall+up/down]
  ratio|r     fold change cuto off value to filter for[1]
  exclude|e   eclude gene with this pattern in name.
  maxnum|m    Only print top this many genes if there more.
  up|u        print up-regulated genes
  down|d      print down-regulated genes
  log|l       expression values are in log
  nofc        do not print fold change. only gene names.
  help|h      help
";

$namcol-=1;
$col-=1;
my$val2=$val;
$val2=$log**abs($val) if $log;
my%maxhit;
die $usage if $help;
my$extra_add=".Upregulated_${val2}x" if $up;
$extra_add=".Downregulated_${val2}x" if $down;
open(EXP,$file) or die "Could not open file\n\n$usage";

my@selected;
my$title;
while (<EXP>) {
  next if m/^\s*$/;
  my@list=split /\s+/;
  next if !$list[$col];
  if ($. == 1 ) {
      $title=$list[$col] if $. ==1;
      push (@selected,"$list[$namcol]\t$list[$col]") if ($. == 1 && !$nofc);
      push (@selected,"$list[$col]") if ($. == 1 && $nofc);
  }else{
      next if ($exclude && $list[$namcol]=~m/"$exclude"/i);
      my$cur_gene= "$list[$namcol]\t$list[$col]";
      $cur_gene= $list[$namcol] if $nofc;
      push (@selected,$cur_gene) if ($list[$col] >= $val && $up);
      push (@selected,$cur_gene) if ($list[$col] <= $val && $down);
      if ($maxnum) {
        $maxhit{$list[$namcol]}=$list[$col] if ($list[$col] >= $val && $up);
        $maxhit{$list[$namcol]}=$list[$col] if ($list[$col] <= $val && $down);
      }

  }
}



$title.=$extra_add;
$title=~s/[\-\>\<]+/_/g;
$title.=".GeneList" if $nofc;
$out="$title.table" if !$out;


if ($maxnum) {
  print "\nlimiting gene hits to $maxnum number\n";
  my@sorted=sort{$maxhit{$a}<=>$maxhit{$b}} keys %maxhit if $down;
  @sorted=sort{$maxhit{$b}<=>$maxhit{$a}} keys %maxhit if $up;
  my@printsort;
  foreach(@sorted){if($nofc){push(@printsort,"$_")}else{push(@printsort,"$_\t$maxhit{$_}")}}
  @selected=($selected[0],@printsort) if scalar@printsort < $maxnum;
  @selected=($selected[0],@printsort[0..$maxnum-1]) if scalar@printsort >= $maxnum;

  print "\nTotal number of $title genes:",scalar@printsort,"\tGenes in the new list:",scalar@selected,"\n";
}

if (scalar@selected > 1){
  open(OUT, ">$out") or die "Unable to write results\n\n$usage";
  print OUT join("\n",@selected) if scalar@selected > 1;
}
else{
  die "\nExiting.No genes qualified the criteria\n";
}
close OUT;