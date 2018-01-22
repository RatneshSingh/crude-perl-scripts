#!/usr/bin/perl
use Data::Dumper;
#perl -e '
use warnings;
use strict;
use Getopt::Long;

our($lev,$excl,@files,$help,$out,$idfier,$sublist,$expression,$split_name,$ncol,$tcol,$get_lev,$verbose);

my$usage="
perl $0 -options....

options:
-l	level to process. [4]
-e	exclusively process level. No parent level in the result.
-i	use gene identifier eg. trinity. Helps if getting error for erroneously assigned gene ids..
-o	output file to save result[Auto]
-s	Name or number of mercator bin to extract names of genes in that bin. use multiple seperated with comma.
-gexp	Gene expression data file to extract expression data related to selected bin in -s.
-glev	Bin level of which genes expression need to be pulled off [1]
-t	column in gene expression table to be used as gene name to match [1]
-d	Split the name in gene expression with this parameter and use -col as name.
-c	col num to be used as name after splitting the name by -d.[1].
	selecting multiple columns allowed by  providing column numbers seperated by comma.
	eg: -c 1,3,4  -> join columns 1,3,4 by delim to create name for comparison.
-h	print help and exit.
-f	all the mercator files. Multiple files could be provide with -f each or together in a row seperated by spaces.
";


## some default values
$get_lev=1;



GetOptions(
  "level|l=i" => \$lev,    ## inclusive level. Print infor on this level and everything below.
  "excl|e" => \$excl,  ## exclusive level flag. Print infor on only this level.
  "mercator|m|f=s" => \@files,
  "out|o=s" => \$out,
  "identifier|i=s" => \$idfier,
  "sub_list|s=s"=> \$sublist,
  "geneexp|gexp=s"=>\$expression,
  "gene_lev|glev=i"=>\$get_lev,
  "delim|d=s"=>\$split_name,
  "col|c=s"=>\$ncol,
  "tcol|t=i"=>\$tcol,
  "help|h" => \$help,
  "verbose|v"=>\$verbose
);


 die "\n\n$usage\n\n" if $help;
die "\n-gexp is useless unless -s flag s set too.... \n$usage\n\n" if ($expression && ! $sublist);

  $lev=$lev?$lev:4;
  #$get_lev=$get_lev?$get_lev:scalar(split(/\./,$sublist));

  $tcol=$tcol?$tcol:1;
  $ncol=$ncol?$ncol:1;
  $out=$out?$out: join ".","Mercator_processed.Level$lev",$excl?"excl.table":"table";


  if (@ARGV) {
    push @files,grep{ -f $_ } @ARGV;
  }

  die "\nNo input files found\n$usage\n\n" if scalar@files < 1;

  my%bins_id;
  my%bins_desc;
  my%fnumgene;
  my%merc_data;
  my%id_store;
  my%gene_desc_all;
  my%gene_desc_lev;
  my $OUT;
  open($OUT, ">",$out);
  foreach my$file(@files){
    print "\nProcessing file: $file";
    #my ($fh);
    next unless ( -e $file);
    open (my $fh,$file);
    my$numgene=0;
    my%all_gene;
    while(<$fh>){
      next if m/^\s*$/;
      my($bid,$bdesc,$gid,$gdesc)=split /\t/;
      substr $bid, 0, 1, "";chop $bid;
      substr $bdesc, 0, 1, "";chop $bdesc;
      substr $gid, 0, 1, "";chop $gid;
      substr $gdesc, 0, 1, "";chop $gdesc;

      if($idfier && $gid !~ /$idfier/){next;}  ## to use gene identifies as filter. Some times gid is assigned from non gene elements.
      if($gid ne ""){
        $numgene++ if !$all_gene{$gid};$all_gene{$gid}=1;
        my@bid=split /\./,$bid;


        my@bdesc=split /\./,$bdesc;

        for my$i(0 .. $#bid){
          my$ibid=$bid[$i];
            $ibid=join ".", @bid[0 .. $i] if $i > 0;
            $ibid=~s/\.$//g;
          my$ibdesc=$bdesc[$i];
            $ibdesc=join ".", @bdesc[0 .. $i] if $i > 0;
	    ### print commands below are for diagnosis purpose only.
	    print "\n print for BDesc info:i:$i\tBDescNumElem:",scalar@bdesc,"\tContentElement:",join ":",@bdesc if $i > $#bdesc;
	    print "\n print for Bid info:i:$i\tBidNumElem:",scalar@bid,"\tContentElement:",join ":",@bid if $i > $#bdesc;
            $ibdesc=~s/\.$//g;

          $bins_id{$ibid}{$file}{num}++;
          $bins_desc{$ibdesc}{$file}{num}++;
          $bins_id{$ibid}{$file}{lev}=$i+1;
          $bins_desc{$ibdesc}{$file}{lev}=$i+1;

          $merc_data{$ibid}{desc}=$ibdesc;
          $merc_data{$ibid}{id}=$ibid;
          $merc_data{$ibid}{lev}=$i+1;
          $merc_data{$ibdesc}{id}=$ibid;
          $merc_data{$ibdesc}{lev}=$i+1;
          $merc_data{$ibdesc}{desc}=$ibdesc;


          ## collect gene names in each bin to isolate them if user demands.
          push(@{$id_store{$ibid}},$gid);
          push(@{$id_store{$ibdesc}},$gid);

          push(@{$gene_desc_all{$gid}},$ibdesc);
          ${$gene_desc_lev{$gid}}[$i+1]=$ibdesc;

        }
      }else{next;#print "\n No hits found for $_";
      }
    }
    $fnumgene{$file}=$numgene;
    close $fh;
  }
  print "\nFinished processing files\n";
  print $OUT join "\t","Level","Bin_ID","Bin_Description",(map{$_."_Num",$_."_Percent"}@files),"\n";
  for my$fbid(sort keys %bins_id){
    my$nline=0;
    next if $merc_data{$fbid}{lev} > $lev;
    next if ($excl && $merc_data{$fbid}{lev} != $lev);  ### print exclusive level if user asked so.
    print $OUT join "\t",$merc_data{$fbid}{lev},$fbid,$merc_data{$fbid}{desc},"";
    foreach my$file(@files){
      #print "\nfbid:$fbid\nfile:$file\nnumgene:$fnumgene{$file}\nbins_id\{$fbid\}\{$file\}\{num\}:$bins_id{$fbid}{$file}{num}\n";

      #if ($bins_id{$fbid}{$file}{lev} == $lev){
        $bins_id{$fbid}{$file}{num}=$bins_id{$fbid}{$file}{num}?$bins_id{$fbid}{$file}{num}:0;
        print $OUT join "\t",$bins_id{$fbid}{$file}{num},sprintf("%.2f",$bins_id{$fbid}{$file}{num}*100/$fnumgene{$file}),"";
        $nline=1;
      #}

    }
    print $OUT "\n" if $nline == 1;
  }
 print "\n";
close $OUT;

#print Dumper %merc_data;
### if expression data is provided, save them as hash to create expression table for specific bins.
my$EXP;
my%exp_data;
if ($expression) {
    open($EXP,"<",$expression) or print "Unable to open Gene xpression file";
    while (<$EXP>) {
	chomp;
        my$line=$_;
        $exp_data{'header'}=$line if $. ==1;
        next if $. ==1;
        my@ext=split(/\t/,$line);
        my$name=$ext[$tcol-1];

        ## split name if user asks.
        if ($split_name) {$name=getname($name,$split_name,$ncol)}

        push(@{$exp_data{lc$name}},$line);
     }
}


### print selected gene names for bin requested.
if ($sublist) {
  my$OUT2;
  my $printed=0;
  my@msublist=split(/\,/,$sublist);
  print '@msublist:' . Dumper \@msublist;
  my $out_file="$out.Bin.".join("_",scalar@msublist < 4?@msublist:scalar@msublist."variousBins","Lev",$get_lev)."expression.list";
  open($OUT2, ">", $out_file) or die "Unable to open file to save list";
  print $OUT2 "BinID\tBinDesc\tGeneID\tGene_Annotation\tGene_All_annotations",$exp_data{'header'}?"\t$exp_data{'header'}":"","\n";



  ### fix for getting genes at certain level.
  #print Dumper @msublist;
  foreach my$esublist(@msublist){
    my$all_ids_at_level=get_levels(\%merc_data,$esublist,$get_lev,1);
    #print Dumper @$all_ids_at_level;
    next if scalar @$all_ids_at_level < 1;
    #print join "***",@$all_ids_at_level;

    foreach my$all_ids(@$all_ids_at_level){
      foreach my$igene(@{$id_store{$all_ids}}){
        my$iname=$igene;
        if ($split_name) {$iname=lc(getname($igene,$split_name,$ncol))}
        foreach my$iexp(@{$exp_data{$iname}}){
          print $OUT2 "$all_ids\t$merc_data{$all_ids}{desc}\t$igene\t${$gene_desc_lev{$iname}}[-1]\t",join (":",@{$gene_desc_all{$iname}}),$iexp?"\t$iexp":"\tNoDesc","\n";
	  $printed++;
        }
      }
    }
  }
  close $OUT2;
  ### delete output file if no gene were found to print.
  if($printed < 1){print "\nNo gene found.Deleting $out_file\n";unlink $out_file}
}








print "\nSummary of DE genes in all files processed:";
foreach my$file(@files){
      print "\n$file : Total DE genes: $fnumgene{$file}";
}
print "\n\n";



#    ' wAnnot.RSEM_express.matrix_2_trimmed_reads_OnLAP.US56.combined_NormBy.TMM_Leaf.genes.counts.matrix.HighDWF2_vs_LowDWF2.DESeq2.DE_results.P0.05_C1.HighDWF2-UP.FC.mercator.results.txt wAnnot.RSEM_express.matrix_2_trimmed_reads_OnLAP.US56.combined_NormBy.TMM_Leaf.genes.counts.matrix.HighDWF2_vs_LowDWF2.DESeq2.DE_results.P0.05_C1.LowDWF2-UP.FC.mercator.results.txt wAnnot.RSEM_express.matrix_2_trimmed_reads_OnLAP.US56.combined_NormBy.TMM_Leaf.genes.counts.matrix.HighStalkVolF2_vs_LowStalkVolF2.DESeq2.DE_results.P0.05_C1.HighStalkVolF2-UP.FC.mercator.results.txt wAnnot.RSEM_express.matrix_2_trimmed_reads_OnLAP.US56.combined_NormBy.TMM_Leaf.genes.counts.matrix.HighStalkVolF2_vs_LowStalkVolF2.DESeq2.DE_results.P0.05_C1.LowStalkVolF2-UP.FC.mercator.results.txt


sub getname{
	my$name=shift;
	my$lim=shift;
	my$sncol=shift;
	$name=~s/^\s*//g;
	$name=~s/\s*$//g;
	my@frags=split /\Q$lim/,$name;
	unshift @frags, "";  ### push an empty block to match the column number to array index
	#print "\nreturned $frags[$sncol-1]  for  $name";
	##check if multiple columns were asked
	#my@mcols;
	#if($sncol=~/,/){
	my@mcols=split /,/,$sncol;
	#}else{@mcols=($sncol)}

	#if(@mcols > 1){$nname=join $lim,@frags[@mcols]}
	my$nname=join $lim,@frags[@mcols];
	print "\nOrg name: $name : ext name:$nname" if $verbose;
	return(join $lim,@frags[@mcols])
	#return ($frags[$sncol-1])

}

sub get_levels{
  my$ref_merc_data=shift;
  my$req_bid=shift;
  my$req_lev=shift;
  my$exclusive=shift;

  #print Dumper %$ref_merc_data;
  my@ibids; ## for collecting bin ids
  my@req_bin=split(/\./,$req_bid);
  #print Dumper $req_bid;
  $req_lev||=scalar@req_bin;

  foreach my$ibid(keys %{$ref_merc_data}){
    my@ibin=split(/\./,$ibid);
    #print Dumper $ibid;
    next if $#ibin < $#req_bin; ### exclude upper level bins.

    next if ( join (".",@req_bin) ne join (".",@ibin[0..$#req_bin]));
    #print Dumper $req_bid,join (".",@req_bin), $ibid,join (".",@ibin[0..$#req_bin]);
    #print Dumper join (".",@req_bin), join (".",@ibin), $#req_bin,$#ibin;
    if ($ibid =~ m/^\Q$req_bid/) {
      push(@ibids,$ibid) if $$ref_merc_data{$ibid}{lev} == $req_lev;
      print "saving $ibid as it matched $req_bid at level $req_lev\n" if $$ref_merc_data{$ibid}{lev} == $req_lev;
    }
  }

  #print Dumper $req_bid,$req_lev;

  return \@ibids;

}
