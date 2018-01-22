use strict;
use warnings;
no autovivification;
    open GRP,$ARGV[0];
    my%list;
    my%sum_grp;
    my%list_sp;
    my%seq;
  print "\nOpening grp file: $ARGV[0]\n";
    $/="\n";
    my $grp=0;
    while(<GRP>){
      s/\s*$|^\s*|^>//g;
      next if /^\s*$/;
      next if /^Number\s+of\s+clusters/;
      $grp = $1 if (/^Cluster\s+([\w\d]+)\s*\;\s+size\=\d+/ || /^(unclustered) points/);
      next if /^Cluster\s+(\d+)\s\;\s+size\=\d+/ || /^(unclustered) points/;
      $list{$_}=$grp;
      my$sp=$_;
      $sp=~s/[_.]+[\S]*//g;
      $sum_grp{$grp}{$sp}+=1;
      $list_sp{$sp}=$sp;

      #print "\nAdded *$_* in group $grp";
    }




    close(GRP);
    print "\nOpening seq file: $ARGV[1]\n";
    $/="\n>";
    open(SEQ,$ARGV[1]);
    while(<SEQ>){
      my($head,@seq)=split(/\n/,$_);
      $head=~s/\s*$|^\s*|^>//g;
      #print "\nProcessing seq *$head*";
      #next if ! exists $list{$head};
      my$seq=join("",@seq);
      $seq=~s/[^\w\-]+//g;
      $seq{$list{$head}}{$head}=$seq;

    }



    for my$group(sort keys %seq){
      print "\n**********Writing seqs for: $group\n";
      open OUT,">",$ARGV[1]."_clust_".$group.".fasta";
      for my$heads(keys %{$seq{$group}}){
        print OUT ">$heads\n$seq{$group}{$heads}\n";
        #print "\n$heads  for group $group\n";
      }
      close(OUT);
    }

    ### print summary of homolog count in each cluster.
    open OUT2,">",$ARGV[1]."_clust_summary.table";
    for my$group(sort keys %seq){
      print "\n**********Writing Summary for cluster: $group\n";
      print OUT2 "\nCluster_$group";
      for my$sps(sort keys %list_sp){
        print OUT2 "\t$sps\:",$sum_grp{$group}{$sps}?$sum_grp{$group}{$sps}:0;

      }

    }
    close(OUT2);
