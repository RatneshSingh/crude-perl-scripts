use strict;
use warnings;
no autovivification;

    my%list;
    my%sum_grp;
    my%list_sp;
    my%seq;
    $/="\n";

    foreach my$file(@ARGV){
      open GRP,$file;
      print "\nOpening grp file: $file\n";
      my $grp;
      while(<GRP>){
        s/\s*$|^\s*|^>//g;
        next if /^\s*$/;
        next if /^Number\s+of\s+clusters/;
        $grp = $1 if (/^Cluster\s+([\w\d]+)\s*\;\s+size\=\d+/ || /^(unclustered) points/);
        next if /^Cluster\s+(\d+)\s\;\s+size\=\d+/ || /^(unclustered) points/;
        $list{$_}=$grp;
        my$sp=$_;


        $sp=~s/[_.]+[\S]*//g;

        $sp = "Zmat" if ($sp eq 'evm' || $sp eq 'novel' || $sp eq 'temp');
        $sum_grp{$file}{$grp}{$sp}+=1;
        $list_sp{$sp}=$sp;

        #print "\nAdded *$_* in group $grp";
      }
      close(GRP);
    }
    ### print summary of homolog count in each cluster.
    open OUT2,">","cluspack_cluster_summary.table";
    print OUT2 "File\tCluster\t",join("\t",sort keys %list_sp);
    foreach my$file(@ARGV){
      for my$group(sort keys %{$sum_grp{$file}}){
        print "\n**********Writing Summary for cluster:$file : $group\n";
        print OUT2 "\n$file\tCluster_$group";
        for my$sps(sort keys %list_sp){
          print OUT2 "\t",$sum_grp{$file}{$group}{$sps}?$sum_grp{$file}{$group}{$sps}:0;

        }

      }
    }
    close(OUT2);

    print "\n Summary were weritten in cluspack_cluster_summary.table\n\n";
