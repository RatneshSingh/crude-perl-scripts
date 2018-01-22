#!/usr/bin/perl 
use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;
our(%exp_data,$grp_uniq,$parent1,$parent2,@other_ind,$grp_cov,$min_read,$deresult_file,$grp_file,$rand_grp,$tmm_file,$log2foldc,$qval,$help,$verbose,$print_unique_alleles);

$qval=0.05;
$log2foldc=1;
$min_read=1;

GetOptions(
  "deresult|d=s"=>\$deresult_file,
  "grp|g=s"=>\$grp_file,
  "random_group|rg=i"=>\$rand_grp,
  "tmm|t=s"=>\$tmm_file,
  "uniq|u"=>\$grp_uniq,
  "poua"=>\$print_unique_alleles,
  "min_read|m=i"=>\$min_read,
  "par1|p1=s"=>\$parent1,
  "par2|p2=s"=>\$parent2,
  "oin=s"=>\@other_ind,
  "grp_cov|c=i"=>\$grp_cov,
  "foldc|f=i"=>\$log2foldc,
  "qval|q=f"=>\$qval,
  "help|h"=>\$help,
  "verbose|v"=>\$verbose
);

my$usage="

perl $0 -options

Options:
  deresult|d     DE result file with Allele specific fold change. Will print results from DE genes only.
               should have 4 columns: Gene:pos:SNP\tLog2FoldChnge\tpvalue\tqvalue 
  grp|g     grp_file    Group infor for each sample.
                2 col: Grp_name\tSample_Name
  random_group|rg create random group of individual.
            e.g. -rg 5 will create 2 groups of randomly selected 5 samples.
  tmm|t    tmm_file. TMM normalized expression matrix of SNP expression Matrix.
  uniq|u    print only alleles specific to one group only.
  poua       print only unique allele and ignore other alleles at the position.
            eg: if a position has A and T alleles, only print unique allele and skip the other.
  par1|p1   print the expression from parent1. eg. -p1 L-LA-Purple
  par2|p2   print the expression from parent2. eg. -p2 L-US56-14
  oin       print the expression from other individuals.
            use multiple times for more than one ind eg. -oin L-7013 -oin L-7012 ect...
  grp_cov|c Print only alleles which has group coverage for members.
            e.g. grp_cov 80 means only print if one of the grp has 80% members with the alleles present.
  min_read|m	minimum number of read group should have[1]
  foldc|f   minimum fold change to filter for in DE result file.
  qval|q    FDR adjusted pvalue limit for a expression to be included in DE result file.
  help|h    help. Print this help and quit.
  verbose|v print some debug info.
";


!$help or die "$usage";

die "Group file required with TMM file.\n$usage\n"  if ($tmm_file && !$grp_file);

die "Atleast one file (DE or TMM) is required.\n$usage\n"  if (!$tmm_file && !$deresult_file);
my(%selected_ase);
if ($deresult_file) {
   open (my$ASE,"<",$deresult_file) or die "\n\nUnable to open ASE file : $deresult_file\n\n";
  while (<$ASE>) {
    next if $. ==1;   ### skip first header
    
    ### 0:id    1:baseMeanA    2:BaseMeanB 3:baseMean  4:log2FC  5:log2FC_SE   6:stat 7:pval   8:padj
    s/^\s+//;
    my@elm=split /\s+/,$_;
    
    my$id=$elm[0];
    my$bmA=$elm[1];
    my$bmB=$elm[2];
    my$bm=$elm[3];
    my$lfc=$elm[4];
    my$lfc_se=$elm[5];
    my$stat=$elm[6];
    my$pval=$elm[7];
    my$padj=$elm[8];
    
    
    
    next if $padj > $qval;
    next if abs($lfc) < $log2foldc;
    
    my($iso,$pos,$base)=split /:/, $id;
    
    my$gene=$iso;
         #$gene=~s/_i\d+//;
 
    
    
    $exp_data{$gene}{$pos}{numallele}++;
    push(@{$exp_data{$gene}{$pos}{base}},$base);
    push(@{$exp_data{$gene}{$pos}{exp}},$lfc);
     
}

### print selected ASE data and fold change info
open(my$ASEO, ">", "$deresult_file\.ASE_q$qval.table") or die "\nUnable to open outfile\n";

foreach my$gene(keys %exp_data){
   foreach my$pos(keys %{$exp_data{$gene}}){
        if ($exp_data{$gene}{$pos}{numallele} > 1) {
            my@foldc=@{$exp_data{$gene}{$pos}{exp}};
            my@base=@{$exp_data{$gene}{$pos}{base}};
            my$avg=avg(@foldc);
            my$print_pos=0;
            foreach my$cexp(@foldc){
                $print_pos=1  if abs($avg-$cexp ) >= $log2foldc;
            }
            print  $ASEO "\n$gene\t$pos\t",join "\t",@foldc,@base  if $print_pos == 1;
            $selected_ase{$gene}{$pos}=1  if $print_pos == 1;
        }
        
        
   
   }
}
print "\nSaved ASE fold change results in $deresult_file\.ASE_q$qval.table\n";
close $ASEO;
}
###################################################################################################################
#exit;
### create a expression table for AS_DE genes from TMM normalized expression table arrange them in order of two groups
open my$grp,"<",$grp_file or die "\nUnable to open the group file $grp_file\n";
my%group_list;

while(<$grp>){
	s/^\s+//;s/\s+$//;
	next if m/^\s*$/;
	my($grp,$sample)=split /\s+/;
	push(@{$group_list{$grp}},$sample) if $sample !~ m/^\s*$/;
	
}

### include parents if author asks.
if ($parent1 || $parent2 || scalar@other_ind > 0) {
    push(@{$group_list{Parent1}},$parent1) if $parent1;
    push(@{$group_list{Parent2}},$parent2) if $parent2;
    push(@{$group_list{Other}},@other_ind) if scalar@other_ind>0;
}


print "\nAnalyzing group:",Dumper \%group_list if !$rand_grp;

#### create sample list in order of group.
my@sample_list;
my@sample_list_withg;
foreach my$grp(sort keys %group_list){foreach my$samp(sort @{$group_list{$grp}}){push(@sample_list,$samp);push(@sample_list_withg,"$grp.$samp");}}




#### read TMM normalized matrix for DE genes
if ($tmm_file) {
    no autovivification;
    my$tmm_out=join ".",$tmm_file,$grp_uniq?"group_uniq.minread$min_read":"exp",$grp_cov?"grpCov$grp_cov\.table":"table";
    open my$tmm,"<",$tmm_file or die "\nUnable to open the group file $tmm_file\n";
    #open(my$TMMO, ">", $tmm_out) or die "\nUnable to open outfile\n";
    my(@samp_names);my %tmm_exp;my %print_this;
    my @rand_group;
    #print $TMMO "Gene_AS\t",join "\t",@sample_list_withg;
    while(<$tmm>){
    	s/^\s+//;s/\s+$//;  ## remove spaces at the start and end
        #my$print_this=0;
    	### process expression data and store it per sample
        if ($.==1){
         @samp_names= split /\s+/ ;
         if($rand_grp){
            %group_list=(); ## empty hash if any previous grp names were in. The script cannot handle more than two groups at the moment.
            @sample_list_withg=();
            
            foreach my$grp("RandomGrp_1","RandomGrp_2"){
               my%hash=();
               my%hash_hd=();
            
              while ((scalar keys %hash) < $rand_grp) {
                  my $r = $samp_names[rand(@samp_names)];
                  $hash{$r} = $r;
                  $hash_hd{$r}=$grp.$r;
               }
              @{$group_list{$grp}}=(keys %hash);
              push(@sample_list_withg,@hash_hd{@{$group_list{$grp}}});
              
              
         }
            
         ### include parents if author asks.
         if ($parent1 || $parent2 || scalar@other_ind > 0) {
            push(@{$group_list{Parent1}},$parent1) if $parent1;
            push(@sample_list_withg,"Parent1.$parent1") if $parent1;
            push(@{$group_list{Parent2}},$parent2) if $parent2;
            push(@sample_list_withg,"Parent2.$parent2") if $parent2;
            if(scalar@other_ind > 0){
              foreach my $oth(@other_ind){
                 push(@{$group_list{Other}},$oth) if $oth;
                 push(@sample_list_withg,"Other.$oth") if $oth;
              }
            }
         }            
         
         
         print "\nAnalyzing Random groups:",Dumper \%group_list;
         
         
         }


         
         next;
        
        }  ## first line has header
        #process_tmm(\%tmm_exp,\@samp_names,$_);
    	my($genei,@exp)= split /\s+/;
        
        die "\n\nThere are more headers than avaibale data. \nPlease Make sure that there is only one header per data point. ID column does not need a header.\n" if !($#samp_names == $#exp);
        
    	my($gene,$pos,$base) = split/:/,$genei;
    
    	next if ($deresult_file && !exists $selected_ase{$gene}{$pos});  ###skip if gene is not listed in DE list and de result file is selected.
    
    	for(my$i=0;$i<=$#samp_names;$i++){
            $tmm_exp{$gene}{$pos}{$base}{$samp_names[$i]}=$exp[$i];   ### put expression values in hash per sample.
            $tmm_exp{$genei}{$samp_names[$i]}=$exp[$i];
            
        }
        ### process the data
        ##### print only lines which are group specific i.e. Present in one group but absent in other group.
        if ($grp_uniq) {
           my@grp_reads;
           my@grp_name;
           my$all_sum=0;
           my@per_nonzero;
           foreach my$grp(sort keys %group_list){
               next if $grp =~ m/^Parent\d/; ## do not include parent in calculations.
               next if $grp =~ m/^Other\d/; ## do not include Othert in calculations.
               my $grp_sum=0;
               my $num_nonzero=0;
               my $num_zero=0;
               
               foreach my$samp(sort @{$group_list{$grp}}){
                   $grp_sum+=$tmm_exp{$genei}{$samp} ;
                   $all_sum+=$tmm_exp{$genei}{$samp} ;
                   if($tmm_exp{$genei}{$samp}>0){$num_nonzero+=1 } else{$num_zero+=1 ;}
               }
               push(@per_nonzero,$num_nonzero*100/($num_nonzero+$num_zero)) if $grp_cov;
               push(@grp_reads,$grp_sum);
               push(@grp_name,$grp);
           }
           print Dumper '@per_nonzero',\@per_nonzero, '@grp_reads',\@grp_reads,'@grp_name', \@grp_name if $verbose;
           ## test filter and mark for printing.
           
           if ($all_sum >= $min_read && ($grp_reads[0] == 0 || $grp_reads[1] ==0)) {
               #$print_this=1;
               $print_this{$grp_name[1]}{$gene}{$pos}{$base}=1 if $grp_reads[0] == 0;
               $print_this{$grp_name[0]}{$gene}{$pos}{$base}=1 if $grp_reads[1] == 0;
               
               if ($grp_cov && $per_nonzero[0] < $grp_cov && $per_nonzero[1] < $grp_cov){
                     delete $print_this{$grp_name[1]}{$gene}{$pos}{$base} if $grp_reads[0] == 0;
                     delete $print_this{$grp_name[0]}{$gene}{$pos}{$base} if $grp_reads[1] == 0;
                }
            }
           
            #if ($print_this == 1){print $TMMO "\n$genei";foreach my$samp(@sample_list){print  $TMMO  "\t$tmm_exp{$genei}{$samp}"};}
        
        } else {   
    
            ### print expression data organized as group
            #print $TMMO "\n$genei";foreach my$samp(@sample_list){print  $TMMO  "\t$tmm_exp{$genei}{$samp}"}
            $print_this{All}{$gene}{$pos}{$base}=1;
        }
    }

    ## print selected genes and positions.
    foreach my$pgrp(keys %print_this){
        #next if $pgrp =~ m/^\s*$/;
       my $tmm_out2=join ".",$tmm_out,$pgrp,$print_unique_alleles?"uniqOnly.table":"table";
       open(my$TMMO,">",$tmm_out2) or die "\nUnable to create file:$tmm_out2 to save expression info\n";
       print $TMMO "Gene_AS\t",join "\t",@sample_list_withg;
       foreach my$pgene(sort keys %{$print_this{$pgrp}}){
           foreach my $ppos(sort keys %{$print_this{$pgrp}{$pgene}}){
               next if scalar (keys %{$print_this{$pgrp}{$pgene}{$ppos}}) < 1;
               
               if($print_unique_alleles){
                  foreach my $pbase(sort keys %{$print_this{$pgrp}{$pgene}{$ppos}}){
                    
                       	print $TMMO "\n$pgene:$ppos:$pbase";
                       	foreach my$samp(@sample_list){print  $TMMO  "\t$tmm_exp{$pgene}{$ppos}{$pbase}{$samp}"};
                  }    
               }
               else{
                  foreach my $pbase(sort keys %{$tmm_exp{$pgene}{$ppos}}){
                       	print $TMMO "\n$pgene:$ppos:$pbase";
                       	foreach my$samp(@sample_list){print  $TMMO  "\t$tmm_exp{$pgene}{$ppos}{$pbase}{$samp}"};
                       
                    }
		    
               }
            }
        }
       print "\nSaved selected TMM results in: $tmm_out2\n";
       close $TMMO;
    }
    
    #print "\nSaved selected TMM results in: $tmm_out \n";
    
    #print Dumper(%print_this);
}





##############################################################################################
##### SUBROUTINES
##############################################################################################

sub avg{ return sum(@_)/scalar@_}
sub sum{my$sum=0; $sum+=$_ foreach @_;return $sum}


sub process_tmm{
    my $rH_tmm_exp=shift;
    my $rA_header=shift;
    my $line=shift;
    
    
	$line=~s/^\s+//;$line=~s/\s+$//;  ## remove spaces at the start and end

	### process expression data and store it per sample 
	#if ($.==1){@samp_names= split /\s+/,$line ; next}  ## first line has header
	my($genei,@exp)= split /\s+/,$line;
    
    die "\n\nThere are more headers than avaibale data. \nPlease Make sure that there is only one header per data point. ID column does not need a header.\n" if !($#$rA_header == $#exp);
    
	my($gene,$pos,$base) = split/:/,$genei;

	#next if ($deresult_file && !exists $selected_ase{$gene}{$pos});  ###skip if gene is not listed in DE list and de result file is selected.

	for(my$i=0;$i<=$#$rA_header;$i++){$rH_tmm_exp->{$gene}->{$pos}->{$base}->{$$rA_header[$i]}=$exp[$i]}   ### put expression values in hash per sample.
	
}

