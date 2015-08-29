#! /usr/bin/perl -w
use strict;
use Getopt::Std;

our($opt_i,$opt_o,$opt_p,$opt_f,$exten,$opt_e);
$opt_p='yes';
$opt_f='genscan';


getopt('i:o:p:f:e');

#$opt_i="genscan_No_name57221.txt";

my$usage="perl Script -options

-i  Genscan/fgenesh OutFile [Required]
-o  Outputfile_name [infile.gff]
-f  genscan or fgenesh [genscan]
-p  yes|no. Parse partial genes[yes]
-e  1|0 1:print extended features. e.g PolyA and TF binding Sites stc[0].

\n";




die "\nUsage: $usage\n\n\n" if !$opt_i;
if(lc $opt_f =~ /genscan/i){open  GENSCAN,"$opt_i" or die "Cannot open file $opt_i\n$usage\n\n\n";}
elsif(lc $opt_f =~/fgenesh/i){open  FGENESH,"$opt_i" or die "Cannot open file $opt_i\n$usage\n\n\n";}
my($sequence_name,$sequence_length,$program,%gff,$fh,%numgenes);

$program='Genscan' if $opt_f=~/genscan/i;
$program='Fgenesh' if $opt_f=~/fgenesh/i;
if (lc $opt_f =~ /genscan/i){

$exten="GS";

while (<GENSCAN>){
next if /^\s*$|^#/;
#die "No exon or gene found in the file" if /'NO EXONS\/GENES PREDICTED IN SEQUENCE'/;
if(/^Sequence\s+(\S+)\s*:\s*(\d+)\s*bp/){$sequence_name=$1;$sequence_length=$2}    # print "\n*-->Processing Sequence $sequence_name";
if(/^\s*([\d\.]+)\s+(\w+)/){

    $_=~s/^\s+//;
    #line has following elements
    # 0:Gn.Ex
    # 1:Type
    # 2:S .
    # 3:Begin ...
    # 4:End .
    # 5:Len
    # 6:Fr
    # 7:Ph
    # 8:I/Ac
    # 9:Do/T
    # 10:CodRg
    # 11:P....
    # 12:Tscr..
    my@line_elements=split(/\s+/,$_);
    my($gene_num,$feature_num)=split(/\./,$line_elements[0]);
    $numgenes{$sequence_name}{$gene_num}=$gene_num; # count the number of genes each sequence has.

    $gff{$sequence_name}{$gene_num}{'strand'}=$line_elements[2];
    if($line_elements[2] eq '+'){ # strand if positive

                # get the start and end of mRNA and gene
                $gff{$sequence_name}{$gene_num}{'mRNA_start'}=$line_elements[3] if $line_elements[1] eq 'Init';
                $gff{$sequence_name}{$gene_num}{'mRNA_start'}=$line_elements[3] if $line_elements[1] eq 'Sngl';
                $gff{$sequence_name}{$gene_num}{'mRNA_end'}=$line_elements[4] if $line_elements[1] eq 'Term';
                $gff{$sequence_name}{$gene_num}{'mRNA_end'}=$line_elements[4] if $line_elements[1] eq 'Sngl';
                #if($line_elements[1] ne 'Prom' && $line_elements[1] ne 'PlyA'){
                if($line_elements[1] =~ /(Intr|Term|Init|Sngl)/){

                    my@cds_coords=($line_elements[3],$line_elements[4]);
                    my$cds_coords=$line_elements[3].'.'.$line_elements[4];
                    $gff{$sequence_name}{$gene_num}{$cds_coords}=$line_elements[7];
                    push(@{$gff{$sequence_name}{$gene_num}{'Exons'}},@cds_coords);
                    #print "CDS-coords Plus strand:@cds_coords\n";
                }
                elsif($line_elements[1] =~/PlyA/){
                  $gff{$sequence_name}{$gene_num}{'PolyAStart'}=$line_elements[3];
                  $gff{$sequence_name}{$gene_num}{'PolyAEnd'}=$line_elements[4];
                }
                elsif($line_elements[1] =~/Prom/){

                  $gff{$sequence_name}{$gene_num}{'TSSStart'}=$line_elements[3];
                  $gff{$sequence_name}{$gene_num}{'TSSEnd'}=$line_elements[4];

                }

   }
    else{ # strand is negative

                # get the start and end of mRNA and gene
                $gff{$sequence_name}{$gene_num}{'mRNA_start'}=$line_elements[4] if $line_elements[1] eq 'Term';
                $gff{$sequence_name}{$gene_num}{'mRNA_start'}=$line_elements[4] if $line_elements[1] eq 'Sngl';
                $gff{$sequence_name}{$gene_num}{'mRNA_end'}=$line_elements[3] if $line_elements[1] eq 'Init';
                $gff{$sequence_name}{$gene_num}{'mRNA_end'}=$line_elements[3] if $line_elements[1] eq 'Sngl';

                if($line_elements[1] =~ /(Intr|Term|Init|Sngl)/){

                    my@cds_coords=($line_elements[4],$line_elements[3]);
                    my$cds_coords=min($line_elements[4],$line_elements[3]).max($line_elements[4],$line_elements[3]);
                    $gff{$sequence_name}{$gene_num}{$cds_coords}=$line_elements[7];
                    push(@{$gff{$sequence_name}{$gene_num}{'Exons'}},@cds_coords);
                    #print "CDS-coords Negative strand:@cds_coords\n";
                }
               elsif($line_elements[1] =~/PlyA/){
                  $gff{$sequence_name}{$gene_num}{'PolyAStart'}=$line_elements[3];
                  $gff{$sequence_name}{$gene_num}{'PolyAEnd'}=$line_elements[4];
               }
               elsif($line_elements[1] =~/Prom/){

                  $gff{$sequence_name}{$gene_num}{'TSSStart'}=$line_elements[3];
                  $gff{$sequence_name}{$gene_num}{'TSSEnd'}=$line_elements[4];

         }

    }
  }
}
}
elsif(lc$opt_f=~/fgenesh/i){

$exten="FG";
######################################################
##### Parsing fgenesh output
######################################################




while (<FGENESH>){
s/^\s+//g;
next if /^\s*$|^#|^FGENESH|^Time|^Number|^Positions|^G\s+Str/;  # add '|\s+TSS\s+|\s+PolA\s+' to exclude promoter and polyA
last if /^Predicted protein/;
#die "No exon or gene found in the file" if /'NO EXONS\/GENES PREDICTED IN SEQUENCE'/;
if(/^Seq\s+name\:\s*(\S+)\s*/){$sequence_name=$1;}    # print "\n*-->Processing Sequence $sequence_name";
if(/^Length\s*of\s*sequence\s*\:\s*(\d+)\s*$/){$sequence_length=$1}
if(/^([\d]+)\s+([\+\-]\s+)/){


    #line has following elements
    # 0:GeneNumber
    # 1:Strand
    # 2:FeatureNum.
    # 3:FeatureType
    # 4:FeatureStart
    # 5:-
    # 6:FeatureEnd
    # 7:Score
    # 8:codon_start_in_Exon
    # 9:-
    # 10:Codon_end_in_exon
    # 11:coding_length
    # 12:

    my@line_elements=split(/\s+/,$_);
    my($gene_num,$feature_num)=($line_elements[0],$line_elements[2]);
    $numgenes{$sequence_name}{$gene_num}=$gene_num; # count the number of genes each sequence has.

    $gff{$sequence_name}{$gene_num}{'strand'}=$line_elements[1];
    if($line_elements[1] eq '+'){ # strand if positive
         if($line_elements[3] =~ /CDSf|CDSo|CDSl/){
                # get the start and end of mRNA and gene
                $gff{$sequence_name}{$gene_num}{'mRNA_start'}=$line_elements[4] if $line_elements[3] eq 'CDSf';
                $gff{$sequence_name}{$gene_num}{'mRNA_start'}=$line_elements[4] if $line_elements[3] eq 'CDSo';
                $gff{$sequence_name}{$gene_num}{'mRNA_end'}=$line_elements[6] if $line_elements[3] eq 'CDSl';
                $gff{$sequence_name}{$gene_num}{'mRNA_end'}=$line_elements[6] if $line_elements[3] eq 'CDSo';


                    my@cds_coords=($line_elements[4],$line_elements[6]);
                    my$cds_coords=$line_elements[4].'.'.$line_elements[6];
                    $gff{$sequence_name}{$gene_num}{$cds_coords}=$line_elements[8]-$line_elements[4];
                    push(@{$gff{$sequence_name}{$gene_num}{'Exons'}},@cds_coords);
                    #print "CDS-coords Plus strand:@cds_coords\n";

         }
         elsif($line_elements[2] =~/PolA/){
                  $gff{$sequence_name}{$gene_num}{'PolyAStart'}=$line_elements[3];
                  $gff{$sequence_name}{$gene_num}{'PolyAEnd'}=$line_elements[3]+6;
         }
         elsif($line_elements[2] =~/TSS/){

                  $gff{$sequence_name}{$gene_num}{'TSSStart'}=$line_elements[3];
                  $gff{$sequence_name}{$gene_num}{'TSSEnd'}=$line_elements[3]+12;

         }

   }
    elsif($line_elements[1] eq '-'){ # strand is negative

                if($line_elements[3] =~ /CDSl|CDSo|CDSf/){
                    # get the start and end of mRNA and gene
                    $gff{$sequence_name}{$gene_num}{'mRNA_start'}=$line_elements[4] if $line_elements[3] eq 'CDSl';
                    $gff{$sequence_name}{$gene_num}{'mRNA_start'}=$line_elements[4] if $line_elements[3] eq 'CDSo';
                    $gff{$sequence_name}{$gene_num}{'mRNA_end'}=$line_elements[6] if $line_elements[3] eq 'CDSf';
                    $gff{$sequence_name}{$gene_num}{'mRNA_end'}=$line_elements[6] if $line_elements[3] eq 'CDSo';

                    my@cds_coords=($line_elements[4],$line_elements[6]);
                    my$cds_coords=min($line_elements[4],$line_elements[6]).max($line_elements[4],$line_elements[6]);
                    $gff{$sequence_name}{$gene_num}{$cds_coords}=$line_elements[6]-$line_elements[10];
                    push(@{$gff{$sequence_name}{$gene_num}{'Exons'}},@cds_coords);
                    #print "CDS-coords Negative strand:@cds_coords\n";
                }
                elsif($line_elements[2] =~ /PolA/){
                  $gff{$sequence_name}{$gene_num}{'PolyAStart'}=$line_elements[3]-6;
                  $gff{$sequence_name}{$gene_num}{'PolyAEnd'}=$line_elements[3];
                }
                elsif($line_elements[2] =~ /TSS/){

                  $gff{$sequence_name}{$gene_num}{'TSSStart'}=$line_elements[3]-12;
                  $gff{$sequence_name}{$gene_num}{'TSSEnd'}=$line_elements[3];

                }
    }
  }
}

}
else{ die "Please provide proper format of input file. e.g. genscan or fgenesh\n\n$usage\n\n"; }





###########################################################
#### END fgenesh parsing
###########################################################








my$outfile=$opt_i.'.gff' if !$opt_o;
$outfile=$opt_o if $opt_o;;
open OUT,">$outfile";
print OUT "\#\#gff-version 3\n";
print OUT "\#\#sequence-region $sequence_name 1 $sequence_length\n";
# collect and sort the number of genes each sequence have.
 foreach my$seqname(keys %numgenes){
    foreach my$genenum(keys %{$numgenes{$seqname}}){
        push(@{$numgenes{$seqname}{'gene_num'}},$genenum)
    }
    if(scalar@{$numgenes{$seqname}{'gene_num'}} >=1 ){
        @{$numgenes{$seqname}{'gene_num'}} = sort{$a<=>$b} @{$numgenes{$seqname}{'gene_num'}};
    }
    else{undef @{$numgenes{$seqname}{'gene_num'}};}

}

 #@genenumarray=sort{$a<=>$b}@genenumarray;



 foreach my$seqname(keys %gff){
    #print "**********************Reading for Sequence:$seqname**************************\n";
    next if !@{$numgenes{$seqname}{'gene_num'}};

    foreach my$genenum(@{$numgenes{$seqname}{'gene_num'}}){

        next if (lc$opt_p eq 'no' && (!$gff{$seqname}{$genenum}{'mRNA_start'} || !$gff{$seqname}{$genenum}{'mRNA_end'})); #in partial genes mRNA_start and mRNA_end is not defined, will throw error.
        #next if (!$gff{$seqname}{$genenum}{'mRNA_start'} || !$gff{$seqname}{$genenum}{'mRNA_end'});
        #print "\n-->Processing Gene:$seqname\tNum: $genenum";
        my $partial='Partial';
        if ($gff{$seqname}{$genenum}{'mRNA_start'} && $gff{$seqname}{$genenum}{'mRNA_end'}){$partial=""}
        if ($opt_e){ print OUT "$seqname\t$program\tgene\t".min(@{$gff{$seqname}{$genenum}{'Exons'}},$gff{$seqname}{$genenum}{'TSSStart'},$gff{$seqname}{$genenum}{'TSSEnd'},$gff{$seqname}{$genenum}{'PolyAEnd'},$gff{$seqname}{$genenum}{'PolyAStart'})."\t".max(@{$gff{$seqname}{$genenum}{'Exons'}},$gff{$seqname}{$genenum}{'TSSStart'},$gff{$seqname}{$genenum}{'TSSEnd'},$gff{$seqname}{$genenum}{'PolyAEnd'},$gff{$seqname}{$genenum}{'PolyAStart'})."\t.\t$gff{$seqname}{$genenum}{'strand'}\t.\tID=$seqname.${partial}Gene.$genenum-${exten}\n" }
            else{print OUT "$seqname\t$program\tgene\t".min(@{$gff{$seqname}{$genenum}{'Exons'}})."\t".max(@{$gff{$seqname}{$genenum}{'Exons'}})."\t.\t$gff{$seqname}{$genenum}{'strand'}\t.\tID=$seqname.${partial}Gene.$genenum-${exten}\n" ;}

            print OUT "$seqname\t$program\tTF_binding_site\t$gff{$seqname}{$genenum}{'TSSStart'}\t$gff{$seqname}{$genenum}{'TSSEnd'}\t.\t$gff{$seqname}{$genenum}{'strand'}\t.\tParent=$seqname.${partial}Gene.$genenum-${exten}\n" if ($gff{$seqname}{$genenum}{'TSSStart'} && $gff{$seqname}{$genenum}{'strand'} eq '+' && $opt_e);

            print OUT "$seqname\t$program\tpolyA\t$gff{$seqname}{$genenum}{'PolyAStart'}\t$gff{$seqname}{$genenum}{'PolyAEnd'}\t.\t$gff{$seqname}{$genenum}{'strand'}\t.\tParent=$seqname.${partial}Gene.$genenum-${exten}\n" if ($gff{$seqname}{$genenum}{'PolyAStart'} && $gff{$seqname}{$genenum}{'strand'} eq '-' && $opt_e);

            print OUT "$seqname\t$program\tmRNA\t".min(@{$gff{$seqname}{$genenum}{'Exons'}})."\t".max(@{$gff{$seqname}{$genenum}{'Exons'}})."\t.\t$gff{$seqname}{$genenum}{'strand'}\t.\tID=$seqname.${partial}mRNA.$genenum-${exten}\;Parent=$seqname.${partial}Gene.$genenum-${exten}\n";

#        }

       # print "\n\n.................Reading for Gene Number:$genenum...............\n";
       # my$num_lines=scalar@{$gff{$seqname}{$genenum}{'Exons'}};
        #print "Start of elements in Sequence:$seqname\tGene number:$genenum\t Has $num_lines lines\n";

        my@cds_array=@{$gff{$seqname}{$genenum}{'Exons'}};
        @cds_array=sort{$a<=>$b}@cds_array;
        my$cds_len=0;
        for(my$i=0;$i<=scalar@cds_array -1;$i=$i+2){

            my$cds_phase;

            if($i==0){$cds_phase=0}
            else{
                for(my$j=0;$j<$i;$j=$j+2){$cds_len+=(($cds_array[$j+1]-$cds_array[$j])+1);}
                $cds_phase=0 if $cds_len%3 == 0;
                $cds_phase=1 if $cds_len%3 == 2;
                $cds_phase=2 if $cds_len%3 == 1;
            }

           # my$phase=$gff{$seqname}{$genenum}{$cds_array[$i].$cds_array[$i+1]};
            #$phase=0 if !$gff{$seqname}{$genenum}{$cds_array[$i].$cds_array[$i+1]};
            print OUT "$seqname\t$program\tCDS\t$cds_array[$i]\t$cds_array[$i+1]\t.\t$gff{$seqname}{$genenum}{'strand'}\t$cds_phase\tID=$seqname.cds.$genenum-${exten}\;Parent=$seqname.${partial}mRNA.$genenum-${exten}\n";

        }

            print OUT "$seqname\t$program\tTF_binding_site\t$gff{$seqname}{$genenum}{'TSSStart'}\t$gff{$seqname}{$genenum}{'TSSEnd'}\t.\t$gff{$seqname}{$genenum}{'strand'}\t.\tParent=$seqname.${partial}Gene.$genenum-${exten}\n" if ($gff{$seqname}{$genenum}{'TSSStart'} && $gff{$seqname}{$genenum}{'strand'} eq '-' && $opt_e);
            print OUT "$seqname\t$program\tpolyA\t$gff{$seqname}{$genenum}{'PolyAStart'}\t$gff{$seqname}{$genenum}{'PolyAEnd'}\t.\t$gff{$seqname}{$genenum}{'strand'}\t.\tParent=$seqname.${partial}Gene.$genenum-${exten}\n" if ($gff{$seqname}{$genenum}{'PolyAStart'} && $gff{$seqname}{$genenum}{'strand'} eq '+' && $opt_e);




        #print "@{$gff{$seqname}{$genenum}}\n End of Sequence:$seqname\t gene:$genenum\n"



    }
 }


sub min{
  no warnings 'uninitialized';
  @_=grep(!/^$/, @_);
  @_= sort { $a <=> $b }@_;
  return $_[0]


}

sub max{
  no warnings 'uninitialized';
  @_=grep(!/^$/, @_);
  @_ = sort { $a <=> $b }@_;
  return $_[-1];


}