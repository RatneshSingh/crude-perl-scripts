
# count the number of nucleotides in a sequence file.
perl -e 'foreach(<>){$N=($_=~tr/ATGCatgc/ATGCatgc/); $Tn+=$N};print "$Tn\n";'

  #To count nucleotides (including Ns) in a sequence file
  perl -e 'foreach(<>){$N=($_=~tr/ATGCatgc/ATGCatgc/); $Tn+=$N};print "$Tn\n";'

  #To count nucleotides (excluding headers, including Ns) in a sequence file
  perl -e 'foreach(<>){next if />/; $N=($_=~tr/ATGCNatgcn/ATGCNatgcn/); $Tn+=$N};print "$Tn\n";'

## count number of nuceotides in each sequence in fasta file
perl  -e '$/="\n>"; while(<>){($header,@sequence)=split(/\n/,$_); $header=~s/>|^\s+//g; $sequence=join("",@sequence); $sequence=~s/\s|>//g; $sequence{$header}=$sequence;} foreach(keys %sequence){print "$_\t".length($sequence{$_})."\n"} '

for file in "Sitalica_164_v2.hardmasked.fa" "Sbicolor_v2.1_255.hardmasked.fa" "Pvirgatum_273_v1.0.hardmasked.fa" "ZMgenome.chromosomes.hardmasked.fa"; do perl  -e '$/="\n>"; while(<>){($header,@sequence)=split(/\n/,$_); $header=
~s/>|^\s+//g; $sequence=join("",@sequence); $sequence=~s/\s|>//g; $sequence{$header}=$sequence;} foreach(keys %sequence){print "$_\t".length($sequence{$_})."\n"} ' < $file>$file.seqlen.txt; done



#To count nucleotides (excluding headers, excluding Ns) in a sequence file
  perl -e 'foreach(<>){next if />/; $N=($_=~tr/ATGCatgc/ATGCatgc/); $Tn+=$N};print "$Tn\n";'


## extract a sequence from multifasta file based on pattern match.
perl -e '$start=0; while(<>){$start=0 if $_=~/^>/;$start=1 if $_=~/Ming1n2n3n4n8kb_contig_228773\D*/; print $_ if $start==1;}'  Ming1n2n3n4n8kb_CLCassembly_contigs.gt1k.fa


  #create hash from fasta sequences
  perl -e '$/="\n>"; while(<>){($header,@sequence)=split(/\n/,$_); $header=~s/>|^\s+//g; $sequence=join("",@sequence); $sequence=~s/\s|>//g; $sequence{$header}=$sequence;}  '

  #### Split multifasta to individual files for each sequence
  perl -e '$/="\n>"; while(<>){($header,@sequence)=split(/\n/,$_); $header=~s/>|^\s+//g; $sequence=join("",@sequence); $sequence=~s/\s|>//g; $shead=$header;$shead=~s/\s+.*$//g;open OUT, ">$shead.fa";  print OUT ">$header\n$sequence\n"; close OUT;} '

#### convert one liner fasta to multiline fasta
  perl
  -e '$/="\n>"; while(<>){($header,@sequence)=split(/\n/,$_); $header=~s/>|^\s+//g; $sequence=join("",@sequence); (/(.{100})/,$sequence);$mseq=join("\n",@seqn);$mseq=~s/\n+/\n/g;  print "\n>$header$mseq"  }  '

#### break sequences on Multi N's and create new sequence with increasing number in heder_count++. Useful for splitting assembled sequences joined by stretch of NNNNNs.
  perl -ne
  '$/="\n>"; ($head,@seq)=split /\n/,$_; $seq=join("",@seq);$seq=~s/>//g; $head=~s/>//;$count=0;foreach $seq_s(split(/N+/,$seq)){print "\n>$head\_$count\n$seq_s";$count++}'

  #create hash from fasta sequences and assemble forward and reverse sequences with CAP3
  perl
  -e '$/="\n>"; while(<>){($header,@sequence)=split(/\n/,$_); $header=~s/>|^\s+//g; $sequence=join("",@sequence); $sequence=~s/\s|>//g; $sequence{$header}=$sequence; @uqh=split(/\s+/,$header);$uniq=$uqh[0];@uqh=();$char=chop($uniq);push(@{$uniqheaders{$uniq}},$header);} foreach(keys %uniqheaders){ system ("rm $_.FnR.fasta"); open OUT,">>$_.FnR.fasta"; foreach $head(@{$uniqheaders{$_}}){print OUT ">$head\n$sequence{$head}\n"} system ("cap3 $_.FnR.fasta")}  '

  #create hash from fasta sequences and assemble forward and reverse sequences with CAP3 blastn and trim vector
  perl
  -e '$/="\n>"; while(<>){($header,@sequence)=split(/\n/,$_); $header=~s/>|^\s+//g; $sequence=join("",@sequence); $sequence=~s/\s|>//g; $sequence{$header}=$sequence; @uqh=split(/\s+/,$header);$uniq=$uqh[0];@uqh=();$char=chop($uniq);push(@{$uniqheaders{$uniq}},$header);} foreach(keys %uniqheaders){ system ("rm $_.FnR.fasta"); open OUT,">>$_.FnR.fasta"; foreach $head(@{$uniqheaders{$_}}){print OUT ">$head\n$sequence{$head}\n"} system ("cap3 $_.FnR.fasta"); system ("rm $_.FnR.fasta.cap.ace  $_.FnR.fasta.cap.contigs.links  $_.FnR.fasta.cap.contigs.qual $_.FnR.fasta.cap.info $_.FnR.fasta.cap.singlets "); system ("blastn -query $_.FnR.fasta.cap.contigs  -db pGEMT.fasta -out $_.FnR.fasta.cap.contigs.vs.pGEMT.blastn"); system ("perl ~/Scripts/Perl/Blast_parser_BIOPERL_V1.0.pl -i $_.FnR.fasta.cap.contigs.vs.pGEMT.blastn -o $_.FnR.fasta.cap.contigs.vs.pGEMT.blastn.table"); system ("perl ~/Scripts/Perl/Blast_Table_Sumup_Redundant_hits_V1.2.pl -b $_.FnR.fasta.cap.contigs.vs.pGEMT.blastn.table -vt -s $_.FnR.fasta.cap.contigs -o $_.FnR.fasta.cap.contigs.vectrimmed.fasta -d \" \" -c 1 -k query") }  '
  < SEQ05- 25 - 11. fasta

# to count the number of reads and coverage from newbler assembled sequences. The headers have the information of reads and assembled length. Average length of reads was considered 500bp.
  perl -ne
  'next if $_!~/^>/;chomp $_;$numread=$length=$coverage=1;$numread=$1 if $_=~/numreads\=(\d+)/;$length=$1 if $_=~/length\=(\d+)/; $coverage=500*$numread/$length;print "$_\t$length\t$numread\t$coverage\n"'
  $file

# create summary of percent identity binning from blast results. Blast results should be parsed through "Blast_parser_BIOPERL_V1.0.pl" and then through "Blast_Table_Sumup_Redundant_hits_V1.3.pl"
  perl
  -e ' while(<>){$data=$1 if $_=~/\%of TotalLength([\s\S]+)/; $frag=$1 if $_=~/Total number of hits passed the filter criteria:(\d+)/; $mb=$2 if $_=~/Total Length of sequences passed the filter criteria:([\d\.]+) kb \/([\d\.]+)Mb/;} print "$frag\t$mb Mb\t$data\n"'
  file
  . blastn
  . table
  . binned

##### pick all gff lines for a named gene "AT5G47970" in this example. replace "AT5G47970" withgene name you want.
  perl -ne 'if(/AT5G47970/){push(@ids,$1) if /ID=([^;]+)/}; $ids=join("|",@ids); print  if ($ids ne ""&& /$ids/i);'

##### pick gene names between given coordinates on a chromosome from gff file.
  perl -ne '$chr="Chr10"; $start=28242067;$end=42201544; next if $_ !~m/$chr/;next if $_ !~m/gene/ @gff=split /\s+/; print if($gff[3] >= $start && $gff[3] <= $end )'

  printf "Chromosome\tGapStart\tGapEnd\tGapLength\tNumGenes\n"&&
while read chr str end ext; do
  let gap=end-str;
  printf "$chr\t$str\t$end\t$gap\t";
  perl -e '
  open GFF,"$ARGV[0]";
  $chr=$ARGV[1];
  $start=$ARGV[2];
  $end=$ARGV[3];
  while (<GFF>){
  next if $_ !~ m/$chr/;
  next if $_ !~m/gene/;
  @gff=split /\s+/;
  print if($gff[3] >= $start && $gff[3] <= $end )
  }'  Sbicolor_v2.1_255_gene.gff3 $chr $str $end |wc -l; done  < temp_table.txt


### identify the location of longest stretch of string in sequence

perl -e '
use warnings;
use strict;


my%sequence;
my@chr=("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10");
#my $len = 5000000;
$/="\n>";
while(<>){
  (my$header,my@sequence)=split(/\n/,$_);
  $header=~s/>|^\s+//g;
  my$sequence=join("",@sequence);
  $sequence=~s/\s|>//g;
  $sequence{$header}=$sequence;
  #push(@chr,$header)
}
foreach my$chr(@chr){
  my $str=$sequence{$chr};
  my $sstr = substr ($str, 0, 1) . (substr ($str, 1) ^ $str);
  my @bestRuns;
  my $match = "\0";
  my $bestRunLen = 2;
  my $scan = 0;
  while (-1 != (my $start = index $sstr, $match, $scan)) {
      my $runLen = length ((substr ($sstr, $start) =~ /(\0+)/)[0]) + 1;

      if ($runLen > $bestRunLen) {
          # new best match
          @bestRuns = $start - 1;
          $bestRunLen = $runLen;
          $match = "\0" x ($bestRunLen - 1);
      } else {
          # another best match
          push @bestRuns, $start - 1;
      }

      $scan = $start + $bestRunLen - 1;
  }

  for (@bestRuns) {
      print "$chr\t" . substr ($str, $_, 1) . "\t$_\t" . ($_ + $bestRunLen - 1) ."\t$bestRunLen\n";
  }
  }

' < Sbicolor_255_v2.0.softmasked.fa




















##### Parse denovo stacks log to identify stacks and loci used.
  grep -A 14 \
  "/usr/local/bin//sstacks" \
  denovo_map.log| grep -v "Parsing" | grep -B 1 -A 5 "stacks compared against the catalog containing" > match_summary.log

  printf "\nSample_Name\tNumberOfStacks\tNumberOfLoci\tNumberOfHaplotypes\n";
perl -e '
        $/="--";
        $num_loci=$samp=$hap=$stacks=0;
        while(<>){
        $num_loci=$samp=$hap=$stacks=0;
        $num_loci=$1 if(/(\d+)\s*matching loci/);
        $samp=$1 if(/Outputing to file \/.+\/([\w\.\d\_]+)/);
        $hap=$1 if(/(\d+) verified/);$stacks = $1 if /(\d+) stacks compared against/;print "\n$samp\t$stacks\t$num_loci\t$hap"}
        ' match_summary.log

#### Parse stacks allele file to get max,min, average number of Reads per loci
for file in *.alleles.tsv ; do printf "\n$file\t";perl -e '
  $total=$max=$numloci=0;$min=1000;
  while(<>){
  @allele=split /\s+/;
  $alcount{$allele[2]}+=$allele[5];  # for individual matches file
  }
  foreach $loci(keys %alcount){
  $numloci++;
  $total+=$alcount{$loci};
  $max=$alcount{$loci} if $max < $alcount{$loci};
  $min=$alcount{$loci} if $min > $alcount{$loci}
  }
  $average=$total/$numloci;
  print "Total:$total\tAverage:$average\tMax:$max\tMin:$min\n"
  ' < $file; done


#### Parse stacks matches file to get max,min, average number of Reads for matched loci
  for file in *.matches.tsv ; do
  printf "\n$file\t"; perl -e '
  $total=$max=$numloci=0;$min=1000;
  while(<>){
  #next if /consensus/;
  @allele=split /\s+/;
  $alcount{$allele[2].$allele[4]}+=$allele[6];
  }
  foreach $loci(keys %alcount){
  $numloci++;
  $total+=$alcount{$loci};
  $max=$alcount{$loci} if $max < $alcount{$loci};
  $min=$alcount{$loci} if $min > $alcount{$loci}
  }
  $average=$total/$numloci;
  print "Total:$total\tAverage:$average\tMax:$max\tMin:$min\n"
  ' < $file; done






for file in `ls|grep -v catalog|grep "alleles.tsv$"`;
do printf "\n\n$file\n";
perl -e '
$total=$max=$numloci=0;$min=1000;
while(<>){
@allele=split /\s+/;
$alcount{$allele[2]}+=$allele[5];
}
foreach $loci(keys %alcount){
$numloci++;
$total+=$alcount{$loci};
$max=$alcount{$loci} if $max < $alcount{$loci};
$min=$alcount{$loci} if $min > $alcount{$loci};
print "Loci:$loci\t$alcount{$loci}\n" if $alcount{$loci} >10000
}
$average=$total/$numloci;
print "\nTotal:$total\tAverage:$average\tMax:$max\tMin:$min\n"'
  < $file;
done

#### calculate ratios of raw reads per haplotype from allele file
  perl -e '
 $total=$max=$numloci=0;$min=1000;
 while(<>){
  @allele=split /\s+/;
  $alcount{$allele[2]}.="$allele[5]\t";
  $altotal{$allele[2]}+=$allele[5]
}

 foreach $loci(keys %alcount){
 next if $altotal{$loci} < 100;
 @ratio=();
  $alcount{$loci}=~s/^\s+||\s+$//g;
  my@nums=split(/\t/,$alcount{$loci});
  next if scalar @nums > 2;
  @sort_nums=sort{$a<=>$b}@nums;
  $skip=$keep=0;
  foreach my$tem(@sort_nums){push(@ratio,$tem/$sort_nums[0]); $skip=1 if $tem < 5; $keep = 1 if ($tem/$sort_nums[0] > 0)}
  next if $skip == 1;
  next if $keep == 0;
  print join("\t",@sort_nums);
  print "\tRatios:\t";
  print join("\t",@ratio);
  print "\n";
}
 ' < Qyu_RAD_index4_TGACCA_L003_R1_001_TTCTC.alleles.tsv

#### calculate ratios of raw reads per haplotype from matches file
  perl -pe '
 $total=$max=$numloci=0;$min=1000;
 my$ratio=1.5;
 my$depth=5;
 my$read=1; ## exclude the haplotype with only read. consider rest
 while(<>){
  @allele=split /\s+/;
  next if $allele[6] <= $read;
  $alcount{$allele[2]}.="$allele[6]\t";
  $altotal{$allele[2]}+=$allele[6]
}

 foreach $loci(keys %alcount){
 next if $altotal{$loci} <= $depth;
 my@ratio=();
  $alcount{$loci}=~s/^\s+//g;
  $alcount{$loci}=~s/\s+$//g;
  my@nums=split(/\t/,$alcount{$loci});
  next if scalar @nums > 2;  ## exclude loci with more than 2 haplotypes.
  next if scalar @nums == 1 ; ## exclude loci with no haplotype/SNP
  #print "\nnums >1  \t $alcount{$loci}" if scalar @nums > 2 ;
  #print "\nnums ==1 \t $alcount{$loci}" if scalar @nums ==1;

  @sort_nums=sort{$a<=>$b}@nums;
  $skip=$keep=0;
  foreach my$tem(@sort_nums){
    push(@ratio,$tem/$sort_nums[0]);
    $skip=1 if $tem < $depth;
    $keep = 1 if ($tem/$sort_nums[0] >= $ratio)
  }
  next if $skip == 1;
  next if $keep == 0;
  print join("\t",@sort_nums);
  print "\tRatios:\t";
  print join("\t",@ratio);
  print "\n";
}
 ' < Qyu_RAD_index4_TGACCA_L003_R1_001_TTCTC.matches.tsv | wc -l

  < Qyu_RAD_index4_TGACCA_L003_R1_001_AGCCC . matches . tsv | wc -l


####################################################################################################
######count the ratio of single dose to multidose markers from joinmap lgf using 1.73 as cuttof ratio

perl -e '
use List::Util qw< sum >;
use Statistics::Distributions qw< chisqrprob >;

$useratio=0;        ## use ratio to determine single dose and multidose
$usechisq=1;        ## use chisq test to determine single dose and multidose
$pat="nnxnp";
$pat="lmxll";
$plim=0.1;       ## to test deviation from 1:1 segregtaion
$indlim=20;         ## min individual should have the marker
$fplim=0.01;     ## alpha value to interpret final result
my$yates=1;
##################################
$chilim=2.77 if $plim <=0.1;
$chilim=3.84 if $plim <=0.05;
$chilim=6.63 if $plim <=0.01;
$chilim=10.83 if $plim <=0.001;
$fchilim=2.77 if $fplim <=0.1;
$fchilim=3.84 if $fplim <=0.05;
$fchilim=6.63 if $fplim <=0.01;
$fchilim=10.83 if $fplim <=0.001;
##################################
if($pat=~/nnxnp/){$one=19;$two=20;$minratio=1.73}
if($pat=~/lmxll/){$one=17;$two=18;$minratio=1.73}
if($pat=~/hkxhk/){$one=17;$two=18;$minratio=6.7}

$sd=$md=$total=$indtotal=0;

while(<>){
next if $.==1;
next if $_!~ /$pat/;
@nums=split(/\s+/);

## plain ratio
$ratio=($nums[$one]>$nums[$two]?$nums[$one]:$nums[$two])/(($nums[$one]<$nums[$two]?$nums[$one]:$nums[$two])+0.000000000000000000000000000000000000001);

# using chi square test.
$total=$nums[$one]+$nums[$two];
next if $total < $indlim;
if($usechisq==1){
$chisq=(($nums[$one] - 0.5*$total)**2)/(0.5*$total)+(($nums[$two] - 0.5*$total)**2)/(0.5*$total);
$chisq=(($nums[$one] - 0.5*$total-0.5)**2)/(0.5*$total)+(($nums[$two] - 0.5*$total-0.5)**2)/(0.5*$total) if $yates=1;;
$md++ if $chisq>$chilim;
$sd++ if $chisq<=$chilim;
}
if($useratio==1){
#print "\n",$chisq>$chilim?"MD":"SD","\t$nums[$one]\t$nums[$two]\t$chisq";
 $md++ if $ratio >= 1.73;
 $sd++ if $ratio < 1.73
}}
$indtotal=$sd+$md;
print "\n$pat\tSD\t$sd\t\tMD\t$md\t\tratio:\t",$sd/($md+$sd),":",$md/($sd+$md),"\n\n";
$auto=(($sd - 0.70*$indtotal)**2) / (0.7*$indtotal)  +  (($md - 0.30*$indtotal)**2)/(0.30*$indtotal);
$allo=(($sd - 0.56*$indtotal)**2)/ (0.56*$indtotal)  +  (($md - 0.44*$indtotal)**2)/(0.44*$indtotal);

$autopval = chisqrprob(1, $auto);
$allopval = chisqrprob(1, $allo);

print "                       \tChi-Squared\tP-Value\tSD-Obs\tSD_Exp\tMD-Obs\tMD_Exp\tResult\n";
print "Test for Autopolyploidy\t$auto\t$autopval\t$sd\t",0.70*$indtotal,"\t$md\t",0.30*$indtotal,"\t",$autopval>=$fplim?"Autopolyploid":"Not Autopolyploid"," at alpha $fplim\n";
print "Test for Allopolyploidy\t$allo\t$allopval\t$sd\t",0.56*$indtotal,"\t$md\t",0.44*$indtotal,"\t",$allopval>=$fplim?"Allopolyploid":"Not Alloopolyploid"," at alpha $fplim\n\n\n";
' < 5716_population_2.min1Ind_JoinMap.lgf


#### calculate occurence of stacks for different level of heterozygosity
 perl -e '
@header=();
my%numhet;
while(<>){
my@data=();
s/^\s+|\s+$//g;
@header=split /\s+/ if $.==1;
@data=split/\s+/ if $.>1;
for(my$i=0;$i<scalar@header;$i++){
next if $data[$i]=~/\-/;
my$count = ($data[$i] =~ s/\//\//g);
$numhet{$header[$i]}{$count+1}++}
}

foreach my$head(keys %numhet){
  foreach my$ploidy(sort{$a<=>$b} keys %{$numhet{$head}}){
  print "$head\t$ploidy\t$numhet{$head}{$ploidy}\n"
}
}

' < haplo_SNPSonly.txt


#### summarize number of sequences went through blastn using parralel script.
for file in *.fa.*.out; do printf "\n$file\tBlastn:";blastn=$(grep -c "Query=" $file);printf "$blastn\tSeq:"; seq=$(grep -c ">" ${file/.out/});printf "$seq"; if [ $seq != $blastn ]; then printf "****Blast run failed for ${file/.out/}*****"; fi; done; echo

### pick top hit based on query coverage per subject. Blast output from  -outfmt '6 std qcovhsp qcovs qlen slen'.
perl -e 'while(<>){($on,$tw,$th,$fo,$fi,$si,$se,$ei,$ni,$te,$el,$tl,$tr,$fr,$fv,$sx,@sv)=split(/\s+/,$_);if($blast{$on}{'qcov'}<$fr){$blast{$on}{'sub'}=$tw;$blast{$on}{'qcov'}=$fr;$blast{$on}{'qlen'}=$fv;$blast{$on}{'subname'}=join " ",@sv;} } foreach $query(keys %blast){print "$query\t$blast{$query}{'sub'}\t$blast{$query}{'qcov'}\t$blast{$query}{'qlen'}\t$blast{$query}{'subname'}\n"}' Ming1n2n3n4n8kb_CLCassembly_contigs.gt1k.fa.vs.Bacto.MitoCat.ChloroCat.Mono5.PineCat.alias.blastn.table > Ming1n2n3n4n8kb_CLCassembly_contigs.gt1k.fa.vs.Bacto.MitoCat.ChloroCat.Mono5.PineCat.alias.blastn.besthitOnqCovPerSubject.table


perl -e 'while(<>){my$line=$_;($on,$tw,$th,$fo,$fi,$si,$se,$ei,$ni,$te,$el,$tl,$tr,$fr,$fv,$sx,@sv)=split(/\s+/,$_);if($blast{$on}{'qcov'}<$fr){$blast{$on}{'sub'}=$tw;$blast{$on}{'qcov'}=$fr;$blast{$on}{'qlen'}=$fv;$blast{$on}{'subname'}=join " ",@sv; $blast{$on}{'all'}=$line }} foreach $query(keys %blast){print "$blast{$query}{'all'}\n"}'

## tophit based on percent id.
perl -e 'while(<>){my$line=$_;($on,$tw,$th,$fo,$fi,$si,$se,$ei,$ni,$te,$el,$tl,$tr,$fr,$fv,$sx,@sv)=split(/\s+/,$_);if($blast{$on}{'perid'}<$th){$blast{$on}{'perid'}=$th;$blast{$on}{'qlen'}=$fv;$blast{$on}{'subname'}=join " ",@sv; $blast{$on}{'all'}=$line }} foreach $query(keys %blast){print "$blast{$query}{'all'}\n"}'

perl -e 'while(<>){my$line=$_;($on,$tw,$th,$fo,$fi,$si,$se,$ei,$ni,$te,$el,$tl,$tr,$fr,$fv,$sx,@sv)=split(/\s+/,$_); if($blast{$on}{'qcov'}<$fr){$blast{$on}{'sub'}=$tw;$blast{$on}{'qcov'}=$fr;$blast{$on}{'qlen'}=$fv;$blast{$on}{'subname'}=join " ",@sv; $blast{$on}{'all'}=$line }} foreach $query(keys %blast){$len{$blast{$query}{'sub'}}{'len'}+=$blast{$query}{'qlen'};$len{$blast{$query}{'sub'}}{'subname'}=$blast{$query}{'subname'}} foreach my $sub(keys %len){print "$len{$sub}{'len'}\t$sub\t$len{$sub}{'subname'}\n"}'


## replace subject column with subject names for drawing figure with proper names instead of genbank number
 perl -e 'while(<>){my$line=$_;($on,$tw,$th,$fo,$fi,$si,$se,$ei,$ni,$te,$el,$tl,$tr,$fr,$fv,$sx,$sv)=split(/\t+/,$_);@sv=split(/\s+/,$sv); $sv[5]="" if  $sv[5]=~/\|/;$sv=join("_",$sv[0],$sv[1],$sv[2],$sv[3],$sv[4],$sv[5]);$sv=~s/[_\,\.\s\-]+$|,//g; print join ("\t",$on,$sv,$th,$fo,$fi,$si,$se,$ei,$ni,$te,$el,$tl,$tr,$fr,$fv,$sx,$tw),"\n"; }' ${seq/.table/.TopHits.BactOnly.table} > ${seq/.table/.TopHits.BactOnly.WidNames.table}

 ### create summary table for percent id similarity for prot tbale
perl -e 'while(<>){($query,$sub,$per,$len,@rest)=split /\s+/;$desc{$sub}=join "_",$rest[12],$rest[13],$rest[14];$desc{$sub}=~s/.+\[//g;$desc{$sub}=~s/\].*//g; $desc{$sub}=~s/[_\,\.\s\-]+$|,//g;$tlen{$desc{$sub}}+=$len;$qlen{$desc{$sub}}+=$rest[10];$tqlen+=$rest[10];$tperlen{$desc{$sub}}+=$per*$len;}foreach $sub(keys %tlen){$avgper=$tperlen{$sub}/$tlen{$sub}; print "\n$sub\t$avgper\t$tlen{$sub}\t$qlen{$sub}\t",$qlen{$sub}*100/$tqlen,"%\n"}'

## add length of query sequence hitting per subject species
perl -e 'while(<>){($bac,$len)=split /\s+/;$sum{$bac}+=$len} foreach my$bact(keys %sum){print "\n$bact\t$sum{$bact}"}print "\n"' Ming1n2n3n4n8kb_bacterial.TopHits.QueryLengthPerBacSpecies.txt|sort -nk2,2 > Ming1n2n3n4n8kb_bacterial.TopHits.QueryLengthPerBacSpecies.Summary.table

## add length of query sequence hitting per subject genus
perl -e 'while(<>){($bac1,$len)=split /\s+/;($bac)=split /_/,$bac1;$sum{$bac}+=$len} foreach my$bact(keys %sum){print "\n$bact\t$sum{$bact}"}print "\n"' Ming1n2n3n4n8kb_bacterial.TopHits.QueryLengthPerBacSpecies.txt|sort -nk2,2 > Ming1n2n3n4n8kb_bacterial.TopHits.QueryLengthPerBacGenus.Summary.table


#### find SSR working
perl -ne '
 $/=">";
 $motiflength=2;
 $minreps=4;
 $max_gap=20;
 $extra=40; ## extra sequence around ssr

 ## imperfect SSR for same motifs
 $regexp = "(([gatc]{2,})\\2{4,}[gatc]{0,$max_gap}\\2{4,})";
 ## perfect SSR
 $regexp = "(([gatc]{2,})\\2{4,})";
 while (<>){
 chomp;
    my ($titleline, $sequence) = split(/\n/,$_,2);
    next unless ($sequence && $titleline);
    $seqcount++;
    my ($id) = $titleline =~ /^(\S+)/; #the ID is the first whitespace-
                                       #delimited item on titleline
    $sequence =~ s/\s//g; #concatenate multi-line sequence
    study($sequence);     #is this necessary?
    while($sequence=~/$regexp/ig){

      print "\n\n$titleline\t",length($sequence),"\n",
      substr($sequence,$-[0]-$extra>=0?$-[0]-$extra:0, $+[0] + $extra <= length($sequence)?$+[0] - $-[0] + 2*$extra:length($sequence)-$-[0]+$extra),
      ":\n$1\t$2\t$-[0]\t$+[0]\n\n";
    }
 }

 ' <Raleigh_CLC_assembly.fa.masked



### find SSR in sequence and design primers..
 perl -ne '
 $/=">";

 $minreps=5;      ### minimum number of repeats
 $min_len=2;      ### minimum length of SSR motif e.g. 2 for di-nucleotide repeats
 $max_gap=20;     ### Maximum gap allowed in imperfect SSRs
 $extra=40;       ### Fetch extra sequence around ssr

 ## create regex for Imperfect Hetero SSR mining.
 #$regexp = "((([gatc]{$min_len,})\\3{$minreps,})([gatc]{0,$max_gap})(([gatc]{$min_len,})\\6{$minreps,}))"; ### $1: full match; $2:MatchBeforeGap;$3: basic motig of SSR $4: gap ;$5:match after gap; $6:second motif;

 ## perfect SSR
 $regexp = "((([gatc]{$min_len,})\\3{$minreps,}))";

 ## create regex for imperfect Homo SSR mining.
 #$regexp = "((([gatc]{$min_len,})\\3{$minreps,})([gatc]{0,$max_gap})(\\3{$minreps,}))";
 while (<>){
 chomp;
    my ($titleline, $sequence) = split(/\n/,$_,2);
    next unless ($sequence && $titleline);
    $seqcount++;
    my ($id) = $titleline =~ /^(\S+)/; #the ID is the first whitespace-
                                       #delimited item on titleline
    $sequence =~ s/\s//g; #concatenate multi-line sequence
    study($sequence);     #is this necessary?
    while($sequence=~/$regexp/ig){

      my $fullmatch=$1; #: full match;
      my $MatchBeforeGap=$2; # :MatchBeforeGap;
      my $motif=$3; ##: basic motig of SSR
      my $gapseq=$4;##: gap ;
      my $aftergap=$5;##:match after gap
      my $sec_motif=$6; ## motif after gap

      print "\n$titleline\t",
      #length($sequence),"\n",
      #substr($sequence,$-[0]-$extra>=0?$-[0]-$extra:0, $+[0] + $extra <= length($sequence)?$+[0] - $-[0] + 2*$extra:length($sequence)-$-[0]+$extra),
      ":\t$fullmatch\t\($motif\)",
      length($gapseq)>0?length($MatchBeforeGap)/length($motif):length($fullmatch)/length($motif),
      length($gapseq)>0?"\.\.$gapseq\.\.":"",
      length($gapseq)>0?join("","\($motif\)",length($aftergap)/length($motif)):"",
      "\t$-[0]\t$+[0]\n\n";

##design primers
my$primer3="/home/ratnesh.singh/softwares/primer3-2.3.6/src/primer3_core";
open PRIMER,">primer3.input";
print PRIMER "SEQUENCE_ID=$titleline
SEQUENCE_TEMPLATE=$sequence
SEQUENCE_TARGET=$-[0],",$+[0]-$-[0],"
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=22
PRIMER_MIN_SIZE=20
PRIMER_MAX_SIZE=24
PRIMER_MAX_NS_ACCEPTED=1
PRIMER_PRODUCT_SIZE_RANGE=150-250
P3_FILE_FLAG=1
SEQUENCE_INTERNAL_EXCLUDED_REGION=$-[0],",$+[0]-$-[0],"
PRIMER_EXPLAIN_FLAG=1
PRIMER_LEFT_NUM_RETURNED=1
PRIMER_RIGHT_NUM_RETURNED=1
PRIMER_PAIR_NUM_RETURNED=1
PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/home/ratnesh.singh/softwares/primer3-2.3.6/src/primer3_config/
\=
";
close PRIMER;

my $cmd="$primer3 < primer3.input";

      system($cmd);
    }
 }

 ' <Raleigh_CLC_assembly.fa.masked

 ###### test if primer sets bind on sequence.
 perl -e '
 use strict;
 my%primer;
 open PRIMER,"4Raleigh_SSR_Primer3results.txt";
      while (<PRIMER>){
      my$primer=$_;
      $primer=~s/^\s+//g;
      $primer=~s/\s+/\t/g;
      $primer=~s/\s+$//g;
      #print "Seq:$titleline\tprimer:$primer\n";
      my($num,$name,$pr1,$pr2,$ext)=split(/\s+/,$primer);
      chomp($num,$name,$pr1,$pr2,$ext);
      foreach($num,$name,$pr1,$pr2,$ext){s/\s+//gi}

     my$rev_pr2=reverse($pr2);
     $rev_pr2=~tr/ATGCatgc/TACGtacg/;

     my$rev_pr1=reverse($pr1);
     $rev_pr1=~tr/ATGCatgc/TACGtacg/;

     $primer{$num}{'name'}=$name;
     $primer{$num}{'num'}=$num;
     $primer{$num}{'pr1'}=$pr1;
     $primer{$num}{'rev_pr1'}=$rev_pr1;
     $primer{$num}{'pr2'}=$pr2;
     $primer{$num}{'rev_pr2'}=$rev_pr2;
}

  $/=">";
  open FASTA,"Raleigh_CLC_assembly.fa.masked.fa";
  my$seqcount=0;
  while (<FASTA>){
    my$seq=$_;
    my ($titleline, $sequence) = split(/\n/,$seq,2);
    next unless ($sequence && $titleline);
    $seqcount++;
    my ($id) = $titleline =~ /^(\S+)/; #the ID is the first whitespace-
                                       #delimited item on titleline
    $sequence =~ s/\s//g; #concatenate multi-line sequence
    #study($sequence);     #is this necessary?

   foreach my$pnum(sort keys %primer){

     #while($sequence=~/$primer{$pnum}{'pr1'}/ig){
     #   print "\n$id\t$primer{$pnum}{'pr1'}\t$primer{$pnum}{'num'}\t$primer{$pnum}{'name'}\t$primer{$pnum}{'pr1'}\t$primer{$pnum}{'pr2'}\t$-[0]\t$+[0]";
     #}
     # while($sequence=~/$primer{$pnum}{'pr2'}/ig){
     #   print "\n$id\t$primer{$pnum}{'pr2'}\t$primer{$pnum}{'num'}\t$primer{$pnum}{'name'}\t$primer{$pnum}{'pr1'}\t$primer{$pnum}{'pr2'}\t$-[0]\t$+[0]";
     #}
     #while($sequence=~/$primer{$pnum}{'rev_pr1'}/ig){
     #   print "\n$id\t$primer{$pnum}{'rev_pr1'}\t$primer{$pnum}{'num'}\t$primer{$pnum}{'name'}\t$primer{$pnum}{'pr1'}\t$primer{$pnum}{'pr2'}\t$-[0]\t$+[0]";
     #}
     #while($sequence=~/$primer{$pnum}{'rev_pr2'}/ig){
     #   print "\n$id\t$primer{$pnum}{'rev_pr2'}\t$primer{$pnum}{'num'}\t$primer{$pnum}{'name'}\t$primer{$pnum}{'pr1'}\t$primer{$pnum}{'pr2'}\t$-[0]\t$+[0]";
     #}
     while($sequence=~/$primer{$pnum}{'pr1'}&& $primer{$pnum}{'rev_pr2'}/ig){
        print "\nBoth\t$id\t$primer{$pnum}{'rev_pr2'}\t$primer{$pnum}{'num'}\t$primer{$pnum}{'name'}\t$primer{$pnum}{'pr1'}\t$primer{$pnum}{'pr2'}\t$-[0]\t$+[0]";
     }
     while($sequence=~/$primer{$pnum}{'pr2'}&& $primer{$pnum}{'rev_pr1'}/ig){
        print "\nBoth\t$id\t$primer{$pnum}{'rev_pr2'}\t$primer{$pnum}{'num'}\t$primer{$pnum}{'name'}\t$primer{$pnum}{'pr1'}\t$primer{$pnum}{'pr2'}\t$-[0]\t$+[0]";
     }

     }
     $/=">";
 }
 print "\n\n";
  '







  ### Count number allele at each locus
  perl -e '
  while(<>){
  next if /consensus|\-/;
  my@all=split /\s+/,$_;

  my$num=$all[3]=~s/\//\//g;

  $alfreq{'pal'}{$num+1}++;
  my$num2=$all[4]=~s/\//\//g;
  $alfreq{'5382'}{$num2+1}++;
  my$num3=$all[5]=~s/\//\//g;
  $alfreq{'0605'}{$num3+1}++;
  #print "\nAllele in Pal:$num1:$num2:$num3:$all[3]\t$all[4]\t$all[5]";
  }
  foreach my$name(keys %alfreq){
  foreach my$number(sort keys %{$alfreq{$name}}){
  print "\n$name\t$number Allele\t$alfreq{$name}{$number}"
  }
  }print "\n"' <batch_1.haplotypes_1.tsv
