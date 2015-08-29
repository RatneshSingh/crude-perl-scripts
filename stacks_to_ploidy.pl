use warnings;
use strict;
use List::Util qw< sum >;
use Statistics::Distributions qw< chisqrprob >;
use Getopt::Long;
our ( $file, $filetype, $useratio, $usechisq, $patr, $plim, $indlim, $fplim, $minratio, $yates, $freq, $as_per, $minbin, $maxbin, $bin, $details, $depth, $help,$verbose );
my $result = GetOptions(
    "joinmap|j|f=s" => \$file,
    "filetype|ft=s" => \$filetype,
    #"useratio|ur|r" => \$useratio,
    #"usechisq|uc|c" => \$usechisq,
    #"marker|m=s"    => \$patr,
    "plim|p=f"      => \$plim,
    "indlim|il=i"   => \$indlim,
    "fplim|fp=f"    => \$fplim,
    "minratio|mr=f" => \$usechisq,
    "yates|yc"      => \$yates,
    "print_frq|pf"  => \$freq,
    "aspercent|ar"  => \$as_per,
    "maxbin|mb=f"   => \$maxbin,
    "minbin|mn=f"   => \$minbin,
    "binsize|bs=f"  => \$bin,
    "details|pd"    => \$details,
    "depth|d=i"     => \$depth,
    "help|h"        => \$help,
    "verbose|v"=>\$verbose

);

my $usage = '
      "joinmap|j|f=s" => file,
      "filetype|ft=s"=>filetype,  matches, joinmap
	  "useratio|ur|r" => useratio,
      "usechisq|uc|c" => usechisq,
      "marker|m" => pat,
      "plim|p=f" => plim,
      "indlim|il=i" => indlim,
      "fplim|fp=f" => fplim,
      "minratio|mr=f" => usechisq,
      "yates|yc" => yates,
      "print_frq|pf" => freq,
      "maxbin|mb=f"   => maxbin,
      "minbin|mn=f"=>minbin,
      "binsize|bs=f"  => bin,
      "aspercent|ap"=>print Segregation as percent present in frequency table
      "details|pd"=>details,
      "help|h" => help
      "verbose|v" => verbose. print chisq comps;
      ';

die "\n\n$usage\n\n" if ( $help || !$file );
$usechisq = 1 if !$useratio;    ## use chisq test to determine single dose and multidose

#$pat=$patr?$patr:"nnxnp";
$plim     = $plim                       ? $plim     : 0.1;         ## to test deviation from 1:1 segregtaion
$indlim   = $indlim                     ? $indlim   : 10;          ## min individual should have the marker
$fplim    = $fplim                      ? $fplim    : 0.05;        ## alpha value to interpret final result
$minratio = $minratio                   ? $minratio : 1.73;
$minbin   = $minbin                     ? $minbin   : 0;
$maxbin   = $maxbin                     ? $maxbin   : 10;
$bin      = $bin                        ? $bin      : 0.25;
$filetype = ( $file =~ m/matches.tsv|alleles.tsv/i ) ? "matches" : "joinmap";
print "\nfiletype detected as :",$file=~m/matches.tsv/i?"Stacks Matches output file":$file=~m/alleles.tsv/i?"Stacks Alleles output file":"Joinmap output file","\n";
##################################

##################################
#$pat=~s/\s+//g;
#if ( $pat =~ /nnxnp/ ) { $one = 19; $two = 20; $minratio = 1.73 }
#elsif ( $pat =~ /lmxll/ ) { $one = 17; $two = 18; $minratio = 1.73 }
#elsif ( $pat =~ /hkxhk/ ) { $one = 17; $two = 18; $minratio = 6.71 }
#else{die "Cannot determine marker type *$pat*\n"}
my %sd;
my %md;
my %par;
my %segrat;
my $total = my $indtotal = 0;

##############################################################################################
### parse matches fiel to identify read ratio
##############################################################################################
if ( $filetype =~ /matches/ ) {
    my $mtotal   = 0;
    my $mmax     = 0;
    my $mnumloci = 0;
    my $min      = 1000;
    my $ratio    = 1.5;
    $depth = $depth ? $depth : 5;

    my %readperalleleeachloci;
    my %readperloci;
    my %segpat;
    my $col1=2;
    my $col2=6;
    my @readperallele;
    my $average;
    my $minreadperallele  = 1;      ## exclude the haplotype with only read. consider rest
    my $maxreadperallele  = 10000;      ## exclude the haplotype with only read. consider rest
    my $minreadperloci=1;
    my $maxreadperloci=10000;
    my %excludedloci;
    open MATCHES, "$file";

    if($file=~/alleles.tsv/){$col1=2;$col2=5;}
    if($file=~/matches.tsv/){$col1=2;$col2=6;}

    ### collect read numbers per loci into hash.
    print "\nReading $file to collect read data\n";
    while (<MATCHES>) {
        my @allele = split /\s+/;

        ## exclude the loci if any allele is less than minimum nuber or more than maximum number of reads.
        $excludedloci{$allele[$col1]}=1 if $allele[$col2] <= $minreadperallele || $allele[$col2] >$maxreadperallele;
        next if exists $excludedloci{$allele[$col1]};

        push(@{$readperalleleeachloci{ $allele[$col1] }},$allele[$col2]);
        $readperloci{ $allele[$col1] } += $allele[$col2];
        push(@readperallele,$allele[$col2]);
    }
    #my@readperloci=map{$readperloci{ $_}} keys %readperloci;

    #print "\nPrinting Frequency table for Read per allele data\n";
    #print_bin(\@readperallele,0,1000,1,"$file.readperallele.table");

    ## calculate average depth for stacks per allele
    $average=mean_selective(\@readperallele,1,30);
    my@ploidylist;
    my @ratiolist;
    my@seglist;
    my@readdistlist;
    my@alleleDosagePloidy;
    print "\nCalulting Ploidy and Read ratio per loci.";
    foreach my $loci ( keys %readperloci ) {
      next if ($readperloci{$loci}<$minreadperloci || $readperloci{$loci}>$maxreadperloci);
      next if @{$readperalleleeachloci{$loci}}<2;
        my@reads=sort{$a<=>$b} @{$readperalleleeachloci{$loci}};
        my$gcf=multigcf(@reads);

        #my($ploidy,$readdist)=test_allele_dosage(\@reads,$average,$fplim,1);

        my@segs=map{$_/$gcf} @reads;
        push(@seglist,join ":",@segs);

        #my@ratios=map{$_/min(@reads)} @reads; ## return reduced ratios to minimum
        push(@ratiolist,max(@reads)/min(@reads)) if scalar@reads ==2;

        #push(@ploidylist,$ploidy);
        #push(@readdistlist,$readdist);
        #push(@alleleDosagePloidy,"$ploidy\*$readdist");
    }

    #print "\nPrinting Frequency table for Read per Loci data\n";
    #print_bin(\@readperloci,0,1000,1,"$file.readperloci.table");
    #
    #print "\nPrinting Frequency table for Ploidy list data\n";
    #print_frequency_table(\@ploidylist,"$file.ploidylist.table");
    #
    #print "\nPrinting Frequency table for Read distribution ratios data\n";
    #print_frequency_table(\@readdistlist,"$file.readdistlist.table");
    #
    #print "\nPrinting Frequency table for Allele Dosage and Ploidy data data\n";
    #print_frequency_table(\@alleleDosagePloidy,"$file.alleleDosagePloidy.table");

    print "\nPrinting Frequency table for Read segregation data\n";
    print_frequency_table(\@seglist,"$file.readseqgregation.table");

    print "\nPrinting Bin table for Read distribution ratios data\n";
    print_bin(\@ratiolist,0,5,0.1,"$file.readratios.table");
}




sub test_read_dosage {
    my $val1  = shift;
    my $val2  = shift;
    my $alpha = shift;

    $alpha = $alpha ? $alpha : 0.05;
    my $segr = "unknown";
    if ( chisqrprob( 1, chisq2( $val1, $val2, 0.5,  0.5 ) ) > $alpha )  { $segr =~ s/unknown//g; $segr .= "1:1\t"; }
    if ( chisqrprob( 1, chisq2( $val1, $val2, 0.75, 0.25 ) ) > $alpha ) { $segr =~ s/unknown//g; $segr .= "3:1\t"; }
    $segr =~ s/^\s+|\s+$//g;
    return $segr;
}

sub test_allele_dosage {
    my $refarray  = shift;
    my $average  = shift;
    my $alpha = shift;
    use Statistics::Distributions qw< chisqrprob >;
    $alpha = $alpha ? $alpha : 0.01;
    my $seg = "unknown";

    my$numal=scalar@$refarray;
    my$sum=sum(@$refarray);
    my$depth=$sum/$average;
    my$ploidy="";

    if ($numal==1) {

        $seg =~ s/unknown//g;$seg .= "1:0,";
        $ploidy.=sprintf("%.0f",$sum/$average).",";
    }
    elsif($numal==2){

      my$chip_11=chisq( $refarray, [1,  1],'pval' ,1);
      my$chip_12=chisq( $refarray, [1,  2],'pval' ,1);
      my$chip_13=chisq( $refarray, [1,  3],'pval' ,1);
      my$chip_23=chisq( $refarray, [2,  3],'pval' ,1);
      my$chip_34=chisq( $refarray, [3,  4],'pval' ,1);
      my$chip_35=chisq( $refarray, [3,  5],'pval' ,1);
      my$chip_45=chisq( $refarray, [4,  5],'pval' ,1);
      my$chip_37=chisq( $refarray, [3,  7],'pval' ,1);



      if ($chip_11 > $alpha || $chip_12 > $alpha || $chip_13 > $alpha || $chip_23 > $alpha || $chip_34 > $alpha || $chip_35> $alpha || $chip_45> $alpha || $chip_37> $alpha  ) {
        $seg =~ s/unknown//g;
        if (max($chip_11 ,$chip_12,$chip_13,$chip_23,$chip_34 ,$chip_35,$chip_45,$chip_37)==$chip_11) {
            $seg .= "1:1,";
            if($depth < 3){$ploidy.="2,"; }
            elsif($depth > 3 && $depth < 5){ $ploidy.="4,";}
            else { $ploidy.=">4,";}
        }
        elsif (max($chip_11 ,$chip_12,$chip_13,$chip_23,$chip_34 ,$chip_35,$chip_45,$chip_37)==$chip_12) {
           $seg .= "1:2,"; $ploidy.="3,";
        }
        elsif (max($chip_11 ,$chip_12,$chip_13,$chip_23,$chip_34 ,$chip_35,$chip_45,$chip_37)==$chip_12 ) {
           $seg .= "1:3,"; $ploidy.="4,";
        }
        elsif (max($chip_11 ,$chip_12,$chip_13,$chip_23,$chip_34 ,$chip_35,$chip_45,$chip_37)==$chip_23) {
           $seg .= "2:3,"; $ploidy.="unknown,";
        }
        elsif (max($chip_11 ,$chip_12,$chip_13,$chip_23,$chip_34 ,$chip_35,$chip_45,$chip_37)==$chip_34) {
           $seg .= "3:4,"; $ploidy.="unknown,";
        }
        elsif (max($chip_11 ,$chip_12,$chip_13,$chip_23,$chip_34 ,$chip_35,$chip_45,$chip_37)==$chip_35) {
           $seg .= "3:5,"; $ploidy.="unknown,";
        }
        elsif (max($chip_11 ,$chip_12,$chip_13,$chip_23,$chip_34 ,$chip_35,$chip_45,$chip_37)==$chip_45) {
           $seg .= "4:5,"; $ploidy.="unknown,";
        }
        elsif (max($chip_11 ,$chip_12,$chip_13,$chip_23,$chip_34 ,$chip_35,$chip_45,$chip_37)==$chip_37) {
           $seg .= "3:7,"; $ploidy.="unknown,";
        }
      }
    }
    elsif($numal==3){
        my$chip_111= chisq( $refarray,  [1, 1, 1] ,'pval' ,1) ;
        my$chip_112=chisq( $refarray,  [1, 1, 2]  ,'pval' ,1) ;
        if ($chip_111 > $alpha || $chip_112 > $alpha) {

           if ($chip_111 > $chip_112 ) {
             $seg =~ s/unknown//g; $seg .= "1:1:1,"; $ploidy.="3,";
           }
           elsif($chip_112 > $chip_112 ){
             $seg =~ s/unknown//g; $seg .= "1:1:2,"; $ploidy.="4,";
           }
        }
    }
    elsif($numal==4){
        if (chisq( $refarray,  [1, 1, 1, 1] ,'pval' ,1) > $alpha  )  { $seg =~ s/unknown//g; $seg .= "1:1:1:1,"; $ploidy.="4,"; }
    }
    else{
        $ploidy="unknown,";
    }
    $seg =~ s/^\s+|\s+$//g;
    print "\n@$refarray\t\tRatios:$seg\tploidy:$ploidy" if $verbose;
    return ($ploidy,$seg);
}
#
##############################################################################################
# print histogram data to create grapgh.
##############################################################################################
sub print_bin {
    my $red_array = shift;
    my $min       = shift;
    my $max       = shift;
    my $sbin      = shift;
#    my $slide=shift;
    my $outfile   = shift;

    $min  = $min  ? $min  : '0';
    $max  = $max  ? $max  : 100;
    $sbin = $sbin ? $sbin : 10;
    $outfile=$outfile?$outfile:"bin.table";
    my %binhash;
    my $btotal;
    open (BIN,">$outfile") or die "Could not open $outfile to save bin result";
    print BIN "\ntotal number of loci:",scalar @$red_array,"\n\nrange\tFreq\t\%";

    ## initialize all the values to zero
    for ( my $i = $min; $i <= $max; $i += $sbin ) {$binhash{$i}=0;}

    foreach my $val ( sort { $a <=> $b } @$red_array ) {
        next if $val=~/^\s*$/;
        $btotal++;
        for ( my $i = $min; $i <= $max; $i += $sbin ) {
            #$binhash{$i}++   if ( $val > $i  && $val <= $i+$sbin );
            #$binhash{$i+$sbin}++ if ( $i+$sbin >= $max && $val > $i+$sbin );

            $binhash{$i}++   if ( $val > $i  && $val <= ($i+$sbin<=$max?$i+$sbin:$max) );
            $binhash{$max}++ if ( $i+$sbin > $max && $val > $max );

        }
    }
    foreach my $range ( sort { $a <=> $b } keys %binhash ) {
      #print "\n range:$range";
        print BIN "\n",
          $range < $max ? round( $range, 2 ) . "-" . round( ( ($range+$sbin<=$max?$range+$sbin:$max) ), 2 ) : ">$max",
          "\t$binhash{$range}\t", round( $binhash{$range} * 100 / $btotal, 2 );
          #print "\n",
          #($range + $sbin) <= $max ? $range. "-" . ($range + $sbin) : ">$max",
          #"\t$binhash{$range}\t", round( $binhash{$range} * 100 / $btotal, 2 );
    }
    print "\n\n";
}
#
##############################################################################################
# print frequency table to create histogram
##############################################################################################
sub print_frequency_table {
    my $red_array = shift;
    my $outfile   = shift;
    $outfile=$outfile?$outfile:"Frequency.table";
    my %binhash;
    my $btotal;
    open (BIN,">$outfile") or die "Could not open $outfile to save bin result";
    print BIN "\ntotal number of loci:",scalar @$red_array,"\n\nrange\tFreq\t\%";
    foreach my $val ( @$red_array ) {
        next if $val=~/^\s*$/;
            $btotal++;
            $binhash{$val}++;
        }
    print BIN "\n \tFrequency\t\%Frequency";
    foreach my $range ( sort keys %binhash ) {
        print BIN "\n$range\t$binhash{$range}\t",round( $binhash{$range} * 100 / $btotal, 2 );
    }
    print "\n\n";
}




##############################################################################################
# return minimum of all the numbers
##############################################################################################
sub min {
    @_ = sort { $a <=> $b } @_;
    return $_[0];

}
##############################################################################################
# return maximum of all the numbers
##############################################################################################
sub max {

    @_ = sort { $a <=> $b } @_;
    return $_[-1];

}
##############################################################################################
# return rounded number to decimal digits
##############################################################################################
sub round {
    my $number   = shift;
    my $decimals = shift;

    #return "0" if $number == 0;
    $decimals = $decimals ? $decimals : 3;
    return substr( $number + ( '0.' . '0' x $decimals . '5' ), 0, $decimals + length( int($number) ) + 1 );
}

##############################################################################################
# return chi-square value for observed values  and expected (values OR decimal ratios (sum 1.0) OR full digit ratios (eg 1:4))
# expected values will converted to decimal ratios and reconverted to expected values in calculations.
# any type of expected values can be used (e.g. 0.25:0.75 OR 1:3 OR 10:30 ).
# returns chisq value, degree of freedom and pvalue.
##############################################################################################
sub chisq {
    use Statistics::Distributions qw< chisqrprob >;
    use List::Util qw< sum >;
    my $ref_array_obs = shift;
    my $ref_array_rat = shift;
    my $return= shift;
    my $yatecor  = shift;


    my $chisq=0;
    my $pval=1;
    my $total=sum(@$ref_array_obs);

    # create ratios from expected values to make it usefull with exp values, ratios with decimal or number
    my $totalexp=sum(@$ref_array_rat);
    my @ratios=map{sprintf("%.5f",$_/$totalexp)} @$ref_array_rat;
    my $degree= scalar @$ref_array_obs -1;

    for(my$i=0;$i<scalar@$ref_array_obs;$i++){
       $chisq += ( ( $$ref_array_obs[$i] - $ratios[$i]*$total )**2 ) / ( $ratios[$i]+0.00000000000000000000001) if !$yatecor;
       $chisq += ( abs( $$ref_array_obs[$i] - $ratios[$i]*$total - 0.5)**2 ) / ( $ratios[$i]+0.00000000000000000000001) if $yatecor;
    }
    print "\n\tObs:@$ref_array_obs\tExp:@$ref_array_obs\tRatios:@$ref_array_rat\tchisq:$chisq\tprob:",sprintf("%.6f",chisqrprob($degree,$chisq)),"\tdegree:$degree"  if $verbose;
    if ($return=~/chisq/i) {return $chisq}
    elsif($return=~/pval/i){return sprintf("%.6f",chisqrprob($degree,$chisq))}
    elsif($return=~/deg/i){return ($degree)}
    elsif($return==1){return ($chisq)}
    elsif($return==2){return ($chisq,sprintf("%.6f",chisqrprob($degree,$chisq)))}
    elsif($return==3){return ($chisq,sprintf("%.6f",chisqrprob($degree,$chisq,$degree)))}
    else{return ($chisq)}
}


##############################################################################################
# return mean of all the numbers between max and minimum
##############################################################################################
sub mean_selective {
	my$ref_array=shift;
	my$min=shift;
	my$max=shift;
	$min=$min?$min:0;
	$max=$max?$max:1e100;
	my$sum=my$num=0;
	foreach my$val(@$ref_array){
		next if ($val < $min || $val >$max );
		$sum+=$val;
		$num++;
	}
	return sprintf("%.2f",$sum/$num);
}

##############################################################################################
# return reduced fraction.
# copied from "http://www.perlmonks.org/?node_id=99531"
##############################################################################################
sub reduce {
  "@{[map 1x$_,@_]}"=~/(1+)\1* \1+$/;map$_/$+[1],@_
}
##############################################################################################
# return gcf greatest common factor.
# copied from "http://www.perlmonks.org/?node_id=99531"
##############################################################################################

sub multigcf {
  my $x = shift;
  $x = gcf($x, shift) while @_;
  return $x;
}
##############################################################################################
# return least common multiple
# copied from "http://www.perlmonks.org/?node_id=99531"
##############################################################################################
sub multilcm {
  my $x = shift;
  $x = lcm($x, shift) while @_;
  return $x;
}
sub gcf {
  my ($x, $y) = @_;
  ($x, $y) = ($y, $x % $y) while $y;
  return $x;
}

sub lcm {
  return($_[0] * $_[1] / gcf($_[0], $_[1]));
}
