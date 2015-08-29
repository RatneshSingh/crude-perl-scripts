use warnings;
use strict;
use List::Util qw< sum >;
use Statistics::Distributions qw< chisqrprob >;
use Getopt::Long;
our ( $file, $filetype, $useratio, $usechisq, $patr, $plim, $indlim, $fplim, $minratio, $yates, $freq, $as_per, $minbin, $maxbin, $bin, $details, $depth, $help );
my $result = GetOptions(
    "joinmap|j|f=s" => \$file,
    "filetype|ft=s" => \$filetype,
    "useratio|ur|r" => \$useratio,
    "usechisq|uc|c" => \$usechisq,
    "marker|m=s"    => \$patr,
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
    "help|h"        => \$help

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
$filetype = ( $file =~ m/matches.tsv|alleles.tsv/ ) ? "matches" : "joinmap";
print "\nfiletype detected as :",$file=~m/matches.tsv/?"Stacks Matches output file":$file=~m/alleles.tsv/?"Stacks Alleles output file":"Joinmap output file","\n";
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
### parse joinmap lgf file to identify ratio of progeny
##############################################################################################

if ( $filetype =~ /joinmap/ ) {
    open FILE, "$file" or die "\n\nCannot open $file \n\n\n";
    while (<FILE>) {
        next if $. == 1;
        my $one = my $two = my $three = 5;
        my @nums = split(/\s+/);
        my $pat = $patr ? $patr : $nums[3];
        $pat =~ s/\<|\>|\s+//g;
        ## set the column number to use and max ratio ratio to define SD depending on marker type
        if    ( $pat =~ /nnxnp/ ) { $one = 19; $two = 20; $minratio = 1.73; $par{$pat}{'col'} = "$one,$two"; $par{$pat}{'rat'} = $minratio; }
        elsif ( $pat =~ /lmxll/ ) { $one = 17; $two = 18; $minratio = 1.73; $par{$pat}{'col'} = "$one,$two"; $par{$pat}{'rat'} = $minratio; }
        elsif ( $pat =~ /hkxhk/ ) { $one = 12; $two = 13; $three = 14; $minratio = 6.71; $par{$pat}{'col'} = "$one,$two,$three"; $par{$pat}{'rat'} = $minratio; }
        else                      { next; }
        my $chisq;
        my $segratio;

        # using chi square test.
        $total = $nums[$one] + $nums[$two] + $nums[$three];

        next if $total < $indlim;
        if ($usechisq) {

            if ( $pat =~ /hkxhk/ ) {
                $chisq = chisq3( $nums[$one], $nums[$two], $nums[$three], 0.25, 0.50, 0.25, $yates );
                $segratio = max( $nums[$one] + $nums[$two], $nums[$two] + $nums[$three] ) / (min( $nums[$one], $nums[$three] )+0.00000000000000000000001);
            }
            else {

                $chisq    = chisq2( $nums[$one], $nums[$two], 0.5, 0.5, $yates );
                $segratio = $nums[$two] / ($nums[$one]+0.00000000000000000000001);
                $segratio = $nums[$one] * 100 / $total if $as_per;
            }
            $md{$pat}++ if chisqrprob( 1, $chisq ) < $plim;
            $sd{$pat}++ if chisqrprob( 1, $chisq ) >= $plim;

        }

        #next if ($nums[$one]==0 || $nums[$two]==0);
        elsif ($useratio) {

            #print "\n",$chisq>$chilim?"MD":"SD","\t$nums[$one]\t$nums[$two]\t$chisq";
            $nums[$one]   += 0.00000000000000000000001 if ( $nums[$one] == 0 );
            $nums[$two]   += 0.00000000000000000000001 if ( $nums[$two] == 0 );
            $nums[$three] += 0.00000000000000000000001 if ( $nums[$three] == 0 );

            my $ratio;

            if ( $pat =~ /hkxhk/ ) {
                $ratio = max( $nums[$one] + $nums[$two], $nums[$two] + $nums[$three] ) / min( $nums[$one] + $nums[$three] );

                #$segratio=($nums[$one]+$nums[$two])*100/$total;
                $segratio = $ratio;
                $segratio = max( $nums[$one] + $nums[$two], $nums[$two] + $nums[$three] ) * 100 / $total if $as_per;
            }
            else {
                $ratio = $nums[$two] / $nums[$one];    # $nums[$one] : $nums[$two] ) + 0.000000000000000000000000000000000000001 );
                $segratio = $ratio;    ## ratio presence vs absense. lmxll and nnxnp in both cases column [one] in homozyg(absense) and column [two] in hetero(presence)
                $segratio = $nums[$one] * 100 / $total if $as_per;
            }

            for ( my $i = $minbin; $i <= $maxbin; $i += $bin ) {
                $segrat{$pat}{$i}++ if ( $segratio > $i && $segratio <= $i + $bin );
                $segrat{$pat}{$i}++ if ( $i >= $maxbin  && $segratio > $maxbin );
            }

            $md{$pat}++ if $ratio >= $minratio;
            $sd{$pat}++ if $ratio < $minratio;

        }
    }
## print results
foreach my $mtype ( keys %sd ) {

my $mtotal = $sd{$mtype} + $md{$mtype};
print "\n========================================================\nAnalysis marker type $mtype\nTotal $mtotal";
print $useratio? "\nCuttof Ratio $par{$mtype}{'rat'}\n" : "\nUsed Chi-sq with plimit $plim\n";
print "Parameters:
Marker type:$mtype
Columns used:$par{$mtype}{'col'}";
print $useratio? "\nCuttof Ratio $par{$mtype}{'rat'}\n" : "\nUsed Chi-sq with plimit $plim\npvalue limit for chi-sq determination of single dose:$plim";
print "\nYates correction:", $yates ? "ON" : "OFF", "min individual should have the marker:$indlim alpha value to interpret final result:$fplim\n";

        print "\n$mtype\tSD\t$sd{$mtype}\t\tMD\t$md{$mtype}\t\tratio:\t", $sd{$mtype} / $mtotal, ":", $md{$mtype} / $mtotal, "\n\n";

        my $auto = chisq2( $sd{$mtype}, $md{$mtype}, 0.7,  0.3,  $yates );
        my $allo = chisq2( $sd{$mtype}, $md{$mtype}, 0.56, 0.44, $yates );
        my $autopval = chisqrprob( 1, $auto );
        my $allopval = chisqrprob( 1, $allo );

        print "                       \tChi-Squared\tP-Value\tSD-Obs\tSD_Exp\tMD-Obs\tMD_Exp\tResult\n";
        print "Test for Autopolyploidy\t$auto\t$autopval\t$sd{$mtype}\t", 0.70 * $mtotal, "\t$md{$mtype}\t", 0.30 * $mtotal, "\t",
          $autopval >= $fplim ? "Autopolyploid" : "Not Autopolyploid",
          " at alpha $fplim\n";
        print "Test for Allopolyploidy\t$allo\t$allopval\t$sd{$mtype}\t", 0.56 * $mtotal, "\t$md{$mtype}\t", 0.44 * $mtotal, "\t",
          $allopval >= $fplim ? "Allopolyploid" : "Not Alloopolyploid",
          " at alpha $fplim\n\n\n";

        print "\n\nSegregationRatio\tFrequency\tPercent";
        foreach my$bins ( sort { $a <=> $b } keys %{ $segrat{$mtype} } ) {
            print "$bins:$bin\n";
            print "\n",
            $bins + $bin < $maxbin ? round( $bins, 2 ) . "-" . round( ( $bins + $bin ), 2 ) : "\>$maxbin",
            "\t$segrat{$mtype}{$bins}\t",
            round( $segrat{$mtype}{$bins} * 100 / $mtotal, 2 );
        }

        print "\n\n";


    }

}
##############################################################################################
### parse matches fiel to identify read ratio
##############################################################################################
elsif ( $filetype =~ /matches/ ) {
    my $mtotal   = 0;
    my $mmax     = 0;
    my $mnumloci = 0;
    my $min      = 1000;
    my $ratio    = 1.5;
    $depth = $depth ? $depth : 5;
    my $read  = 1;      ## exclude the haplotype with only read. consider rest
    my %alcount;
    my %altotal;
    my %segpat;
    my $col1=2;
    my $col2=6;
    my @depth;
    my $average=1;
    open MATCHES, "$file";

    if($file=~/alleles.tsv/){$col1=2;$col2=5;}
    if($file=~/matches.tsv/){$col1=2;$col2=6;}


    while (<MATCHES>) {
        my @allele = split /\s+/;
        next if $allele[$col2] <= $read;
        $alcount{ $allele[$col1] } .= "$allele[$col2]\t";
        $altotal{ $allele[$col1] } += $allele[$col2];
        push(@depth,$allele[$col2]);
    }

    print_bin(\@depth,0,1000,1);

    ## calculate average depth for stacks
    foreach my $loci ( keys %alcount ) {
        $average=mean_selective(\@depth,1,200);
    }



    my @all_ratios;

    foreach my $loci ( keys %alcount ) {
        next if $altotal{$loci} <= $depth; ### exclude loci if sequencing depth is lower than threshold.
        my @ratio = ();
        $alcount{$loci} =~ s/^\s+|\s+$//g;
        #$alcount{$loci} =~ s/\s+$//g;
        my @nums = split( /\t/, $alcount{$loci} );
        #next if scalar @nums > 2;     ## exclude loci with more than 2 haplotypes.
        #next if scalar @nums == 1;    ## exclude loci with no haplotype/SNP
                                      #print "\nnums >1  \t $alcount{$loci}" if scalar @nums > 2 ;
                                      #print "\nnums ==1 \t $alcount{$loci}" if scalar @nums ==1;
        ##test for segregation pattern
        my $seg = test_read_dosage( $nums[0], $nums[1], $plim );
        my($ploidy,$seqregation)=test_allele_dosage(\@nums,$average,0.05);
        $segpat{$seg}++;
        print "\n$nums[0]:$nums[1]\t$seg" if $details;

        my @sort_nums = sort { $a <=> $b } @nums;
        push( @all_ratios, $sort_nums[1] / $sort_nums[0] );
        my $skip = my $keep = 0;




        ## decide if to keep or skip the data
        foreach my $tem (@sort_nums) {
            push( @ratio, $tem / $sort_nums[0] );

            $skip = 1 if $tem < $depth;
            $keep = 1 if ( $tem / $sort_nums[0] >= $ratio );
        }
        next if $skip == 1;
        next if $keep == 0;
        if ($details) {
            print "\n";
            print join( "\t", @sort_nums );
            print "\tRatios:\t";
            print join( "\t", @ratio );
            print "\n";
        }
    }

    ###print the frequency of read ratios
    print "\n\nMajorSNP:MinorSNP\tFrequency\t\%_of_loci";
    foreach my $segtype ( keys %segpat ) {my$st=$segtype;$st=~s/\s+/\,/g; print "\n$st\t$segpat{$segtype}\t",round($segpat{$segtype}*100/@all_ratios,2) }
    print "\n\n";

    print_bin( \@all_ratios, $minbin, $maxbin, $bin );

}





##############################################################################################
# return chi-square value for 2 numbers and expected ratio
##############################################################################################
sub chisq2 {
    my $val1     = shift;
    my $val2     = shift;
    my $rat1     = shift;
    my $rat2     = shift;
    my $yatecor  = shift;
    my $indtotal = $val1 + $val2;

    my $chisq2 = ( ( $val1 - $rat1 * $indtotal )**2 ) / ( $rat1 * $indtotal +0.00000000000000000000001) + ( ( $val2 - $rat2 * $indtotal )**2 ) / ( $rat2 * $indtotal+0.00000000000000000000001 );

    $chisq2 = ( ( abs( $val1 - $rat1 * $indtotal ) - 0.5 )**2 ) / ( $rat1 * $indtotal +0.00000000000000000000001) + ( ( abs( $val2 - $rat2 * $indtotal ) - 0.5 )**2 ) / ( $rat2 * $indtotal+0.00000000000000000000001 ) if $yatecor;

    return ($chisq2);

}

sub chisq3 {
    my $val1     = shift;
    my $val2     = shift;
    my $val3     = shift;
    my $rat1     = shift;
    my $rat2     = shift;
    my $rat3     = shift;
    my $yatecor  = shift;
    my $indtotal = $val1 + $val2 + $val3;

    my $chisq3 =
      ( ( $val1 - $rat1 * $indtotal )**2 ) / ( $rat1 * $indtotal ) +
      ( ( $val2 - $rat2 * $indtotal )**2 ) / ( $rat2 * $indtotal ) +
      ( ( $val3 - $rat3 * $indtotal )**2 ) / ( $rat3 * $indtotal );

    $chisq3 =
      ( ( abs( $val1 - $rat1 * $indtotal ) - 0.5 )**2 ) / ( $rat1 * $indtotal ) +
      ( ( abs( $val2 - $rat2 * $indtotal ) - 0.5 )**2 ) / ( $rat2 * $indtotal ) +
      ( ( abs( $val3 - $rat3 * $indtotal ) - 0.5 )**2 ) / ( $rat3 * $indtotal )
      if $yatecor;

    return ($chisq3);

}

#####
sub test_seg {
    my $val1  = shift;
    my $val2  = shift;
    my $alpha = shift;

    $alpha = $alpha ? $alpha : 0.05;
    my $seg = "unknown";
    if ( chisqrprob( 1, chisq2( $val1, $val2, 0.5,  0.5 ) ) > $alpha )  { $seg =~ s/unknown//g; $seg .= "1:1\t"; }
    if ( chisqrprob( 1, chisq2( $val1, $val2, 0.75, 0.25 ) ) > $alpha ) { $seg =~ s/unknown//g; $seg .= "3:1\t" }
    if ( chisqrprob( 1, chisq2( $val1, $val2, 0.83, 0.17 ) ) > $alpha ) { $seg =~ s/unknown//g; $seg .= "5:1\t" }
    if ( chisqrprob( 1, chisq2( $val1, $val2, 0.92, 0.08 ) ) > $alpha ) { $seg =~ s/unknown//g; $seg .= "11:1\t" }
    if ( chisqrprob( 1, chisq2( $val1, $val2, 0.97, 0.03 ) ) > $alpha ) { $seg =~ s/unknown//g; $seg .= "35:1\t" }
    $seg =~ s/^\s+|\s+$//g;
    return $seg;

}

sub test_read_dosage {
    my $val1  = shift;
    my $val2  = shift;
    my $alpha = shift;

    $alpha = $alpha ? $alpha : 0.05;
    my $segr = "unknown";
    if ( chisqrprob( 1, chisq2( $val1, $val2, 0.5,  0.5 ) ) > $alpha )  { $segr =~ s/unknown//g; $segr .= "1:1\t"; }
    if ( chisqrprob( 1, chisq2( $val1, $val2, 0.75, 0.25 ) ) > $alpha ) { $segr =~ s/unknown//g; $segr .= "3:1\t" }
    #if ( chisqrprob( 1, chisq2( $val1, $val2, 0.83, 0.17 ) ) > $alpha ) { $seg =~ s/unknown//g; $seg .= "5:1\t" }
    #if ( chisqrprob( 1, chisq2( $val1, $val2, 0.92, 0.08 ) ) > $alpha ) { $seg =~ s/unknown//g; $seg .= "11:1\t" }
    #if ( chisqrprob( 1, chisq2( $val1, $val2, 0.97, 0.03 ) ) > $alpha ) { $seg =~ s/unknown//g; $seg .= "35:1\t" }

    $segr =~ s/^\s+|\s+$//g;
    return $segr;

}

sub test_allele_dosage {
    my $refarray  = shift;
    my $average  = shift;
    my $alpha = shift;
    use Statistics::Distributions qw< chisqrprob >;
    $alpha = $alpha ? $alpha : 0.05;
    my $seg = "unknown";

    my$numal=scalar@$refarray;
    my$sum=sum(@$refarray);
    my$ploidy=0;
    if ($numal==1) {
        $ploidy=sprintf("%.0f",$$refarray[0]/$average);
    }
    elsif($numal==2){

        my@dosage=map{$_/$average} sort{$a<=>$b}@$refarray;
        if ( chisqrprob( 1, (chisq( $refarray, [1,  1] ) )) > $alpha )  { $seg =~ s/unknown//g; $seg .= "1:1\t"; $ploidy=2; }
        if ( chisqrprob( 1, (chisq( $refarray, [1,  2] ) )) > $alpha )  { $seg =~ s/unknown//g; $seg .= "1:2\t"; $ploidy=3;}
        if ( chisqrprob( 1, (chisq( $refarray, [2,  2] ) )) > $alpha && $sum/$average >=2.7 )  { $seg =~ s/unknown//g; $seg .= "2:2\t"; $ploidy=4;}
        if ( chisqrprob( 1, (chisq( $refarray, [1,  3] ) )) > $alpha )  { $seg =~ s/unknown//g; $seg .= "1:3\t"; $ploidy=4;}
    }
    elsif($numal==3){
        my@dosage=map{$_/$average} sort{$a<=>$b}@$refarray;
        if ( chisqrprob( 2, (chisq( $refarray, [1, 1, 1] ) )) > $alpha )  { $seg =~ s/unknown//g; $seg .= "1:1:1\t"; $ploidy=3;}
        if ( chisqrprob( 2, (chisq( $refarray, [1, 1, 2] ) )) > $alpha && $sum/$average >=2.7 )  { $seg =~ s/unknown//g; $seg .= "1:1:2\t"; $ploidy=4;}
    }
    elsif($numal==4){
        if ( chisqrprob( 2, (chisq( $refarray, [1, 1, 1, 1] ) )) > $alpha )  { $seg =~ s/unknown//g; $seg .= "1:1:1:1\t"; $ploidy=4;}
    }
    else{
        $ploidy=">4";
    }
    $seg =~ s/^\s+|\s+$//g;

    return ($ploidy,$seg);

}



sub print_bin {
    my $red_array = shift;
    my $min       = shift;
    my $max       = shift;
    my $sbin      = shift;
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
# return chi-square value for 2 numbers and expected ratio
##############################################################################################
sub chisq2ratio {
    my $val1     = shift;
    my $val2     = shift;
    my $rat1     = shift;
    my $rat2     = shift;
    my $yatecor  = shift;
    my $indtotal = $val1 + $val2;

    my $chisq2 = ( ( $val1 - $rat1 * $indtotal )**2 ) / ( $rat1 * $indtotal +0.00000000000000000000001) + ( ( $val2 - $rat2 * $indtotal )**2 ) / ( $rat2 * $indtotal+0.00000000000000000000001 );

    $chisq2 = ( ( abs( $val1 - $rat1 * $indtotal ) - 0.5 )**2 ) / ( $rat1 * $indtotal +0.00000000000000000000001) + ( ( abs( $val2 - $rat2 * $indtotal ) - 0.5 )**2 ) / ( $rat2 * $indtotal+0.00000000000000000000001 ) if $yatecor;

    return ($chisq2);

}
##############################################################################################
# return chi-square value for 3 numbers and expected ratio
##############################################################################################
sub chisq3ratio {
    my $val1     = shift;
    my $val2     = shift;
    my $val3     = shift;
    my $rat1     = shift;
    my $rat2     = shift;
    my $rat3     = shift;
    my $yatecor  = shift;
    my $indtotal = $val1 + $val2 + $val3;

    my $chisq3 =
      ( ( $val1 - $rat1 * $indtotal )**2 ) / ( $rat1 * $indtotal ) +
      ( ( $val2 - $rat2 * $indtotal )**2 ) / ( $rat2 * $indtotal ) +
      ( ( $val3 - $rat3 * $indtotal )**2 ) / ( $rat3 * $indtotal );

    $chisq3 =
      ( ( abs( $val1 - $rat1 * $indtotal ) - 0.5 )**2 ) / ( $rat1 * $indtotal ) +
      ( ( abs( $val2 - $rat2 * $indtotal ) - 0.5 )**2 ) / ( $rat2 * $indtotal ) +
      ( ( abs( $val3 - $rat3 * $indtotal ) - 0.5 )**2 ) / ( $rat3 * $indtotal )
      if $yatecor;

    return ($chisq3);

}
###############################################################################################
## return chi-square value for 2 observed and expected
###############################################################################################
#sub chisq2exp {
#    my $ref_array_obs = shift;
#    my $ref_array_exp = shift;
#    my $yatecor  = shift;
#    my $chisq2=0;
#
#
#
#    for(my$i=0;$i<scalar@$ref_array_obs;$i++){
#       $chisq2 += ( ( $$ref_array_obs[$i] - $$ref_array_exp[$i] )**2 ) / ( $$ref_array_exp[$i]+0.00000000000000000000001) if !$yatecor;
#       $chisq2 += ( ( $$ref_array_obs[$i] - $$ref_array_exp[$i] - 0.5)**2 ) / ( $$ref_array_exp[$i]+0.00000000000000000000001) if $yatecor;
#
#    }
#    return ($chisq2);
#}

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
    my $yatecor  = shift;


    my $chisq=0;
    my $pval=1;
    my $total=sum(@$ref_array_obs);

    # create ratios from expected values to make it usefull with exp values, ratios with decimal or number
    my $totalexp=sum(@$ref_array_rat);
    my @ratios=map{sprintf("%.2f",$_/$totalexp)} @$ref_array_rat;
    my $degree= scalar @$ref_array_obs -1;

    for(my$i=0;$i<scalar@$ref_array_obs;$i++){
       $chisq += ( ( $$ref_array_obs[$i] - $ratios[$i]*$total )**2 ) / ( $ratios[$i]+0.00000000000000000000001) if !$yatecor;
       $chisq += ( abs( $$ref_array_obs[$i] - $ratios[$i]*$total - 0.5)**2 ) / ( $ratios[$i]+0.00000000000000000000001) if $yatecor;
    }
    return ($chisq,$degree,chisqrprob($degree,$chisq));
}
##############################################################################################
# return chi-square obsue for 3 numbers and expected expio
##############################################################################################
sub chisq3exp {
    my $obs1     = shift;
    my $obs2     = shift;
    my $obs3     = shift;
    my $exp1     = shift;
    my $exp2     = shift;
    my $exp3     = shift;
    my $yatecor  = shift;
    my $indtotal = $obs1 + $obs2 + $obs3;

    my $chisq3 =
      ( ( $obs1 - $exp1 * $indtotal )**2 ) / ( $exp1 * $indtotal ) +
      ( ( $obs2 - $exp2 * $indtotal )**2 ) / ( $exp2 * $indtotal ) +
      ( ( $obs3 - $exp3 * $indtotal )**2 ) / ( $exp3 * $indtotal );

    $chisq3 =
      ( ( abs( $obs1 - $exp1 * $indtotal ) - 0.5 )**2 ) / ( $exp1 * $indtotal ) +
      ( ( abs( $obs2 - $exp2 * $indtotal ) - 0.5 )**2 ) / ( $exp2 * $indtotal ) +
      ( ( abs( $obs3 - $exp3 * $indtotal ) - 0.5 )**2 ) / ( $exp3 * $indtotal )
      if $yatecor;

    return ($chisq3);

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
