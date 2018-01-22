#!/usr/bin/perl
use warnings;
use strict;

our ( @header, $filter_for_ratio, $skip, @all_alleles, @nnnp_alleles, @lmll_alleles );

my $usage = "
usage:
$0 VCF_file_name



";

print "$usage" if !$ARGV[0];

open( VCF, $ARGV[0] ) or die "Unable to find vcf file $ARGV[0]";
$filter_for_ratio = 1;
my $mname    = "LA_Purple";
my $fname    = "US56_14_4";
my $f1_name  = "9202";
my $pop_type = "F2";
my $alpha    = 0.0000000000000000000000000004;

open(OUT,">$ARGV[0].allele_ratios_for_plotiing.table");

### read VCF file file to collect information
while ( my $line = <VCF> ) {

    ### collect header information
    if ( $line =~ /^#CHROM/ ) { @header = split( /\s+/, $line ); s/\s+//g foreach @header; print "\nLine Header:", join ":", @header; next }

    ### collect marker information
    my @markers = split( /\s+/, $line );
    s/\s+//g foreach @markers;

    ### first 4 colums are useless for this purpose so process after that
    my%hmarkers;
    for ( my $i = 4; $i < scalar @markers; $i++ ) {
        $hmarkers{ $header[$i] } = $markers[$i];
    }

    ###
    #next if ( $hmarkers{$mname} eq $hmarkers{$fname} );    ## skip both parents have same marker
    #next if ( $hmarkers{$mname} eq "-" || $hmarkers{$fname} eq "_" || $hmarkers{$f1_name} eq "_");    ## skip if one of the parentof F1 progeny has missing marker
    #next if ( sort_order($hmarkers{$f1_name}) ne sort_order(non_redundant_allele($hmarkers{$mname} , $hmarkers{$fname})));
    my @name    = @header[ 4 .. $#markers ];
    my @alleles = @markers[ 4 .. $#markers ];

    my $jmline = allele_to_joinmap( \@name, \@alleles, $mname, $fname );
    my($segval,$seg_al,$seg_rat,$alrat)=estimate_segregation(\%hmarkers, \@name, \@alleles, $mname, $fname, $f1_name );
    #print "\n$mname:$hmarkers{$mname}\t$fname:$hmarkers{$fname}\t$f1_name:$hmarkers{$f1_name}\t$segval\t\t$seg_al\t\t$seg_rat\t\t$alrat";
    print OUT "\n$alrat" if $alrat;
    $alrat=~s/^\s+//g;
    
    my($homo_rat,$het_rat)=split /\s+/,$alrat;
    my($chisq, $degree, $chisqrprob)=chisq([$homo_rat,$het_rat],[1,3],1) if ($homo_rat && $het_rat);
    
     print "\n$mname:$hmarkers{$mname}\t$fname:$hmarkers{$fname}\t$f1_name:$hmarkers{$f1_name}\t$seg_al\t$seg_rat\t$chisqrprob",join("\t",) if ($homo_rat && $het_rat && $chisqrprob > 0.00000001);
    
    
    
    
    #next if $jmline =~ /xx/i;
    #my $type = $1 if $jmline =~ /\<([\S]+)\>/;
    ##### test the required ratio (1:1 in this case) and skip
    #if ($filter_for_ratio) {
    #    $skip = undef;
    #    my @geno = split( /\s+/, $jmline );
    #    shift(@geno);
    #    my %genot;
    #    foreach (@geno) { $genot{$_}++ }
    #    my @genoty;
    #    push( @genoty, values %genot );
    #    next if scalar @genoty < 2;
    #    my @ratio = ( 1, 1 );
    #    my ( $chisq, $degree, $pval ) = chisq( \@genoty, \@ratio );
    #    print "\n@genoty,$chisq,$degree,$pval\n";
    #    $skip = 1 if ( $pval < $alpha );
    #}
    #next if $skip;
    #
    #print "\n$type\t@alleles[-3..-1]";
    #
    #push( @all_alleles, join( "\t", "$markers[0]\_$markers[1]", $jmline ) );
    #
    #push( @nnnp_alleles, join( "\t", "$markers[0]\_$markers[1]", $jmline ) ) if $jmline =~ /nnxnp/;
    #
    #push( @lmll_alleles, join( "\t", "$markers[0]\_$markers[1]", $jmline ) ) if $jmline =~ /lmxll/;

}

### new header after removing mother anf father
#my @new_header;
#for ( my $i = 4; $i < scalar @header; $i++ ) {
#    push( @new_header, $header[$i] ) if ( $header[$i] ne $mname && $header[$i] ne $fname && $header[$i] ne "9202" );
#}
#
#print "\nNew Headers:\n", join ":", @new_header, "\n";
#
#open OUT, ">$ARGV[0].pval$alpha.loc" or die "Unable to open output file.";
#print OUT "name=$ARGV[0].loc";
#print OUT "\npopt=$pop_type";
#print OUT "\nnloc=", scalar @all_alleles;
#print OUT "\nnind=", scalar @new_header;
#print OUT join "\n", "\n", @all_alleles, "\n";
#print OUT join "\n", "\n", "individual names:", @new_header;
#
#### print lmxll
#open LMLL, ">lmxll.$mname.$ARGV[0].pval$alpha.loc" or die "Unable to open output file.";
#print LMLL "name=$ARGV[0].loc";
#print LMLL "\npopt=$pop_type";
#print LMLL "\nnloc=", scalar @lmll_alleles;
#print LMLL "\nnind=", scalar @new_header;
#print LMLL join "\n", "\n", @lmll_alleles, "\n";
#print LMLL join "\n", "\n", "individual names:", @new_header;
#
#close LMLL;
#
#### print nnxnp
#open NNNP, ">nnxnp.$fname.$ARGV[0].pval$alpha.loc" or die "Unable to open output file.";
#print NNNP "name=$ARGV[0].loc";
#print NNNP "\npopt=$pop_type";
#print NNNP "\nnloc=", scalar @nnnp_alleles;
#print NNNP "\nnind=", scalar @new_header;
#print NNNP join "\n", "\n", @nnnp_alleles, "\n";
#print NNNP join "\n", "\n", "individual names:", @new_header;
#
#close NNNP;

#######################################################################

sub allele_to_joinmap {
    my $refarr_header  = shift;
    my $refarr_alleles = shift;
    my $mother_name    = shift;
    my $father_name    = shift;
    my %markers;
    my $jmap_line;
    for ( my $i = 0; $i < scalar @$refarr_alleles; $i++ ) { $markers{ $$refarr_header[$i] } = $$refarr_alleles[$i]; }

    my ( $refHash_genotypes, $mar_type ) = allele_type( $markers{$mother_name}, $markers{$father_name} );

    $jmap_line .= "$mar_type";
    for ( my $i = 0; $i < scalar @$refarr_alleles; $i++ ) {
        next if ( $$refarr_header[$i] eq $mname || $$refarr_header[$i] eq $fname || $$refarr_header[$i] eq "9202" );
        $jmap_line .= "\t";
        $jmap_line .= $$refHash_genotypes{ allele_to_digital( $$refarr_alleles[$i] ) };
    }

    return ($jmap_line);

}

sub dna_to_digital {
    my $allele   = shift;
    my $maternal = shift;
    my $paternal = shift;

}

sub allele_type_lmxll {
    my $maternal = shift;
    my $paternal = shift;
    my ( $mal, $pal, $f1pattern, $marker_type, %f2_genotypes );

    ( $mal, $pal ) = allele_to_digital( $maternal, $paternal );
    $f1pattern = "$mal $pal";

    if ( $f1pattern eq "0/1 0/0" || $f1pattern eq "1/0 0/0" || $f1pattern eq "0/1 1/1" || $f1pattern eq "1/0 1/1" ) {
        $marker_type         = "<lmxll>";
        $f2_genotypes{"0/0"} = "ll";
        $f2_genotypes{"0/1"} = "lm";
        $f2_genotypes{"1/0"} = "lm";
        $f2_genotypes{"1/1"} = "ll";

    }
    elsif ( $f1pattern eq "0/0 0/1" || $f1pattern eq "0/0 1/0" || $f1pattern eq "1/1 0/1" || $f1pattern eq "1/1 1/0" ) {
        $marker_type         = "<nnxnp>";
        $f2_genotypes{"0/0"} = "nn";
        $f2_genotypes{"0/1"} = "np";
        $f2_genotypes{"1/0"} = "np";
        $f2_genotypes{"1/1"} = "nn";

    }
    elsif ( $f1pattern eq "0/1 0/1" || $f1pattern eq "1/0 1/0" || $f1pattern eq "1/0 0/1" || $f1pattern eq "0/1 1/0" ) {
        $marker_type         = "<hkxhk>";
        $f2_genotypes{"0/0"} = "hh";
        $f2_genotypes{"0/1"} = "hk";
        $f2_genotypes{"1/0"} = "hk";
        $f2_genotypes{"1/1"} = "kk";

    }
    else {
        $marker_type         = "<XXxXX>";
        $f2_genotypes{"0/0"} = "xx";
        $f2_genotypes{"0/1"} = "xx";
        $f2_genotypes{"1/0"} = "xx";
        $f2_genotypes{"1/1"} = "xx";

        #print STDERR "\nSkipping line due to non-informative marker type: $f1pattern\n";# if $verbose;
    }

    return ( \%f2_genotypes, $marker_type );
}



sub allele_type_ahb {
    my $maternal = shift;
    my $paternal = shift;
    my ( $mal, $pal, $f1pattern, $marker_type, %f2_genotypes );

    ( $mal, $pal ) = allele_to_digital( $maternal, $paternal );
    $f1pattern = "$mal $pal";

    if ( $f1pattern eq "0/1 0/0" || $f1pattern eq "1/0 0/0" || $f1pattern eq "0/1 1/1" || $f1pattern eq "1/0 1/1" ) {
        $marker_type         = "<lmxll>";
        $f2_genotypes{"0/0"} = "ll";
        $f2_genotypes{"0/1"} = "lm";
        $f2_genotypes{"1/0"} = "lm";
        $f2_genotypes{"1/1"} = "ll";

    }
    elsif ( $f1pattern eq "0/0 0/1" || $f1pattern eq "0/0 1/0" || $f1pattern eq "1/1 0/1" || $f1pattern eq "1/1 1/0" ) {
        $marker_type         = "<nnxnp>";
        $f2_genotypes{"0/0"} = "nn";
        $f2_genotypes{"0/1"} = "np";
        $f2_genotypes{"1/0"} = "np";
        $f2_genotypes{"1/1"} = "nn";

    }
    elsif ( $f1pattern eq "0/1 0/1" || $f1pattern eq "1/0 1/0" || $f1pattern eq "1/0 0/1" || $f1pattern eq "0/1 1/0" ) {
        $marker_type         = "<hkxhk>";
        $f2_genotypes{"0/0"} = "hh";
        $f2_genotypes{"0/1"} = "hk";
        $f2_genotypes{"1/0"} = "hk";
        $f2_genotypes{"1/1"} = "kk";

    }
    else {
        $marker_type         = "<XXxXX>";
        $f2_genotypes{"0/0"} = "xx";
        $f2_genotypes{"0/1"} = "xx";
        $f2_genotypes{"1/0"} = "xx";
        $f2_genotypes{"1/1"} = "xx";

        #print STDERR "\nSkipping line due to non-informative marker type: $f1pattern\n";# if $verbose;
    }

    return ( \%f2_genotypes, $marker_type );
}




sub allele_to_digital {
    my $maternal = shift;
    my $paternal = shift;
    my ( $mal, $pal, $f1pattern, $marker_type, %f2_genotypes );

    ### if alleles are not in 1/1 format. convert it to 0 1 format
    if ( !$paternal ) {
        if    ( length($maternal) == 1 ) { $mal = "0/0" }
        elsif ( length($maternal) == 2 ) { $mal = "0/1" }
        else                             { $mal = "x/x" }
        return $mal;
    }
    elsif ( $maternal eq $paternal ) {

        if    ( length($maternal) == 1 ) { $mal = "0/0" }
        elsif ( length($maternal) == 2 ) { $mal = "0/1" }
        else                             { $mal = "x/x" }
        if    ( length($paternal) == 1 ) { $pal = "0/0" }
        elsif ( length($paternal) == 2 ) { $pal = "0/1" }
        else                             { $pal = "x/x" }

    }
    else {
        if    ( length($maternal) == 1 && length($paternal) == 1 ) { $mal = "0/0"; $pal = "1/1" }
        elsif ( length($maternal) == 1 && length($paternal) == 2 ) { $mal = "0/0"; $pal = "0/1" }
        elsif ( length($maternal) == 2 && length($paternal) == 1 ) { $mal = "0/1"; $pal = "1/1" }
        elsif ( length($maternal) == 2 && length($paternal) == 2 ) { $pal = "0/0"; $pal = "1/1" }
        else                                                       { $mal = "x/x"; $pal = "x/x" }
    }

    return ( $mal, $pal );

}

sub estimate_segregation {
    my $hMarkers = shift;  ### ref of hash with name as key and allele as value.
    my $aName    = shift;  ### ref of array of sample names
    my $aAlleles = shift;  ### ref of array of alleles in same order as names in $aNames
    my $mname    = shift;  ### mother name
    my $fname    = shift;  ### father name
    my $f1_name  = shift;  ### f1 progeny name.
    
    my %segregant;
    
    foreach my$name(keys %{$hMarkers}){
      next if ($name eq $mname || $name eq $fname || $name eq $f1_name || $$hMarkers{$name} eq '-');
       $segregant{$$hMarkers{$name}}++;
    }
    
    
    my($segval,$seg_al,$seg_rat,$alrat);
    my@alrat;
    foreach my$allel(sort{length($a) <=> length($b)} keys %segregant){
      
      $segval.="$allel:$segregant{$allel};"; ### AT:61
      $seg_al.="$allel:";                    ### AT:GC;
      $seg_rat.="$segregant{$allel}:";       ### 1:3
      $alrat.="$segregant{$allel}\t";
      push(@alrat,$segregant{$allel});
      
      
    }
    
    if (scalar@alrat ==1) {$alrat=join("\t",0,@alrat)}
    elsif(scalar@alrat ==2){$alrat=join("\t",@alrat)}
    else{$alrat="\t"}
    
    return ($segval,$seg_al,$seg_rat,$alrat);
    

}

sub non_redundant_allele{
  my$allele=join("",@_);
  my%nonred;

  my$all=non_redun($allele);
  return sort_order($all);
   
   
}

sub non_redun{
  my$string=shift;
  my%nonred;
  my$all;
  my@chars=split("",$string);
  foreach(@chars){$nonred{$_}=1;}
  foreach(keys %nonred){$all.=$_}
  return ($all);
  
}

sub sort_order{
  my$allele=shift;
  return join("", sort{$a cmp $b} split //,$allele);
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

    my $chisq2 =
      ( ( $val1 - $rat1 * $indtotal )**2 ) / ( $rat1 * $indtotal + 0.00000000000000000000001 ) +
      ( ( $val2 - $rat2 * $indtotal )**2 ) / ( $rat2 * $indtotal + 0.00000000000000000000001 );

    $chisq2 =
      ( ( abs( $val1 - $rat1 * $indtotal ) - 0.5 )**2 ) / ( $rat1 * $indtotal + 0.00000000000000000000001 ) +
      ( ( abs( $val2 - $rat2 * $indtotal ) - 0.5 )**2 ) / ( $rat2 * $indtotal + 0.00000000000000000000001 )
      if $yatecor;

    return ($chisq2);

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
    my $yatecor       = shift;

    my $chisq = 0;
    my $pval  = 1;
    my $total = sum(@$ref_array_obs);

    # create ratios from expected values and convert to expected numbers.
    my $totalexp = sum(@$ref_array_rat);
    my @ratios   = map { sprintf( "%.2f", $_ / $totalexp ) } @$ref_array_rat;
    my $degree   = scalar @$ref_array_obs - 1;

    for ( my $i = 0; $i < scalar @$ref_array_obs; $i++ ) {
        $chisq += ( ( $$ref_array_obs[$i] - $ratios[$i] * $total )**2 ) / ( $ratios[$i] + 0.00000000000000000000001 ) if !$yatecor;
        $chisq += ( abs( $$ref_array_obs[$i] - $ratios[$i] * $total - 0.5 )**2 ) / ( $ratios[$i] + 0.00000000000000000000001 ) if $yatecor;
    }
    
    
    
    
    return ( $chisq, $degree, chisqrprob( $degree, $chisq ) );
}
