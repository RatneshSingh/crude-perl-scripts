use strict;
use Getopt::Long;

our($diff,$fold,$cu,$ct,$su,$st,$min,$model,$verbose,$exp_model,$ratios,$out);

my $uasge="

$0 -f cuffnorm_fpkm_file

options:
-fold|d  fold difference to filter for. Use - for downregulated genes.
-cu           column containing control untreated samples [1].
-ct           column containing control treated samples [2].
-su           column containing sample untreated samples[3].
-st           column containing sample treated samples[4].
-min|m        exclude gene if sum expression raw values dont add up to this number.
-model|exp|e  show genes only which fits this model. e.g '1 2 1 1'
-ratios       Print log2 fold change for all the combination of samples.
-out|o        Output file for saving result. Default is STDOUT
-verbose|v    Print the progress of run.

";





my$options=GetOptions(
  "diff|f=s"=>\$diff,
  "fold|d=i" =>\$fold,
  "control_untreated|cu=i"=>\$cu,
  "control_treated|ct=i"=>\$ct,
  "sample_untreated|su=i"=>\$su,
  "sample_treated|st=i"=>\$st,
  "min_value|min|m=i"=>\$min,
  "model|exp|e=s"=>\$model,
  "ratios|r"=>\$ratios,
  "out|o=s"=>\$out,
  "verbose|v"=>\$verbose

);

## default columns for samples and controls
$cu=$cu?$cu:1;
$ct=$ct?$ct:2;
$su=$su?$su:3;
$st=$st?$st:4;
$min=$min?$min:4;
#$model=$model?$model:"1\t1\t1\t2";
open DIFF, "$diff" or die "unable to open $diff\n\n$uasge";
#$fold=$fold?$fold:2;

if ($out) {
   open(OUTPUT, '>', $out) or die;
} else {
   *OUTPUT = *STDOUT;
}


while (<DIFF>) {
    chomp;
    next if /^\s*$/;
    my $line = $_;
    $line =~ s/^\s+//g;
    my@header;
    if ( $. == 1 ) {


      @header = split( /\s+/, $line );
      print "Comparing gene expression between\nControl_untreated:$header[$cu]\nControl_treated:$header[$ct]\nSample_untreated:$header[$su]\nSample_untreated:$header[$st]\n" if !$ratios;
     my $head= join("\t",@header[0,$cu,$ct,$su,$st]);

      if ($ratios) {
        print "Calculating Log2 fold change between all the combination of samples \n";
          my @result;

          foreach my $a (1..$#header) {
            foreach my $b ($a+1..$#header) {
              push @result, $header[$a]."->".$header[$b];
            }
          }
          $head=join("\t","gene",@result);
      }
      print OUTPUT "$head\n";
      next;
    }

    my @line = split( /\s+/, $line );
    next if ($line[$cu]+$line[$ct]+$line[$su]+$line[$st])<$min;

    if ($ratios) {
          my @result;

          foreach my $a (1..$#line) {
            foreach my $b ($a+1..$#line) {
              push @result, round(fold_diff($line[$b],$line[$a],2),2);
            }
          }
          print OUTPUT join("\t","$line[0]",@result),"\n";
    }elsif ($model) {

        print OUTPUT join("\t",@line[0,$cu,$ct,$su,$st],"\n") if match_expression("$model","$line[$cu]\t$line[$ct]\t$line[$su]\t$line[$st]");
    }elsif($fold < 0 ){
        print OUTPUT join("\t",@line[0,$cu,$ct,$su,$st],"\n") if exp_diff($line[$cu],$line[$ct],$line[$su],$line[$st],2) <= logN($fold,2);
    }elsif($fold>0){
        print OUTPUT join("\t",@line[0,$cu,$ct,$su,$st],"\n") if abs(exp_diff($line[$cu],$line[$ct],$line[$su],$line[$st],2)) >= logN($fold,2);
    }else{
      print OUTPUT join("\t",@line[0,$cu,$ct,$su,$st],"\n") if abs(fold_diff($line[$cu],$line[$ct],2)) > abs(logN($fold,2))|| abs(fold_diff($line[$su],$line[$st],2) > abs(logN($fold,2)));
    }

}



############################################
sub logN{
  my$val=shift;
  my$base=shift;
  #print "\nLog of $val with base $base:",log($val)/log($base);
  return (log($val)/log($base)) if $val > 0;
  return (undef) if $val < 0;

}

sub round{
    my $number = shift;
	my $decimals=shift;
	$decimals=$decimals?$decimals:3;
    substr( $number + ( '0.' . '0' x $decimals . '5' ), 0, $decimals + length(int($number)) + 1 );
}


sub pvalue {
    my $val1 = shift;
    my $val2 = shift;
    return 1 if ($val1+$val2) < 1;
    ## calculate simple statistics ( a-b/sqrt(a+b) ) and pvalue ((0.0098*($zval**4))-(0.1303*($zval**3))+(0.6536*($zval**2))-(1.472*$zval)+1.2601;)
    my $zval = ( $val1 > $val2 ? $val1 - $val2 : $val2 - $val1 ) / sqrt( $val1 + $val2 );
    my $pval = ( 0.0098 * ( $zval**4 ) ) - ( 0.1303 * ( $zval**3 ) ) + ( 0.6536 * ( $zval**2 ) ) - ( 1.472 * $zval ) + 1.2601;
    return $pval;
}

#sub get_ratio{
#    my $control_untreated = shift;
#    my $control_treated = shift;
#    my $sample_untreated = shift;
#    my $sample_treated = shift;
#    my $base = shift;
#
#    $control_untreated = $control_untreated>1?$control_untreated:1;
#    $control_treated = $control_treated>1?$control_treated:1;
#    $sample_untreated = $sample_untreated>1?$sample_untreated:1;
#    $sample_treated = $sample_treated>1?$sample_treated:1;
#
#
#    return (( $sample_treated/$sample_untreated ),( $control_treated/$control_untreated ));
#
#    ## calculate fold difference
#    if ($base) {
#
#      return logN($fold_diff,$base);
#
#    }
#    else{
#
#      return ( $fold_diff );
#
#    }
#
#}



sub fold_diff {
    my $val1 = shift;
    my $val2 = shift;
    my$base=shift;
    ## calculate fold difference
    return logN(( ( $val1 > 1 ? $val1 : 1 ) / ( $val2 > 1 ? $val2 : 1 ) ),$base) if $base;
    return ( ( $val1 > 1 ? $val1 : 1 ) / ( $val2 > 1 ? $val2 : 1 ) );

}


sub exp_diff {
    my $control_untreated = shift;
    my $control_treated = shift;
    my $sample_untreated = shift;
    my $sample_treated = shift;
    my $base = shift;

    $control_untreated = $control_untreated>1?$control_untreated:1;
    $control_treated = $control_treated>1?$control_treated:1;
    $sample_untreated = $sample_untreated>1?$sample_untreated:1;
    $sample_treated = $sample_treated>1?$sample_treated:1;


    my$fold_diff=( $sample_treated/$sample_untreated )/( $control_treated/$control_untreated );

    ## calculate fold difference
    if ($base) {

      return logN($fold_diff,$base);

    }
    else{

      return ( $fold_diff );

    }
}

sub match_expression{
  my$model = shift;
  my$expr=shift;
  my$pval=shift;
  #my$cu =shift;
  #my$ct=shift;
  #my$su=shift;
  #my$st=shift;

  $pval=$pval?$pval:0.05;

  my@model=split(/\s+/,$model);
  my@expr=split(/\s+/,$expr);

  ## test if model fits the expression
  return 1 if (chisq(\@expr,\@model,'pval') >=$pval);

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

    my$ref_array_exp=[];
    my $chisq=0;
    my $pval=1;
    my $total=sum(@$ref_array_obs);

    # create ratios from expected values to make it usefull with exp values, ratios with decimal or number
    my $totalexp=sum(@$ref_array_rat);
    my @ratios=map{sprintf("%.5f",$_/$totalexp)} @$ref_array_rat;
    my $degree= scalar @$ref_array_obs -1;

    for(my$i=0;$i<scalar@$ref_array_obs;$i++){
      push(@$ref_array_exp,$ratios[$i]*$total);
       $chisq += ( ( $$ref_array_obs[$i] - $ratios[$i]*$total )**2 ) / ( $ratios[$i]+0.00000000000000000000001) if !$yatecor;
       $chisq += ( abs( $$ref_array_obs[$i] - $ratios[$i]*$total - 0.5)**2 ) / ( $ratios[$i]+0.00000000000000000000001) if $yatecor;
    }
    print "\n***\tObs:@$ref_array_obs\tExp:@$ref_array_exp\tRatios:@$ref_array_rat\tchisq:$chisq\tprob:",sprintf("%.6f",chisqrprob($degree,$chisq)),"\tdegree:$degree"  if $verbose;
    if ($return=~/chisq/i) {return $chisq}
    elsif($return=~/pval/i){return sprintf("%.6f",chisqrprob($degree,$chisq))}
    elsif($return=~/deg/i){return ($degree)}
    elsif($return==1){return ($chisq)}
    elsif($return==2){return ($chisq,sprintf("%.6f",chisqrprob($degree,$chisq)))}
    elsif($return==3){return ($chisq,sprintf("%.6f",chisqrprob($degree,$chisq,$degree)))}
    else{return ($chisq)}
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
