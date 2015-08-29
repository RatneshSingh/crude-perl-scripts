#!/usr/local/bin/perl -w
# Change above line to point to your perl binary
use GD::Graph::bars;
use GD::Graph::lines;
use GD::Graph::xypoints;
use GD::Graph::xylines;
use GD::Simple;
use List::Util qw(sum);
use strict;
use Getopt::Long;

our($map,$fst,$out_file,$plot_column,$step,$line,$window);
$plot_column=18;

my$result=GetOptions(
  "map=s"=>\$map,
  "fst=s"=>\$fst,
  "out=s"=>\$out_file,
  "window=i"=>\$window,
  "step=i"=>\$step,
  "plot=i"=>\$plot_column,
  "line"=>\$line,
);

my$usage="
$0 -options...

-map  Blast mapping file produced my linkage mapping script.
-fst  Pairwise Fst file produced by Stacks.
-out  Output file to save results [modified fst filename]
-window Window size for estimating Fst in sliding window.
-step step size of window sliding along the chromosome.
-plot Plot what fst to plot [18].
      8:Fst
      14:Corrected Fst
      17:Amova Fst
      18:Corrected Amova Fst
-line Draw line connecting points.
      Works well for average Fst with -step options.
      Creates zigzag for non-averaged Fsts.
";

die "\n$usage\n" if (!$map || !$fst);

if (defined $plot_column && $plot_column!=8 && $plot_column!=14 && $plot_column!=17 && $plot_column!=18){
print "\nOnly 8,14,17,18 are the options acepted for -plot options. Setting -plot to 18\n";
$plot_column=18;
}

#$line=1 if $step;
open(MAP, "$map") or die "Unable to find BlastMapping file:$map\n$usage\n";
open(FST, "$fst") or die "Unable to open FST file :$fst\nUsage:\n$usage\n";

my$tempout=$fst;
$tempout=~s/.*\///ig;
$tempout=~s/.tsv/.WidChrLocs.tsv/gi;

$out_file=$out_file?$out_file:$tempout;
my%chrinfo;
my%blast;
while (<MAP>) {
  s/^\s+|\s+$//g;
  my@line=split /\t/;
  next if $line[0]=~/^\s*$/;
  next if $line[6]=~/^\s*$/;
  $blast{$line[0]}{'chr'}=$line[6];
  $blast{$line[0]}{'loc'}=$line[7];
  #print "\n$line[6]\t$blast{$line[0]}{'chr'}";
  $chrinfo{$line[6]}{'chrlen'} = $chrinfo{$line[6]}{'chrlen'} ? $chrinfo{$line[6]}{'chrlen'} : 1;
  $chrinfo{$line[6]}{'chrlen'} = $line[7] if  $chrinfo{$line[6]}{'chrlen'} < $line[7];

}

my%fst;
my%data_x;
my%data_y;
my%average;
open(TEXT, ">$out_file");
print TEXT join "\t",qw(Marker Chromosome Location Fst CorrectedFst AmovaFst CorrectedAmovaFst) if !$step;
while (<FST>) {
  s/^\s+|\s+$//g;
  my@line=split /\t+|\s+/;
  next if $line[1]=~/^\s*$/;
  #print "\n$blast{$line[0]}{'chr'}";
  next if !defined $blast{$line[1]}{'chr'};
  next if $blast{$line[1]}{'chr'}=~/chrun|chrsy/i;
  my$loc=$blast{$line[1]}{'loc'};
  my$chr=$blast{$line[1]}{'chr'};
  my$fst=$line[$plot_column];
  print TEXT "\n";
  print TEXT join "\t",$line[1],$chr,$loc,@line[8,14,17,18];


  ####collecting for average
  #my$wind=(int($blast{$line[1]}{'loc'}/$step)+1) if $step;

  ## collect fst in bins of overlapping coords.
  my$cur_bin=abs((int(($loc-$window)/$step))-1);
  my$steps_per_win=$window/$step;
  my$bin_start=$cur_bin*$step+$window-$step;
  $bin_start = $bin_start > 0?$bin_start:0;
  my$bin_end=$loc+$window+1;
  $bin_end= $bin_end < $chrinfo{$chr}{'chrlen'} ? $bin_end : $chrinfo{$chr}{'chrlen'};
  #print "\nbinstart\t$bin_start .. $bin_end. \t for Loc\t$loc";

  for (my$i=$bin_start;$i<=$bin_end;$i+=$step){

    push(@{$average{$chr}{$i}},$fst) if ($i > $loc && $i < $loc + $window ) ;;
    #print "\n\t$loc is added to group bin:$i ..",$i+$step if ($i > $loc && $i < $loc + $window ) ;

  }

  #$average{$chr}{$wind}=() if !$average{$chr}{$wind};
  #push(@{$average{$chr}{$wind}},$line[$plot_column]) if $step;




### create dataset for plotting if average is not to be calculated
  if (!$step) {
    push( @{  $data_x{$chr}  }, $blast{$line[1]}{'loc'});
    push( @{  $data_y{$chr}  }, $line[$plot_column]);
  }
}

###calculate average if window was provided
if ($step) {
  open(OUT2,">$out_file.Fst.Smoothed$step.bp.tsv");
  print OUT2 "Chromosome\tWindowLocation\tAverage_Fst";
  foreach my$chr(keys %average){
    @{  $data_y{$chr}  }=();
    @{  $data_x{$chr}  }=();
    foreach my$window(sort {$a<=> $b} keys %{$average{$chr}}){
      my$mean_fst=mean( @{$average{$chr}{$window}} );
      next if $mean_fst <= 0;
      #print "\nChr $chr\tWindow $window\tMean fst $mean_fst";
      push( @{  $data_y{$chr}  },  $mean_fst  );
      push( @{  $data_x{$chr}  },  $window    );
      print OUT2 "\n";
      print OUT2 join "\t",$chr,$window,$mean_fst;


    }
  }
}



### plot the data
###plotData
my$fwidth=1200;
my$fheight=300;

my$gap=50;

my$Wtiled=$fwidth;
my$Htiled=$fheight * scalar (keys %data_x) + $gap * scalar(keys %data_x);
### composite images coords

my $myimage=new GD::Image($Wtiled, $Htiled);
my$white = $myimage->colorAllocate(255,255,255);
$myimage->filledRectangle(0,0,$Wtiled,$Htiled,$white);
my%mygraph;
my$count=0;

my$font_file="/usr/share/fonts/liberation/LiberationSans-Regular.ttf";
foreach my$chr(sort{only_nums($a,100)<=> only_nums($b,100) } keys %data_x){

    my $graph = GD::Graph::xypoints->new($fwidth, $fheight);
     $graph = GD::Graph::xylines->new($fwidth, $fheight) if $line;

     #$graph->set( 'y_number_format' => \&y_format );
    $graph->set_title_font($font_file,18);
    $graph->set_x_label_font($font_file,12);
    $graph->set_y_label_font($font_file,12);
    $graph->set_x_axis_font($font_file,10);
    $graph->set_y_axis_font($font_file,10);
    $graph->set_values_font($font_file,10);
    $graph->set( y_max_value => 1,
                y_min_value => 0,
                'y_number_format' => \&y_format,
                'x_number_format' => \&x_format,
                 marker_size => 3,
                 markers     => [7],
                      ) or die $graph->error;
    #$graph->set( x_max_value => $plot_xmax ) or die $graph->error if defined $plot_xmax;
    #$graph->set( x_min_value => $plot_xmin ) or die $graph->error if defined $plot_xmin;
    #$graph->set( y_min_value => $plot_ymin ) or die $graph->error if defined $plot_ymin;




    $graph->set(
        x_label     => $chr,
        y_label     => 'Fst',
        title       => $step?"Average Fst value in Sliding window of $window bp, sliding by $step bp along the Chromosomes $chr":"Fst value along the Chromosomes $chr",
    ) or warn $graph->error;

    my@data=(\@{$data_x{$chr}},\@{$data_y{$chr}})  ;
    my$gd=$graph->plot(\@data) or die $graph->error;

    my$dstX=0;
    my$dstY=$count*$fheight+$count*$gap;
    $myimage->GD::Image::copy($gd,$dstX,$dstY,0,0,$fwidth,$fheight);

    $count++;

}








# Open a file for writing
open(PICTURE, ">$out_file.png") or die("Cannot open file for writing image");

# Make sure we are writing to a binary stream
binmode PICTURE;

# Convert the image to PNG and print it to the file PICTURE
print PICTURE $myimage->png;
close PICTURE;





sub mean {
    return sum(@_)/@_;
}


sub only_nums{ my $string=shift;
    my $return_if_non_num=shift;
    $string=~s/\D+//g; # expand all numbers to 4 digits
   return $return_if_non_num if length$string<1;
   return $string;
}


sub round{
    my $number = shift;
	my $decimals=shift;
	$decimals=$decimals?$decimals:3;
    substr( $number + ( '0.' . '0' x $decimals . '5' ), 0, $decimals + length(int($number)) + 1 );
}


sub y_format
    {
        my $value = shift;
        my $ret;

        $ret=round($value,2);

        return $ret;
    }

sub x_format
    {
        my $value = shift;
        my $ret;

        $ret=round($value/1000000,2);

        return "$ret Mb";
    }
