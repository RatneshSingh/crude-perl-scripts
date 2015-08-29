#!/usr/bin/perl
use warnings;
use GD::Graph::pie;
use GD::Graph::bars;
use GD::Graph::lines;
use GD::Graph::area;

use Chart::Clicker;
use Chart::Clicker::Data::Series;
use Chart::Clicker::Data::DataSet;

use Getopt::Long;
use strict;

our ( %GO, %total_lines, %plot_names, %plot_values, %image, $help, $file, $opt_print_graphics, $opt_save_annot );

my $opt_GOnum      = 40;
my $opt_width      = 900;
my $opt_height     = 900;
my $opt_graph_type = 'pie';
my $font_type      = 'arial';
my $font_size      = 50;

my $result = GetOptions(
  "input|i=s"      => \$file,
  "show_top|sto=i" => \$opt_GOnum,   # plot this many in graph from sorted list.
  "width|wd=i"     => \$opt_width,
  "height|ht=i"    => \$opt_height,
  "graph_type|gt=s" => \$opt_graph_type,
  "font_type|ft=s"  => \$font_type,
  "font_size|fs=i"  => \$font_size,
  "print_graphics|pg"=>\$opt_print_graphics,
  "b2g_annot|b2g"=>\$opt_save_annot, # save for Blast2Go as annotation file
  "help|h"          => \$help
);

my $usage = "
perl script options...

  input|i     file. Output file from IPRscan in raw format.
  show_top|sto    How many data to plot from sorted list [10].
  width|wd       Graph width [800]
  height|ht      Graph height [800]
  graph_type|gt  Graph type to show e.g. bar, line, area, pie [pie]
  font_type|ft   Font type for text. e.g. 'verdana', 'arial', gdMediumBoldFont. ['ariel']
  font_size|fs   Font size for text [50]

Output options
  print_graphics|pg Save graphics for results.
  b2g_annot|b2g Save annotation for Blast2Go program.
  help|h      help

";
die $usage if ($help||!$file);

open ANNOT,">$file.annot" if $opt_save_annot;
#'C:\Users\ratnesh.singh\Documents\My Projects\Assembled_Pineapple_Sequences\iprscan_runs\test_F153_run61_mira_assembly_asPairedSolexa_out.padded.IPRSCAN.raw';
open IPRSCAN, "$file";
my$seq_num=0;
while (<IPRSCAN>) {
  next if /^\s*$/;
  my @line = split /\t/, $_;
  foreach (@line) { s/^\s+//g; s/\s+$//g; }
  $seq_num++;
  #    print $line[13]."\n" if $line[13];
  if ( $line[13] ) {

    foreach( $line[13] =~ /Molecular Function:([\w\,\s]+)\((GO:\d+)\)/ ) {

      #print $1."\t".$2."\n";
      $GO{'Molecular Function'}{$2}{'count'} = 0
        if !$GO{'Molecular Function'}{$2}{'count'};
      $GO{'Molecular Function'}{$2}{'GO'} = $1;
      $GO{'Molecular Function'}{$2}{'count'}++;
      $total_lines{'Molecular Function'}++;

      if($opt_save_annot){print ANNOT"Seq$seq_num\t$2\t$1\n"}

    }
    foreach ( $line[13] =~ /Biological Process:([\w\,\s]+)\((GO:\d+)\)/ ) {

      #print $1."\t".$2."\n";
      $GO{'Biological Process'}{$2}{'count'} = 0
        if !$GO{'Biological Process'}{$2}{'count'};
      $GO{'Biological Process'}{$2}{'GO'} = $1;
      $GO{'Biological Process'}{$2}{'count'}++;
      $total_lines{'Biological Process'}++;
     if($opt_save_annot){print ANNOT"Seq$seq_num\t$2\t$1\n"}

    }
    foreach( $line[13] =~ /Cellular Component:([\w\,\s]+)\((GO:\d+)\)/ ) {

      #print $1."\t".$2."\n";
      $GO{'Cellular Component'}{$2}{'count'} = 0
        if !$GO{'Cellular Component'}{$2}{'count'};
      $GO{'Cellular Component'}{$2}{'GO'} = $1;
      $GO{'Cellular Component'}{$2}{'count'}++;
      $total_lines{'Cellular Component'}++;
      if($opt_save_annot){print ANNOT"Seq$seq_num\t$2\t$1\n"}
    }
    $total_lines{'total'}++;
  }

}

open OUT, ">$file.summary.table";
print OUT "Total lines procssed: $total_lines{'total'}\n";

foreach my $category ( keys %GO ) {
  next if $category =~ /^\s*$/;

  @{ $plot_names{$category} } =
    sort { $GO{$category}{$b}{'count'} <=> $GO{$category}{$a}{'count'} }
    keys %{ $GO{$category} };

  print OUT "\n\n\n$category\t$total_lines{$category}\n";

  #  print "\nselected names: @{ $plot_names{$category} }";

  foreach my $go_anno ( keys %{ $GO{$category} } ) {

    print OUT "$GO{$category}{$go_anno}{'GO'}\t"
      . $GO{$category}{$go_anno}{'count'} . "\n";

    #push( @{ $plot_names{$category} },  $GO{$category}{$go_anno}{'GO'} );
    #push( @{ $plot_values{$category} }, $GO{$category}{$go_anno}{'count'} );
  }
}

# sort and collect top 20 hits

#############################################
# Draw graph for the given data

# set data in array format. one array for name and other for value. both array should have same number of elements
if($opt_print_graphics){
foreach my $category ( keys %plot_names ) {
  next if $category =~ /^\s*$/;
  my @plot_name;
  my @plot_value;
  for ( my $i = 0 ; $i <= $opt_GOnum ; $i++ ) {

    push( @plot_name, $GO{$category}{ ${ $plot_names{$category} }[$i] }{'GO'} );
    push( @plot_value,
      $GO{$category}{ ${ $plot_names{$category} }[$i] }{'count'} );
  }
  my @data = ( \@plot_name, \@plot_value );
  my $mygraph;


##################################################
# print "plotting following data\n".join ("\t",@plot_name)."\n".join ("\t",@plot_value)."\n";
  ###########################################################################
  # for drawing Pie chart
  if ( lc $opt_graph_type eq 'pie' ) {
    $mygraph = GD::Graph::pie->new( max( $opt_width, $opt_height ),
      max( $opt_width, $opt_height ) );
    $mygraph->set(
      title => "Pie chart showing top $opt_GOnum GO annotations for $category",
      '3d'  => 0,
      start_angle => 90,
      label => $category,
      suppress_angle => 5, # suppress the slices less than this angle from getting the label
      axislabelclr => 'black'
    ) or warn $mygraph->error;

    $mygraph->set_label_font( $font_type, $font_size );    # set label font
    $mygraph->set_value_font( $font_type, $font_size );    # set value font
    $mygraph->set_title_font( $font_type, $font_size );

  }
  ###########################################################################
  # for drawing Area graph
  elsif ( lc $opt_graph_type eq 'area' ) {
    $mygraph = GD::Graph::area->new( $opt_width, $opt_height );
    $mygraph->set(
      x_label => $category,
      y_label => 'Number of elements',
      title   => "Top $opt_GOnum GO annotations for $category",
    ) or warn $mygraph->error;
  }
  ###########################################################################
  # for drawing line graph
  elsif ( lc $opt_graph_type eq 'lines' ) {
    $mygraph = GD::Graph::lines->new( $opt_width, $opt_height );
    $mygraph->set(
      x_label => $category,
      y_label => 'Number of elements',
      title   => "Top $opt_GOnum GO annotations for $category",

      # Draw datasets in 'solid', 'dashed' and 'dotted-dashed' lines
      line_types => [ 1, 2, 4 ],

      # Set the thickness of line
      line_width => 2,

      # Set colors for datasets
      dclrs => [
        'white',   'lgray',  'gray',    'dgray',   'black',  'lblue',
        'blue',    'dblue',  'gold',    'lyellow', 'yellow', 'dyellow',
        'lgreen',  'green',  'dgreen',  'lred',    'red',    'dred',
        'lpurple', 'purple', 'dpurple', 'lorange', 'orange', 'pink',
        'dpink',   'marine', 'cyan',    'lbrown',  'dbrown'
      ],
    ) or warn $mygraph->error;
    $mygraph->set_x_label_font($font_type,$font_size);
    $mygraph->set_y_label_font($font_type,$font_size);
    $mygraph->set_x_axis_font($font_type,$font_size);
    $mygraph->set_y_axis_font($font_type,$font_size);
    $mygraph->set_values_font($font_type,$font_size);
    $mygraph->set_legend_font($font_type,$font_size);
    $mygraph->set_legend(@plot_name);
    $mygraph->set_legend_font($font_type,$font_size);
  }
###########################################################################
  # for drawing bar graph
  elsif ( lc $opt_graph_type eq 'bar' ) {

    $mygraph = GD::Graph::bars->new( $opt_width, $opt_height );
    $mygraph->set(
      x_label => $category,
      y_label => 'Number of elements',
      title   => "Top $opt_GOnum GO annotations for $category",
    ) or warn $mygraph->error;

    $mygraph->set_x_label_font($font_type,$font_size);
    $mygraph->set_y_label_font($font_type,$font_size);
    $mygraph->set_x_axis_font($font_type,$font_size);
    $mygraph->set_y_axis_font($font_type,$font_size);
    $mygraph->set_values_font($font_type,$font_size);
    $mygraph->set_legend_font($font_type,$font_size);
    $mygraph->set_legend(@plot_name);
    $mygraph->set_legend_font($font_type,$font_size);


  }
#####################################################################
  # plot the data on selected type of graph
  my $myimage = $mygraph->plot( \@data ) or die $mygraph->error;

  #open OUT, ">$file.$category.jpg";
  #print OUT $myimage->jpeg;
  #close OUT;

  open OUT, ">$file.$category.png";
  print OUT $myimage->png;
  close OUT;
}
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

############################################################################################
# Draw chart from provided data
############################################################################################
sub draw_chart{

    my$Ref_array_keys=shift;
    my$Ref_array_values=shift;
  # build the chart
  my $chart = Chart::Clicker->new;

  # build the series (static here, will usually be supplied arrayrefs from elsewhere)
  my $series = Chart::Clicker::Data::Series->new(
      keys    =>      \$Ref_array_keys,
      values  =>      \$Ref_array_values,
  );

  # build the dataset
  my $dataset = Chart::Clicker::Data::DataSet->new(
      series  =>      [ $series ],
  );

  # add the dataset to the chart
  $chart->add_to_datasets($dataset);

  # write the chart to a file
  $chart->write_output('chart.png');

}
